#include "config.h"
#include "globals.h"
#include "ads/trie/trie.h"
#include "ads/calc_utils.h"
#include "ads/sax/sax.h"
#include "ads/sfa/dft.h"
#include "ads/spartan/spartan.h"
#include "ads/pisa/pisa.h"
#include "ads/inmemory_index_engine.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <pthread.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct symbolic_trie_node {
    struct symbolic_trie_node *children[8];
    struct symbolic_trie_node *parent;
    sax_type *min_word;
    sax_type *max_word;
    sax_type **words;
    file_position_type *positions;
    int size;
    int capacity;
    int split_dimension;
    uint64_t used_dimensions;
    unsigned char leaf;
    unsigned char mbb_valid;
    unsigned char split_exhausted;
} symbolic_trie_node;

struct symbolic_trie_index {
    symbolic_trie_node *root;
    int dimensions;
    int bound_dimensions;
    unsigned long split_exhausted_leaves;
    unsigned long node_count;
};

/* A full entropy pass over every dimension of a very large root is
 * prohibitively expensive (64 complete scans at the current maximum word
 * length).  A deterministic, evenly-spaced sample retains the distribution
 * information needed for split selection while keeping the expensive full
 * pass to the one partition pass that cannot be avoided. */
#define TRIE_ROOT_SPLIT_SAMPLE_SIZE 1000000L

/* The index settings describe the full split word (normally 64 dimensions).
 * Lower bounds intentionally see only the configured bound prefix.  A local
 * shallow copy is safe for concurrent queries and reuses the representation-
 * specific SFA/SPARTAN/PISA/SAX bound implementations unchanged. */
static float trie_lower_bound(const struct symbolic_trie_index *trie, isax_index *index,
                              const ts_type *transform, sax_type *sax_min,
                              sax_type *sax_max, int dimensions, float bsf) {
    isax_index shadow_index = *index;
    isax_index_settings shadow_settings = *index->settings;
    shadow_settings.n_segments = dimensions;
    shadow_index.settings = &shadow_settings;
    return messi_minidist_range_raw(&shadow_index, (float *) transform, sax_min, sax_max,
                                    shadow_settings.max_sax_cardinalities, bsf);
}

/* Leaf records have one concrete word, unlike an MBR range.  Going through
 * the individual-word dispatcher enables its existing 16-dimensional SIMD
 * paths for SAX, SFA/PISA, and SPARTAN. */
static float trie_record_lower_bound(const struct symbolic_trie_index *trie,
                                     isax_index *index, const ts_type *transform,
                                     sax_type *word, float bsf) {
    isax_index shadow_index = *index;
    isax_index_settings shadow_settings = *index->settings;
    shadow_settings.n_segments = trie->bound_dimensions;
    shadow_index.settings = &shadow_settings;
    return messi_minidist_raw(&shadow_index, (float *) transform, word,
                               shadow_settings.max_sax_cardinalities, bsf);
}

typedef struct {
    unsigned long checked_nodes;
    unsigned long lower_bounds;
    unsigned long exact_distances;
    unsigned long long total_microseconds;
    float approximate_distance;
} trie_query_stats;

typedef struct {
    float lower_bound;
    int record_index;
} trie_leaf_candidate;

/* One scratch heap per batch-query worker.  It is retained across leaves and
 * avoids per-record allocation in best-first leaf refinement. */
typedef struct {
    trie_leaf_candidate *candidates;
    int capacity;
} trie_query_scratch;

/* FFTW plan creation is not thread-safe on all supported builds. */
static pthread_mutex_t trie_fftw_plan_lock = PTHREAD_MUTEX_INITIALIZER;

static symbolic_trie_node *trie_node_create(int dimensions, symbolic_trie_node *parent,
                                            uint64_t used_dimensions) {
    symbolic_trie_node *node = calloc(1, sizeof(*node));
    if (node == NULL) return NULL;
    node->min_word = malloc((size_t) dimensions);
    node->max_word = malloc((size_t) dimensions);
    if (node->min_word == NULL || node->max_word == NULL) {
        free(node->min_word); free(node->max_word); free(node); return NULL;
    }
    node->parent = parent;
    node->used_dimensions = used_dimensions;
    node->split_dimension = -1;
    node->leaf = 1;
    return node;
}

static void trie_node_destroy(symbolic_trie_node *node) {
    if (node == NULL) return;
    for (int i = 0; i < 8; ++i) trie_node_destroy(node->children[i]);
    for (int i = 0; i < node->size; ++i) free(node->words[i]);
    free(node->words); free(node->positions); free(node->min_word); free(node->max_word); free(node);
}

static unsigned long trie_node_count(const symbolic_trie_node *node) {
    if (node == NULL) return 0;
    unsigned long count = 1;
    for (int i = 0; i < 8; ++i) count += trie_node_count(node->children[i]);
    return count;
}

static unsigned long long trie_monotonic_microseconds(void) {
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    return (unsigned long long) time.tv_sec * 1000000ULL + (unsigned long long) time.tv_nsec / 1000ULL;
}

static void trie_print_query_stats(int query_index, const struct symbolic_trie_index *trie,
                                   const trie_query_stats *stats, float distance) {
    printf("%3d: %10d %10d %6lu %13lu %14d %12d %12d %20.3f %12.3f %12llu\n",
           query_index, 0, 0, trie->node_count, stats->checked_nodes, 0, 0, 0,
           stats->approximate_distance, distance, stats->total_microseconds);
}

static void trie_save_query_stats(const struct symbolic_trie_index *trie,
                                  const trie_query_stats *stats, float distance) {
    total_tree_nodes = (int) trie->node_count;
    checked_nodes = (int) stats->checked_nodes;
    BYTES_ACCESSED = 0;
    loaded_nodes = 0;
    loaded_records = 0;
    APPROXIMATE = stats->approximate_distance;
    total_querying_time = (double) stats->total_microseconds;
    total_init_time = 0.0;
    total_tree_pass_time = 0.0;
    TOTAL_PQ_INSERT_TIME = 0;
    TOTAL_PQ_REMOVE_TIME = 0;
    TOTAL_LB_DIST_CALC_TIME = 0;
    TOTAL_REAL_DIST_CALC_TIME = 0;
    total_time = (double) stats->total_microseconds;
    LBDcalculationnumber = stats->lower_bounds;
    RDcalculationnumber = stats->exact_distances;
    SAVE_STATS(distance)
}

static void trie_node_update_mbb(symbolic_trie_node *node, const sax_type *word, int dimensions) {
    if (!node->mbb_valid) {
        memcpy(node->min_word, word, (size_t) dimensions);
        memcpy(node->max_word, word, (size_t) dimensions);
        node->mbb_valid = 1;
        return;
    }
    for (int d = 0; d < dimensions; ++d) {
        if (word[d] < node->min_word[d]) node->min_word[d] = word[d];
        if (word[d] > node->max_word[d]) node->max_word[d] = word[d];
    }
}

static int trie_scratch_reserve(trie_query_scratch *scratch, int capacity) {
    if (scratch == NULL) return 0;
    if (capacity <= scratch->capacity) return 1;
    trie_leaf_candidate *candidates = realloc(scratch->candidates,
                                              sizeof(*candidates) * (size_t) capacity);
    if (candidates == NULL) return 0;
    scratch->candidates = candidates;
    scratch->capacity = capacity;
    return 1;
}

static void trie_heap_sift_down(trie_leaf_candidate *heap, int size, int root) {
    while (1) {
        int child = root * 2 + 1;
        if (child >= size) return;
        if (child + 1 < size && heap[child + 1].lower_bound < heap[child].lower_bound) ++child;
        if (heap[root].lower_bound <= heap[child].lower_bound) return;
        trie_leaf_candidate temporary = heap[root];
        heap[root] = heap[child];
        heap[child] = temporary;
        root = child;
    }
}

static void trie_heap_build(trie_leaf_candidate *heap, int size) {
    for (int i = size / 2 - 1; i >= 0; --i) trie_heap_sift_down(heap, size, i);
}

static trie_leaf_candidate trie_heap_pop(trie_leaf_candidate *heap, int *size) {
    trie_leaf_candidate result = heap[0];
    --*size;
    if (*size > 0) {
        heap[0] = heap[*size];
        trie_heap_sift_down(heap, *size, 0);
    }
    return result;
}

static int trie_leaf_append(symbolic_trie_node *node, sax_type *word, file_position_type position,
                            int dimensions) {
    if (node->size == node->capacity) {
        int capacity = node->capacity == 0 ? 32 : node->capacity * 2;
        sax_type **words = realloc(node->words, sizeof(*words) * (size_t) capacity);
        if (words == NULL) return 0;
        node->words = words;
        file_position_type *positions = realloc(node->positions, sizeof(*positions) * (size_t) capacity);
        if (positions == NULL) return 0;
        node->positions = positions; node->capacity = capacity;
    }
    trie_node_update_mbb(node, word, dimensions);
    node->words[node->size] = word;
    node->positions[node->size++] = position;
    return 1;
}

static int trie_choose_split(const symbolic_trie_node *node, int dimensions,
                             const double *variances) {
    int best = -1;
    double best_score = -1.0;
    for (int d = 0; d < dimensions; ++d) {
        if (node->used_dimensions & (UINT64_C(1) << d)) continue;
        int counts[8] = {0};
        for (int i = 0; i < node->size; ++i) ++counts[node->words[i][d] >> 5];
        double entropy = 0.0;
        int occupied = 0;
        for (int b = 0; b < 8; ++b) if (counts[b] != 0) {
            double p = (double) counts[b] / node->size;
            entropy -= p * log(p); ++occupied;
        }
        if (occupied < 2) continue;
        double score = entropy / log(8.0);
        if (variances != NULL) score *= variances[d];
        if (score > best_score) { best_score = score; best = d; }
    }
    return best;
}

static int trie_split_leaf(struct symbolic_trie_index *trie, isax_index *index, symbolic_trie_node *node) {
    int dimension = trie_choose_split(node, trie->dimensions, index->settings->symbolic_variances);
    if (dimension < 0) { node->split_exhausted = 1; ++trie->split_exhausted_leaves; return 0; }
    uint64_t used = node->used_dimensions | (UINT64_C(1) << dimension);
    symbolic_trie_node *children[8] = {0};
    for (int i = 0; i < node->size; ++i) {
        int bucket = node->words[i][dimension] >> 5;
        if (children[bucket] == NULL && (children[bucket] = trie_node_create(trie->dimensions, node, used)) == NULL) return 0;
        if (!trie_leaf_append(children[bucket], node->words[i], node->positions[i], trie->dimensions)) return 0;
        node->words[i] = NULL;
    }
    free(node->words); free(node->positions);
    node->words = NULL; node->positions = NULL; node->size = node->capacity = 0;
    node->leaf = 0; node->split_dimension = dimension;
    memcpy(node->children, children, sizeof(children));
    return 1;
}

/* Return 1 after a successful split, 0 if no usable dimension remains, and
 * -1 on allocation failure.  This is intentionally used only for the root:
 * after it has been partitioned, ordinary task-parallel subtree construction
 * handles the much smaller leaves with the exact entropy routine above. */
static int trie_split_root_sampled_parallel(struct symbolic_trie_index *trie,
                                            isax_index *index,
                                            symbolic_trie_node *node) {
    const long sample_size = node->size < TRIE_ROOT_SPLIT_SAMPLE_SIZE
                                 ? node->size : TRIE_ROOT_SPLIT_SAMPLE_SIZE;
    const int dimensions = trie->dimensions;
    const int workers = maxquerythread > 0 ? maxquerythread : 1;
    unsigned long *histograms = calloc((size_t) workers * dimensions * 8,
                                       sizeof(*histograms));
    if (histograms == NULL) return -1;

    const double selection_start = messi_monotonic_seconds();
#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
    {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        unsigned long *local = histograms + (size_t) worker * dimensions * 8;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (long sample = 0; sample < sample_size; ++sample) {
            long record = (long) (((unsigned long long) sample * node->size) / sample_size);
            const sax_type *word = node->words[record];
            for (int d = 0; d < dimensions; ++d) ++local[d * 8 + (word[d] >> 5)];
        }
    }

    int dimension = -1;
    double best_score = -1.0;
    for (int d = 0; d < dimensions; ++d) {
        if (node->used_dimensions & (UINT64_C(1) << d)) continue;
        unsigned long counts[8] = {0};
        for (int worker = 0; worker < workers; ++worker) {
            const unsigned long *local = histograms + (size_t) worker * dimensions * 8 + d * 8;
            for (int bucket = 0; bucket < 8; ++bucket) counts[bucket] += local[bucket];
        }
        double entropy = 0.0;
        int occupied = 0;
        for (int bucket = 0; bucket < 8; ++bucket) if (counts[bucket] != 0) {
            double p = (double) counts[bucket] / sample_size;
            entropy -= p * log(p);
            ++occupied;
        }
        if (occupied < 2) continue;
        double score = entropy / log(8.0);
        if (index->settings->symbolic_variances != NULL) score *= index->settings->symbolic_variances[d];
        if (score > best_score) { best_score = score; dimension = d; }
    }
    free(histograms);
    const double selection_end = messi_monotonic_seconds();
    if (dimension < 0) {
        node->split_exhausted = 1;
        ++trie->split_exhausted_leaves;
        return 0;
    }

    unsigned long *counts = calloc((size_t) workers * 8, sizeof(*counts));
    unsigned long *offsets = calloc((size_t) workers * 8, sizeof(*offsets));
    sax_type *local_min = malloc((size_t) workers * 8 * dimensions);
    sax_type *local_max = calloc((size_t) workers * 8 * dimensions, sizeof(*local_max));
    if (counts == NULL || offsets == NULL || local_min == NULL || local_max == NULL) {
        free(counts); free(offsets); free(local_min); free(local_max);
        return -1;
    }
    memset(local_min, UCHAR_MAX, (size_t) workers * 8 * dimensions);

#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
    {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        unsigned long *local = counts + (size_t) worker * 8;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int i = 0; i < node->size; ++i) ++local[node->words[i][dimension] >> 5];
    }

    unsigned long bucket_sizes[8] = {0};
    for (int bucket = 0; bucket < 8; ++bucket) {
        unsigned long offset = 0;
        for (int worker = 0; worker < workers; ++worker) {
            offsets[(size_t) worker * 8 + bucket] = offset;
            offset += counts[(size_t) worker * 8 + bucket];
        }
        bucket_sizes[bucket] = offset;
    }

    symbolic_trie_node *children[8] = {0};
    const uint64_t used = node->used_dimensions | (UINT64_C(1) << dimension);
    for (int bucket = 0; bucket < 8; ++bucket) if (bucket_sizes[bucket] != 0) {
        children[bucket] = trie_node_create(dimensions, node, used);
        if (children[bucket] == NULL ||
            (children[bucket]->words = malloc(sizeof(*children[bucket]->words) * bucket_sizes[bucket])) == NULL ||
            (children[bucket]->positions = malloc(sizeof(*children[bucket]->positions) * bucket_sizes[bucket])) == NULL) {
            for (int b = 0; b < 8; ++b) if (children[b] != NULL) {
                free(children[b]->words); free(children[b]->positions);
                free(children[b]->min_word); free(children[b]->max_word); free(children[b]);
            }
            free(counts); free(offsets); free(local_min); free(local_max);
            return -1;
        }
        children[bucket]->capacity = children[bucket]->size = (int) bucket_sizes[bucket];
    }

#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
    {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        unsigned long *local_offsets = offsets + (size_t) worker * 8;
        sax_type *min_word = local_min + (size_t) worker * 8 * dimensions;
        sax_type *max_word = local_max + (size_t) worker * 8 * dimensions;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int i = 0; i < node->size; ++i) {
            sax_type *word = node->words[i];
            int bucket = word[dimension] >> 5;
            unsigned long position = local_offsets[bucket]++;
            children[bucket]->words[position] = word;
            children[bucket]->positions[position] = node->positions[i];
            sax_type *min = min_word + bucket * dimensions;
            sax_type *max = max_word + bucket * dimensions;
            for (int d = 0; d < dimensions; ++d) {
                if (word[d] < min[d]) min[d] = word[d];
                if (word[d] > max[d]) max[d] = word[d];
            }
        }
    }

    for (int bucket = 0; bucket < 8; ++bucket) if (children[bucket] != NULL) {
        for (int d = 0; d < dimensions; ++d) {
            sax_type min = UCHAR_MAX, max = 0;
            for (int worker = 0; worker < workers; ++worker) {
                const sax_type local_low = local_min[((size_t) worker * 8 + bucket) * dimensions + d];
                const sax_type local_high = local_max[((size_t) worker * 8 + bucket) * dimensions + d];
                if (local_low < min) min = local_low;
                if (local_high > max) max = local_high;
            }
            children[bucket]->min_word[d] = min;
            children[bucket]->max_word[d] = max;
        }
        children[bucket]->mbb_valid = 1;
    }
    free(counts); free(offsets); free(local_min); free(local_max);
    free(node->words); free(node->positions);
    node->words = NULL; node->positions = NULL; node->size = node->capacity = 0;
    node->leaf = 0; node->split_dimension = dimension;
    memcpy(node->children, children, sizeof(children));
    const double partition_end = messi_monotonic_seconds();
    fprintf(stderr,
            ">>> trie root split timing: sample=%ld selection=%.3fs partition=%.3fs dimension=%d\n",
            sample_size, selection_end - selection_start, partition_end - selection_end, dimension);
    return 1;
}

static enum response trie_word_from_ts(isax_index *index, const ts_type *ts, sax_type *word,
                                       ts_type *transform, fftw_workspace *fftw) {
    int d = index->settings->n_segments;
    if (index->settings->function_type == 3) {
        if (paa_from_ts((ts_type *) ts, transform, index->settings) != SUCCESS) return FAILURE;
        sax_from_paa(transform, word, index->settings);
    } else if (index->settings->function_type == 4) {
        memcpy(fftw->ts, ts, sizeof(ts_type) * index->settings->timeseries_size);
        if (fft_from_ts(index, d, index->settings->n_coefficients != 0, fftw) != SUCCESS) return FAILURE;
        memcpy(transform, fftw->transform, sizeof(ts_type) * d);
        sfa_from_fft(index, transform, word);
    } else if (index->settings->function_type == 5) {
        if (pca_from_ts(index, ts, transform) != SUCCESS) return FAILURE;
        spartan_from_pca(index, transform, word);
    } else if (index->settings->function_type == 6) {
        if (pisa_pca_from_ts(index, ts, transform, fftw) != SUCCESS) return FAILURE;
        sfa_from_fft(index, transform, word);
    } else return FAILURE;
    return SUCCESS;
}

static int trie_insert(struct symbolic_trie_index *trie, isax_index *index, sax_type *word,
                       file_position_type position) {
    symbolic_trie_node *node = trie->root;
    while (!node->leaf) {
        trie_node_update_mbb(node, word, trie->dimensions);
        int bucket = word[node->split_dimension] >> 5;
        if (node->children[bucket] == NULL) {
            node->children[bucket] = trie_node_create(trie->dimensions, node, node->used_dimensions);
            if (node->children[bucket] == NULL) return 0;
        }
        node = node->children[bucket];
    }
    if (!trie_leaf_append(node, word, position, trie->dimensions)) return 0;
    if (node->size > index->settings->max_leaf_size && !node->split_exhausted) trie_split_leaf(trie, index, node);
    return 1;
}

/* Each invocation owns its node and may therefore split it without locking.
 * The only shared state is the OpenMP task scheduler. */
static void trie_build_subtree(struct symbolic_trie_index *trie, isax_index *index,
                               symbolic_trie_node *node, int depth) {
    if (!node->leaf || node->size <= index->settings->max_leaf_size || node->split_exhausted) return;
    if (!trie_split_leaf(trie, index, node)) return;
    for (int i = 0; i < 8; ++i) if (node->children[i] != NULL) {
#ifdef _OPENMP
#pragma omp task firstprivate(i) if (depth < 6)
#endif
        trie_build_subtree(trie, index, node->children[i], depth + 1);
    }
}

enum response symbolic_trie_build(isax_index *index, const char *path, long ts_num,
                                  int filetype_int, int apply_znorm) {
    if (index == NULL || index->settings == NULL || ts_num <= 0 || filetype_int) return FAILURE;
    double build_start = messi_monotonic_seconds();
    FILE *file = fopen(path, "rb");
    if (file == NULL) return FAILURE;
    rawfile = malloc((size_t) ts_num * index->settings->timeseries_size * sizeof(*rawfile));
    if (rawfile == NULL || fread(rawfile, sizeof(*rawfile),
        (size_t) ts_num * index->settings->timeseries_size, file) !=
        (size_t) ts_num * index->settings->timeseries_size) { fclose(file); free(rawfile); rawfile = NULL; return FAILURE; }
    fclose(file);
    double read_end = messi_monotonic_seconds();
    struct symbolic_trie_index *trie = calloc(1, sizeof(*trie));
    if (trie == NULL || (trie->root = trie_node_create(index->settings->n_segments, NULL, 0)) == NULL) {
        free(trie); free(rawfile); rawfile = NULL; return FAILURE;
    }
    trie->dimensions = index->settings->n_segments;
    trie->bound_dimensions = index->settings->trie_bound_dimensions > 0
                                 ? index->settings->trie_bound_dimensions
                                 : trie->dimensions;
    /* Materialize all full words in parallel.  Subsequent splits move only
     * these word pointers and positions, never raw time series. */
    trie->root->words = calloc((size_t) ts_num, sizeof(*trie->root->words));
    trie->root->positions = malloc((size_t) ts_num * sizeof(*trie->root->positions));
    trie->root->capacity = trie->root->size = (int) ts_num;
    int failed = trie->root->words == NULL || trie->root->positions == NULL;
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#endif
    {
        fftw_workspace fftw = {0};
        ts_type *transform = malloc(sizeof(*transform) * (size_t) trie->dimensions);
        pthread_mutex_lock(&trie_fftw_plan_lock);
        fftw_workspace_init(&fftw, index->settings->timeseries_size);
        pthread_mutex_unlock(&trie_fftw_plan_lock);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (long i = 0; i < ts_num; ++i) {
            if (failed) continue;
            ts_type *ts = rawfile + (size_t) i * index->settings->timeseries_size;
            if (apply_znorm) znorm(ts, index->settings->timeseries_size);
            sax_type *word = malloc((size_t) trie->dimensions);
            if (word == NULL || transform == NULL || trie_word_from_ts(index, ts, word, transform, &fftw) != SUCCESS) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                failed = 1;
                free(word);
            } else {
                trie->root->words[i] = word;
                trie->root->positions[i] = (file_position_type) ((size_t) i * index->settings->timeseries_size);
            }
        }
        free(transform);
        pthread_mutex_lock(&trie_fftw_plan_lock);
        fftw_workspace_destroy(&fftw);
        pthread_mutex_unlock(&trie_fftw_plan_lock);
    }
    if (failed) { trie_node_destroy(trie->root); free(trie); free(rawfile); rawfile = NULL; return FAILURE; }
    double transform_end = messi_monotonic_seconds();
    for (long i = 0; i < ts_num; ++i) trie_node_update_mbb(trie->root, trie->root->words[i], trie->dimensions);
    double root_mbb_end = messi_monotonic_seconds();

    /* Split the initially enormous root before entering the task region.  The
     * sampled splitter uses OpenMP worksharing itself; nested worksharing from
     * inside the later omp single region would otherwise be serialized. */
    if (trie->root->size > index->settings->max_leaf_size) {
        int root_split = trie_split_root_sampled_parallel(trie, index, trie->root);
        if (root_split < 0) {
            trie_node_destroy(trie->root); free(trie); free(rawfile); rawfile = NULL;
            return FAILURE;
        }
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#pragma omp single
#endif
    {
        if (trie->root->leaf) {
            trie_build_subtree(trie, index, trie->root, 0);
        } else {
            for (int i = 0; i < 8; ++i) if (trie->root->children[i] != NULL) {
#ifdef _OPENMP
#pragma omp task firstprivate(i)
#endif
                trie_build_subtree(trie, index, trie->root->children[i], 1);
            }
        }
    }
    trie->node_count = trie_node_count(trie->root);
    double split_end = messi_monotonic_seconds();
    index->trie = trie;
    index->total_records = ts_num;
    fprintf(stderr,
            ">>> trie build timing: read=%.3fs transform=%.3fs root-mbb=%.3fs "
            "split=%.3fs total=%.3fs\n",
            read_end - build_start, transform_end - read_end,
            root_mbb_end - transform_end, split_end - root_mbb_end,
            split_end - build_start);
    return SUCCESS;
}

static float trie_scan_leaf_best_first(isax_index *index, const symbolic_trie_node *node,
                                       const ts_type *query, const ts_type *transform, float bsf,
                                       trie_query_stats *stats, trie_query_scratch *scratch) {
    if (!trie_scratch_reserve(scratch, node->size)) {
        /* Allocation failure is not a correctness failure: retain the former
         * streaming order rather than dropping candidates. */
        for (int i = 0; i < node->size; ++i) {
            if (stats != NULL) ++stats->lower_bounds;
            else __sync_fetch_and_add(&LBDcalculationnumber, 1);
            if (trie_record_lower_bound(index->trie, index, transform, node->words[i], bsf) < bsf) {
                if (stats != NULL) ++stats->exact_distances;
                else __sync_fetch_and_add(&RDcalculationnumber, 1);
                float d = ts_ed((ts_type *) query, rawfile + node->positions[i],
                                index->settings->timeseries_size, bsf);
                if (d < bsf) bsf = d;
            }
        }
        return bsf;
    }

    int candidate_count = 0;
    for (int i = 0; i < node->size; ++i) {
        if (stats != NULL) ++stats->lower_bounds;
        else __sync_fetch_and_add(&LBDcalculationnumber, 1);
        float lower_bound = trie_record_lower_bound(index->trie, index, transform,
                                                     node->words[i], bsf);
        if (lower_bound < bsf) {
            scratch->candidates[candidate_count].lower_bound = lower_bound;
            scratch->candidates[candidate_count++].record_index = i;
        }
    }

    trie_heap_build(scratch->candidates, candidate_count);
    while (candidate_count > 0 && scratch->candidates[0].lower_bound < bsf) {
        trie_leaf_candidate candidate = trie_heap_pop(scratch->candidates, &candidate_count);
        if (stats != NULL) ++stats->exact_distances;
        else __sync_fetch_and_add(&RDcalculationnumber, 1);
        float d = ts_ed((ts_type *) query,
                        rawfile + node->positions[candidate.record_index],
                        index->settings->timeseries_size, bsf);
        if (d < bsf) bsf = d;
    }
    return bsf;
}

static float trie_search_node(isax_index *index, const symbolic_trie_node *node, const ts_type *query,
                              const ts_type *transform, float bsf, trie_query_stats *stats,
                              const symbolic_trie_node *skip_leaf, trie_query_scratch *scratch) {
    if (node == skip_leaf) return bsf;
    if (stats != NULL) ++stats->checked_nodes;
    float lower = trie_lower_bound(index->trie, index, transform, node->min_word, node->max_word,
                                   index->trie->dimensions, bsf);
    if (lower >= bsf) return bsf;
    if (!node->leaf) {
        typedef struct { const symbolic_trie_node *node; float bound; } trie_child_bound;
        trie_child_bound ordered[8];
        int n = 0;
        for (int i = 0; i < 8; ++i) if (node->children[i] != NULL) {
            ordered[n].node = node->children[i];
            ordered[n++].bound = trie_lower_bound(index->trie, index, transform,
                node->children[i]->min_word, node->children[i]->max_word,
                index->trie->dimensions, bsf);
        }
        for (int i = 1; i < n; ++i) {
            trie_child_bound current = ordered[i]; int j = i - 1;
            while (j >= 0 && ordered[j].bound > current.bound) { ordered[j + 1] = ordered[j]; --j; }
            ordered[j + 1] = current;
        }
        for (int i = 0; i < n; ++i) if (ordered[i].bound < bsf)
            bsf = trie_search_node(index, ordered[i].node, query, transform, bsf, stats, skip_leaf, scratch);
        return bsf;
    }
    return trie_scan_leaf_best_first(index, node, query, transform, bsf, stats, scratch);
}

static const symbolic_trie_node *trie_seed_leaf(isax_index *index, const sax_type *query_word,
                                                 const ts_type *transform) {
    const symbolic_trie_node *node = index->trie->root;
    while (!node->leaf) {
        int bucket = query_word[node->split_dimension] >> 5;
        const symbolic_trie_node *next = node->children[bucket];
        if (next == NULL) {
            float best_bound = FLT_MAX;
            for (int i = 0; i < 8; ++i) if (node->children[i] != NULL) {
                float bound = trie_lower_bound(index->trie, index, transform,
                                                node->children[i]->min_word,
                                                node->children[i]->max_word,
                                                index->trie->dimensions, best_bound);
                if (bound < best_bound) {
                    best_bound = bound;
                    next = node->children[i];
                }
            }
        }
        if (next == NULL) return NULL;
        node = next;
    }
    return node;
}

static void trie_search_task(isax_index *index, const symbolic_trie_node *node, const ts_type *query,
                             const ts_type *transform, float *shared_bsf, pthread_mutex_t *bsf_lock,
                             int depth, const symbolic_trie_node *skip_leaf) {
    if (node == skip_leaf) return;
    __sync_fetch_and_add(&checked_nodes, 1);
    float bsf;
    pthread_mutex_lock(bsf_lock); bsf = *shared_bsf; pthread_mutex_unlock(bsf_lock);
    if (trie_lower_bound(index->trie, index, transform, node->min_word, node->max_word,
                         index->trie->dimensions, bsf) >= bsf) return;
    if (!node->leaf) {
        for (int i = 0; i < 8; ++i) if (node->children[i] != NULL) {
#ifdef _OPENMP
#pragma omp task firstprivate(i) if (depth < 6)
#endif
            trie_search_task(index, node->children[i], query, transform, shared_bsf, bsf_lock, depth + 1,
                             skip_leaf);
        }
        return;
    }
    for (int i = 0; i < node->size; ++i) {
        pthread_mutex_lock(bsf_lock); bsf = *shared_bsf; pthread_mutex_unlock(bsf_lock);
        __sync_fetch_and_add(&LBDcalculationnumber, 1);
        if (trie_record_lower_bound(index->trie, index, transform, node->words[i], bsf) < bsf) {
            __sync_fetch_and_add(&RDcalculationnumber, 1);
            float d = ts_ed((ts_type *) query, rawfile + node->positions[i], index->settings->timeseries_size, bsf);
            if (d < bsf) { pthread_mutex_lock(bsf_lock); if (d < *shared_bsf) *shared_bsf = d; pthread_mutex_unlock(bsf_lock); }
        }
    }
}

static float trie_parallel_exact_search(isax_index *index, const ts_type *query,
                                        const ts_type *transform, float bsf,
                                        const symbolic_trie_node *skip_leaf) {
    pthread_mutex_t bsf_lock = PTHREAD_MUTEX_INITIALIZER;
    float shared_bsf = bsf;
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#pragma omp single
#endif
    trie_search_task(index, index->trie->root, query, transform, &shared_bsf, &bsf_lock, 0, skip_leaf);
    pthread_mutex_destroy(&bsf_lock);
    return shared_bsf;
}

enum response symbolic_trie_query_file(isax_index *index, const char *path, int query_count,
                                       int filetype_int, int apply_znorm, float minimum_distance) {
    if (index == NULL || path == NULL || query_count < 0 || filetype_int) return FAILURE;
    FILE *file = fopen(path, "rb");
    if (file == NULL) return FAILURE;
    int ts_length = index->settings->timeseries_size;
    ts_type *query = malloc(sizeof(*query) * (size_t) ts_length);
    ts_type *transform = malloc(sizeof(*transform) * (size_t) index->settings->n_segments);
    sax_type *word = malloc((size_t) index->settings->n_segments);
    trie_query_scratch scratch = {0};
    fftw_workspace fftw = {0}; fftw_workspace_init(&fftw, ts_length);
    if (query == NULL || transform == NULL || word == NULL) { fclose(file); free(query); free(transform); free(word); fftw_workspace_destroy(&fftw); return FAILURE; }
    for (int i = 0; i < query_count; ++i) {
        if (fread(query, sizeof(*query), (size_t) ts_length, file) != (size_t) ts_length ||
            (apply_znorm && (znorm(query, ts_length), 0)) ||
            trie_word_from_ts(index, query, word, transform, &fftw) != SUCCESS) {
            fclose(file); free(query); free(transform); free(word); fftw_workspace_destroy(&fftw); return FAILURE;
        }
        /* The transform drives lower bounds; the word selects the seed path. */
        trie_query_stats stats = {0};
        unsigned long long query_start = trie_monotonic_microseconds();
        const symbolic_trie_node *seed_leaf = trie_seed_leaf(index, word, transform);
        float bsf = seed_leaf == NULL ? FLT_MAX :
                trie_search_node(index, seed_leaf, query, transform, FLT_MAX, &stats, NULL, &scratch);
        stats.approximate_distance = bsf;
        float distance;
        if (maxquerythread > 1) {
            LBDcalculationnumber = RDcalculationnumber = 0;
            checked_nodes = 0;
            if (minimum_distance < bsf) bsf = minimum_distance;
            distance = trie_parallel_exact_search(index, query, transform, bsf, seed_leaf);
            stats.checked_nodes += checked_nodes;
            stats.lower_bounds += LBDcalculationnumber;
            stats.exact_distances += RDcalculationnumber;
            LBDcalculationnumber = RDcalculationnumber = 0;
            checked_nodes = 0;
        } else {
            distance = trie_search_node(index, index->trie->root, query, transform,
                                        minimum_distance < bsf ? minimum_distance : bsf, &stats, seed_leaf, &scratch);
        }
        stats.total_microseconds = trie_monotonic_microseconds() - query_start;
        trie_save_query_stats(index->trie, &stats, distance);
        if (i == 0) PRINT_STATS_HEADER();
        trie_print_query_stats(i, index->trie, &stats, distance);
    }
    fclose(file); free(query); free(transform); free(word); free(scratch.candidates);
    fftw_workspace_destroy(&fftw); return SUCCESS;
}

enum response symbolic_trie_query_file_batch(isax_index *index, const char *path, int query_count,
                                             int filetype_int, int apply_znorm, float minimum_distance) {
    if (index == NULL || path == NULL || query_count < 0 || filetype_int) return FAILURE;
    FILE *file = fopen(path, "rb");
    int length = index->settings->timeseries_size, dimensions = index->settings->n_segments;
    ts_type *queries = malloc((size_t) query_count * length * sizeof(*queries));
    float *distances = malloc((size_t) query_count * sizeof(*distances));
    trie_query_stats *stats = calloc((size_t) query_count, sizeof(*stats));
    if (file == NULL || queries == NULL || distances == NULL || stats == NULL ||
        fread(queries, sizeof(*queries), (size_t) query_count * length, file) != (size_t) query_count * length) {
        if (file) fclose(file); free(queries); free(distances); free(stats); return FAILURE;
    }
    fclose(file);
    int failed = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#endif
    {
        fftw_workspace fftw = {0};
        ts_type *transform = malloc((size_t) dimensions * sizeof(*transform));
        sax_type *word = malloc((size_t) dimensions);
        trie_query_scratch scratch = {0};
        pthread_mutex_lock(&trie_fftw_plan_lock);
        fftw_workspace_init(&fftw, length);
        pthread_mutex_unlock(&trie_fftw_plan_lock);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int i = 0; i < query_count; ++i) {
            ts_type *query = queries + (size_t) i * length;
            if (apply_znorm) znorm(query, length);
            if (transform == NULL || word == NULL || trie_word_from_ts(index, query, word, transform, &fftw) != SUCCESS) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                failed = 1;
            } else {
                /* One worker owns one query: avoid nested per-query tasks. */
                unsigned long long query_start = trie_monotonic_microseconds();
                const symbolic_trie_node *seed_leaf = trie_seed_leaf(index, word, transform);
                float bsf = seed_leaf == NULL ? FLT_MAX :
                        trie_search_node(index, seed_leaf, query, transform, FLT_MAX, &stats[i], NULL, &scratch);
                stats[i].approximate_distance = bsf;
                if (minimum_distance < bsf) bsf = minimum_distance;
                distances[i] = trie_search_node(index, index->trie->root, query, transform, bsf,
                                                &stats[i], seed_leaf, &scratch);
                stats[i].total_microseconds = trie_monotonic_microseconds() - query_start;
            }
        }
        free(transform); free(word); free(scratch.candidates);
        pthread_mutex_lock(&trie_fftw_plan_lock);
        fftw_workspace_destroy(&fftw);
        pthread_mutex_unlock(&trie_fftw_plan_lock);
    }
    if (!failed) {
        if (query_count > 0) PRINT_STATS_HEADER();
        for (int i = 0; i < query_count; ++i) {
            trie_save_query_stats(index->trie, &stats[i], distances[i]);
            trie_print_query_stats(i, index->trie, &stats[i], distances[i]);
        }
        fflush(stdout);
    }
    free(queries); free(distances); free(stats); return failed ? FAILURE : SUCCESS;
}

query_result symbolic_trie_exact_search(isax_index *index, const ts_type *query,
                                        const ts_type *transform, float bsf) {
    query_result result = { FLT_MAX, NULL, 0 };
    if (index == NULL || index->trie == NULL) return result;
    if (maxquerythread <= 1) {
        trie_query_scratch scratch = {0};
        result.distance = trie_search_node(index, index->trie->root, query, transform, bsf,
                                           NULL, NULL, &scratch);
        free(scratch.candidates);
    }
    else result.distance = trie_parallel_exact_search(index, query, transform, bsf, NULL);
    return result;
}

void symbolic_trie_destroy(isax_index *index) {
    if (index != NULL && index->trie != NULL) {
        trie_node_destroy(index->trie->root);
        free(index->trie); index->trie = NULL;
    }
}
