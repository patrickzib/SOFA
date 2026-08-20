#include "config.h"
#include "globals.h"
#include "ads/trie/trie.h"
#include "ads/calc_utils.h"
#include "ads/sfa/dft.h"
#include "ads/spartan/spartan.h"
#include "ads/pisa/pisa.h"
#include "ads/inmemory_index_engine.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
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
    unsigned long split_exhausted_leaves;
    unsigned long node_count;
};

typedef struct {
    unsigned long checked_nodes;
    unsigned long lower_bounds;
    unsigned long exact_distances;
    unsigned long long total_microseconds;
    float approximate_distance;
} trie_query_stats;

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

static enum response trie_word_from_ts(isax_index *index, const ts_type *ts, sax_type *word,
                                       ts_type *transform, fftw_workspace *fftw) {
    int d = index->settings->n_segments;
    if (index->settings->function_type == 4) {
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
    FILE *file = fopen(path, "rb");
    if (file == NULL) return FAILURE;
    rawfile = malloc((size_t) ts_num * index->settings->timeseries_size * sizeof(*rawfile));
    if (rawfile == NULL || fread(rawfile, sizeof(*rawfile),
        (size_t) ts_num * index->settings->timeseries_size, file) !=
        (size_t) ts_num * index->settings->timeseries_size) { fclose(file); free(rawfile); rawfile = NULL; return FAILURE; }
    fclose(file);
    struct symbolic_trie_index *trie = calloc(1, sizeof(*trie));
    if (trie == NULL || (trie->root = trie_node_create(index->settings->n_segments, NULL, 0)) == NULL) {
        free(trie); free(rawfile); rawfile = NULL; return FAILURE;
    }
    trie->dimensions = index->settings->n_segments;
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
    for (long i = 0; i < ts_num; ++i) trie_node_update_mbb(trie->root, trie->root->words[i], trie->dimensions);
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#pragma omp single
#endif
    trie_build_subtree(trie, index, trie->root, 0);
    trie->node_count = trie_node_count(trie->root);
    index->trie = trie;
    index->total_records = ts_num;
    return SUCCESS;
}

static float trie_search_node(isax_index *index, const symbolic_trie_node *node, const ts_type *query,
                              const ts_type *transform, float bsf, trie_query_stats *stats,
                              const symbolic_trie_node *skip_leaf) {
    if (node == skip_leaf) return bsf;
    if (stats != NULL) ++stats->checked_nodes;
    float lower = messi_minidist_range_raw(index, (float *) transform, node->min_word, node->max_word,
                                           index->settings->max_sax_cardinalities, bsf);
    if (lower >= bsf) return bsf;
    if (!node->leaf) {
        typedef struct { const symbolic_trie_node *node; float bound; } trie_child_bound;
        trie_child_bound ordered[8];
        int n = 0;
        for (int i = 0; i < 8; ++i) if (node->children[i] != NULL) {
            ordered[n].node = node->children[i];
            ordered[n++].bound = messi_minidist_range_raw(index, (float *) transform,
                node->children[i]->min_word, node->children[i]->max_word,
                index->settings->max_sax_cardinalities, bsf);
        }
        for (int i = 1; i < n; ++i) {
            trie_child_bound current = ordered[i]; int j = i - 1;
            while (j >= 0 && ordered[j].bound > current.bound) { ordered[j + 1] = ordered[j]; --j; }
            ordered[j + 1] = current;
        }
        for (int i = 0; i < n; ++i) if (ordered[i].bound < bsf)
            bsf = trie_search_node(index, ordered[i].node, query, transform, bsf, stats, skip_leaf);
        return bsf;
    }
    for (int i = 0; i < node->size; ++i) {
        if (stats != NULL) ++stats->lower_bounds;
        else __sync_fetch_and_add(&LBDcalculationnumber, 1);
        if (messi_minidist_raw(index, (float *) transform, node->words[i], index->settings->max_sax_cardinalities, bsf) < bsf) {
            if (stats != NULL) ++stats->exact_distances;
            else __sync_fetch_and_add(&RDcalculationnumber, 1);
            float d = ts_ed((ts_type *) query, rawfile + node->positions[i], index->settings->timeseries_size, bsf);
            if (d < bsf) bsf = d;
        }
    }
    return bsf;
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
                float bound = messi_minidist_range_raw(index, (float *) transform,
                                                        node->children[i]->min_word,
                                                        node->children[i]->max_word,
                                                        index->settings->max_sax_cardinalities,
                                                        best_bound);
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
    if (messi_minidist_range_raw(index, (float *) transform, node->min_word, node->max_word,
                                 index->settings->max_sax_cardinalities, bsf) >= bsf) return;
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
        if (messi_minidist_raw(index, (float *) transform, node->words[i], index->settings->max_sax_cardinalities, bsf) < bsf) {
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
                trie_search_node(index, seed_leaf, query, transform, FLT_MAX, &stats, NULL);
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
                                        minimum_distance < bsf ? minimum_distance : bsf, &stats, seed_leaf);
        }
        stats.total_microseconds = trie_monotonic_microseconds() - query_start;
        trie_save_query_stats(index->trie, &stats, distance);
        if (i == 0) PRINT_STATS_HEADER();
        trie_print_query_stats(i, index->trie, &stats, distance);
    }
    fclose(file); free(query); free(transform); free(word); fftw_workspace_destroy(&fftw); return SUCCESS;
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
                        trie_search_node(index, seed_leaf, query, transform, FLT_MAX, &stats[i], NULL);
                stats[i].approximate_distance = bsf;
                if (minimum_distance < bsf) bsf = minimum_distance;
                distances[i] = trie_search_node(index, index->trie->root, query, transform, bsf,
                                                &stats[i], seed_leaf);
                stats[i].total_microseconds = trie_monotonic_microseconds() - query_start;
            }
        }
        free(transform); free(word);
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
    if (maxquerythread <= 1) result.distance = trie_search_node(index, index->trie->root, query,
                                                                  transform, bsf, NULL, NULL);
    else result.distance = trie_parallel_exact_search(index, query, transform, bsf, NULL);
    return result;
}

void symbolic_trie_destroy(isax_index *index) {
    if (index != NULL && index->trie != NULL) {
        trie_node_destroy(index->trie->root);
        free(index->trie); index->trie = NULL;
    }
}
