#include "config.h"
#include "globals.h"
#include "ads/trie/trie.h"
#include "ads/trie/radial_bound.h"

#define TRIE_MAX_FANOUT 256
#include "ads/calc_utils.h"
#include "ads/lower_bound_simd.h"
#include "ads/sax/sax.h"
#include "ads/sfa/dft.h"
#include "ads/spartan/pca.h"
#include "ads/spartan/spartan.h"
#include "ads/pisa/pisa.h"
#include "ads/inmemory_index_engine.h"
#include "ads/build_progress.h"

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

/* Keep large diagnostic counts compact and easy to scan in benchmark logs. */
static const char *trie_format_compact_count(unsigned long long value,
                                              char buffer[32]) {
    if (value >= 1000000000ULL)
        snprintf(buffer, 32, "%.3g G", (double) value / 1e9);
    else if (value >= 1000000ULL)
        snprintf(buffer, 32, "%.3g M", (double) value / 1e6);
    else if (value >= 1000ULL)
        snprintf(buffer, 32, "%.3g K", (double) value / 1e3);
    else
        snprintf(buffer, 32, "%llu", value);
    return buffer;
}

typedef struct {
    uint64_t low;
    uint64_t high;
} trie_dimension_mask;

typedef struct {
    int offset;
    int size;
} trie_leaf_cluster;

typedef struct symbolic_trie_node {
    struct symbolic_trie_node *children[TRIE_MAX_FANOUT];
    struct symbolic_trie_node *parent;
    sax_type *min_word;
    sax_type *max_word;
    sax_type **words;
    file_position_type *positions;
    int size;
    int capacity;
    int split_dimension;
    uint16_t split_fanout;
    /* The trie may split on up to 128 transform dimensions. */
    trie_dimension_mask used_dimensions;
    unsigned char leaf;
    unsigned char mbb_valid;
    unsigned char split_exhausted;
    trie_leaf_cluster *clusters;
    sax_type *cluster_min_words;
    sax_type *cluster_max_words;
    /* A ball in the exact (already normalized, when requested) raw-series
     * space.  Unlike the symbolic k-means centres used to form clusters,
     * these are certified for pruning by the triangle inequality. */
    ts_type *cluster_raw_centroids;
    float *cluster_raw_radii;
    /* Per-record centroid radii in the existing IVF-cluster record order. */
    float *record_raw_radii;
    int cluster_count;
} symbolic_trie_node;

struct symbolic_trie_index {
    symbolic_trie_node *root;
    /* All record words are slices of this one allocation.  Leaf splits move
     * pointers into it, so teardown has one free rather than one per series. */
    sax_type *word_arena;
    int dimensions;
    int bound_dimensions;
    int split_dimensions;
    int fanout;
    int dynamic_alphabet;
    unsigned long split_exhausted_leaves;
    unsigned long node_count;
    unsigned long clustered_leaves;
    unsigned long long cluster_count;
    unsigned long radial_leaves;
    unsigned long long radial_records;
};

static int trie_dimension_is_used(const symbolic_trie_node *node, int dimension) {
    return dimension < 64
               ? (node->used_dimensions.low & (UINT64_C(1) << dimension)) != 0
               : (node->used_dimensions.high & (UINT64_C(1) << (dimension - 64))) != 0;
}

static trie_dimension_mask trie_dimension_mark_used(trie_dimension_mask used,
                                                    int dimension) {
    if (dimension < 64) used.low |= UINT64_C(1) << dimension;
    else used.high |= UINT64_C(1) << (dimension - 64);
    return used;
}

static int trie_bucket(const sax_type *word, int dimension, int fanout) {
    int bits = 0;
    for (int value = fanout; value > 1; value >>= 1) ++bits;
    return word[dimension] >> (8 - bits);
}

static int trie_dimension_fanout(const struct symbolic_trie_index *trie,
                                 const isax_index *index, int dimension) {
    if (!trie->dynamic_alphabet) return trie->fanout;
    const int bits = index->settings->root_bit_cardinalities[dimension];
    return bits > 0 ? 1 << bits : 1;
}

/* A full entropy pass over every dimension of a very large root is
 * prohibitively expensive (64 complete scans at the current maximum word
 * length).  A deterministic, evenly-spaced sample retains the distribution
 * information needed for split selection while keeping the expensive full
 * pass to the one partition pass that cannot be avoided. */
#define TRIE_ROOT_SPLIT_SAMPLE_SIZE 1000000L
#define TRIE_PCA_PROJECTION_BLOCK_RECORDS 4096U

typedef struct trie_query_stats {
    unsigned long checked_nodes;
    unsigned long lower_bounds;
    unsigned long exact_distances;
    unsigned long long total_microseconds;
    unsigned long long mbr_bound_microseconds;
    unsigned long long record_bound_microseconds;
    unsigned long long exact_distance_microseconds;
    unsigned long long traversal_microseconds;
    unsigned long long frontier_microseconds;
    unsigned long long queue_microseconds;
    unsigned long long candidate_heap_microseconds;
    unsigned long long synchronization_microseconds;
    unsigned long cluster_bounds;
    unsigned long cluster_pruned;
    unsigned long cluster_records_pruned;
    unsigned long cluster_symbolic_pruned;
    unsigned long cluster_symbolic_records_pruned;
    unsigned long cluster_raw_ball_bounds;
    unsigned long cluster_raw_ball_pruned;
    unsigned long cluster_raw_ball_records_pruned;
    unsigned long radial_candidates;
    unsigned long radial_pruned;
    float approximate_distance;
} trie_query_stats;
static unsigned long long trie_monotonic_microseconds(void);

/* The index settings describe the full split word (normally 64 dimensions).
 * Lower bounds intentionally see only the configured bound prefix.  A local
 * shallow copy is safe for concurrent queries and reuses the representation-
 * specific SFA/SPARTAN/PISA/SAX bound implementations unchanged. */
static float trie_lower_bound(const struct symbolic_trie_index *trie, isax_index *index,
                              const ts_type *transform, sax_type *sax_min,
                              sax_type *sax_max, int dimensions, float bsf,
                              trie_query_stats *stats) {
    unsigned long long start = 0;
    if (profile_query_phases && stats != NULL) start = trie_monotonic_microseconds();
    isax_index shadow_index = *index;
    isax_index_settings shadow_settings = *index->settings;
    shadow_settings.n_segments = dimensions;
    shadow_index.settings = &shadow_settings;
    float result = messi_minidist_range_raw(&shadow_index, (float *) transform, sax_min, sax_max,
                                            shadow_settings.max_sax_cardinalities, bsf);
    if (start != 0) stats->mbr_bound_microseconds += trie_monotonic_microseconds() - start;
    return result;
}

typedef struct {
    float lower_bound;
    int record_index;
} trie_leaf_candidate;

#define TRIE_RADIAL_AUTO_SAMPLE_CANDIDATES 65536UL
#define TRIE_RADIAL_AUTO_REJECTION_DENOMINATOR 4UL

typedef struct {
    pthread_mutex_t lock;
    unsigned long candidates;
    unsigned long pruned;
    /* 0 samples, 1 keeps the radial bound, -1 bypasses it. */
    int decision;
} trie_radial_auto_gate;

/* One scratch heap per batch-query worker.  It is retained across leaves and
 * avoids per-record allocation in best-first leaf refinement. */
typedef struct {
    trie_leaf_candidate *candidates;
    int capacity;
    float record_lb_table[MESSI_RECORD_LB_MAX_DIMENSIONS][256];
    int record_lb_table_ready;
    trie_radial_auto_gate *radial_auto_gate;
} trie_query_scratch;

static int trie_radial_auto_gate_init(const isax_index *index,
                                      trie_radial_auto_gate *gate) {
    if (index == NULL || gate == NULL ||
        !index->settings->trie_leaf_ivf_radial_bound_auto) return 0;
    gate->candidates = 0;
    gate->pruned = 0;
    gate->decision = 0;
    return pthread_mutex_init(&gate->lock, NULL) == 0;
}

static void trie_radial_auto_gate_destroy(trie_radial_auto_gate *gate,
                                          int initialized) {
    if (initialized) pthread_mutex_destroy(&gate->lock);
}

static inline int trie_radial_auto_gate_allows(const trie_query_scratch *scratch) {
    return scratch != NULL && scratch->radial_auto_gate != NULL &&
           __atomic_load_n(&scratch->radial_auto_gate->decision,
                           __ATOMIC_ACQUIRE) >= 0;
}

/* This runs once per scanned IVF range only during the bounded calibration
 * prefix.  No adaptive counters or synchronization remain after a decision. */
static void trie_radial_auto_gate_observe(trie_radial_auto_gate *gate,
                                          unsigned long candidates,
                                          unsigned long pruned) {
    if (gate == NULL || candidates == 0 ||
        __atomic_load_n(&gate->decision, __ATOMIC_ACQUIRE) != 0) return;
    pthread_mutex_lock(&gate->lock);
    if (__atomic_load_n(&gate->decision, __ATOMIC_RELAXED) == 0) {
        gate->candidates += candidates;
        gate->pruned += pruned;
        if (gate->candidates >= TRIE_RADIAL_AUTO_SAMPLE_CANDIDATES) {
            const unsigned long required =
                (gate->candidates + TRIE_RADIAL_AUTO_REJECTION_DENOMINATOR - 1) /
                TRIE_RADIAL_AUTO_REJECTION_DENOMINATOR;
            __atomic_store_n(&gate->decision,
                             gate->pruned >= required ? 1 : -1,
                             __ATOMIC_RELEASE);
        }
    }
    pthread_mutex_unlock(&gate->lock);
}

/* Leaf records have one concrete word, unlike an MBR range.  The query-local
 * symbol table covers every supported record-bound prefix (16--64 dims). */
static float trie_record_lower_bound(const struct symbolic_trie_index *trie,
                                     isax_index *index, const ts_type *transform,
                                     sax_type *word, float bsf, float mbr_suffix,
                                     const trie_query_scratch *scratch) {
    if (scratch != NULL && scratch->record_lb_table_ready) {
        const float record_limit = bsf - mbr_suffix;
        if (record_limit <= 0.0f) return mbr_suffix;
        const float distance = messi_record_lb_table_sum(scratch->record_lb_table, word,
            trie->bound_dimensions, record_limit, index->settings->SIMD_flag != 0);
        return mbr_suffix + distance;
    }
    isax_index shadow_index = *index;
    isax_index_settings shadow_settings = *index->settings;
    shadow_settings.n_segments = trie->bound_dimensions;
    shadow_index.settings = &shadow_settings;
    return messi_minidist_raw(&shadow_index, (float *) transform, word,
                               shadow_settings.max_sax_cardinalities, bsf);
}

static float trie_record_mbr_suffix(const struct symbolic_trie_index *trie,
                                    isax_index *index, const ts_type *transform,
                                    const sax_type *min_word, const sax_type *max_word,
                                    const trie_query_scratch *scratch) {
    if (!index->settings->trie_record_mbr_suffix_bound || scratch == NULL ||
        !scratch->record_lb_table_ready ||
        trie->dimensions <= trie->bound_dimensions) return 0.0f;
    isax_index shadow_index = *index;
    isax_index_settings shadow_settings = *index->settings;
    shadow_settings.n_segments = trie->dimensions;
    shadow_index.settings = &shadow_settings;
    float suffix = 0.0f;
    (void) messi_minidist_range_raw_partitioned(&shadow_index, (float *) transform,
                                                (sax_type *) min_word, (sax_type *) max_word,
                                                shadow_settings.max_sax_cardinalities,
                                                FLT_MAX, trie->bound_dimensions, &suffix);
    return suffix;
}

static void trie_prepare_record_lb_table(const struct symbolic_trie_index *trie,
                                         const isax_index *index,
                                         const ts_type *transform,
                                         trie_query_scratch *scratch) {
    if (scratch == NULL) return;
    scratch->record_lb_table_ready =
        messi_build_record_lb_table(index, transform, trie->bound_dimensions,
                                    scratch->record_lb_table);
}

/* FFTW plan creation is not thread-safe on all supported builds. */
static pthread_mutex_t trie_fftw_plan_lock = PTHREAD_MUTEX_INITIALIZER;

static symbolic_trie_node *trie_node_create(int dimensions, symbolic_trie_node *parent,
                                            trie_dimension_mask used_dimensions) {
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
    node->split_fanout = 0;
    node->leaf = 1;
    return node;
}

static void trie_node_destroy(symbolic_trie_node *node) {
    if (node == NULL) return;
    for (int i = 0; i < TRIE_MAX_FANOUT; ++i) trie_node_destroy(node->children[i]);
    free(node->words); free(node->positions); free(node->min_word); free(node->max_word);
    free(node->clusters); free(node->cluster_min_words); free(node->cluster_max_words);
    free(node->cluster_raw_centroids); free(node->cluster_raw_radii);
    free(node->record_raw_radii); free(node);
}

static unsigned long trie_node_count(const symbolic_trie_node *node) {
    if (node == NULL) return 0;
    unsigned long count = 1;
    for (int i = 0; i < TRIE_MAX_FANOUT; ++i) count += trie_node_count(node->children[i]);
    return count;
}

/* Compute split quality from the completed tree so parallel construction does
 * not need shared diagnostic counters.  Child record counts are recovered
 * from their subtrees because internal nodes release their leaf arrays. */
static unsigned long trie_collect_split_diagnostics(const symbolic_trie_node *node,
                                                    unsigned long *internal_nodes,
                                                    unsigned long long *nonempty_children,
                                                    unsigned long long *child_slots,
                                                    double *largest_child_share) {
    if (node == NULL) return 0;
    if (node->leaf) return (unsigned long) node->size;

    unsigned long total = 0, largest = 0;
    int occupied = 0;
    for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL) {
        unsigned long child_records = trie_collect_split_diagnostics(node->children[i], internal_nodes,
                                                                       nonempty_children, child_slots,
                                                                       largest_child_share);
        total += child_records;
        if (child_records != 0) {
            ++occupied;
            if (child_records > largest) largest = child_records;
        }
    }
    ++*internal_nodes;
    *nonempty_children += (unsigned long long) occupied;
    *child_slots += node->split_fanout;
    if (total != 0) *largest_child_share += (double) largest / (double) total;
    return total;
}

static unsigned long long trie_monotonic_microseconds(void) {
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    return (unsigned long long) time.tv_sec * 1000000ULL + (unsigned long long) time.tv_nsec / 1000ULL;
}

static void trie_print_query_stats(int query_index, const struct symbolic_trie_index *trie,
                                   const trie_query_stats *stats, float distance,
                                   unsigned long long cumulative_microseconds) {
    fprintf(stderr, "%3d: %8lu %13lu %20.3f %12.3f %12.3f\n",
           query_index, trie->node_count, stats->checked_nodes,
           stats->approximate_distance, distance,
           (double) cumulative_microseconds / 1000.0);
}

static void trie_save_query_stats(const struct symbolic_trie_index *trie,
                                  const trie_query_stats *stats, float distance,
                                  unsigned long long cumulative_microseconds) {
    total_tree_nodes = (int) trie->node_count;
    checked_nodes = (int) stats->checked_nodes;
    BYTES_ACCESSED = 0;
    loaded_nodes = 0;
    loaded_records = 0;
    APPROXIMATE = stats->approximate_distance;
    total_querying_time = (double) stats->total_microseconds;
    total_init_time = 0.0;
    total_tree_pass_time = (double) stats->traversal_microseconds;
    TOTAL_PQ_INSERT_TIME = 0;
    TOTAL_PQ_REMOVE_TIME = 0;
    TOTAL_MBR_DIST_CALC_TIME = (unsigned long) stats->mbr_bound_microseconds;
    TOTAL_RECORD_LB_DIST_CALC_TIME = (unsigned long) stats->record_bound_microseconds;
    TOTAL_LB_DIST_CALC_TIME = TOTAL_MBR_DIST_CALC_TIME + TOTAL_RECORD_LB_DIST_CALC_TIME;
    TOTAL_REAL_DIST_CALC_TIME = (unsigned long) stats->exact_distance_microseconds;
    TOTAL_TRIE_FRONTIER_TIME = (unsigned long) stats->frontier_microseconds;
    TOTAL_TRIE_QUEUE_TIME = (unsigned long) stats->queue_microseconds;
    TOTAL_TRIE_HEAP_TIME = (unsigned long) stats->candidate_heap_microseconds;
    TOTAL_TRIE_SYNC_TIME = (unsigned long) stats->synchronization_microseconds;
    total_time = (double) cumulative_microseconds;
    LBDcalculationnumber = stats->lower_bounds;
    RDcalculationnumber = stats->exact_distances;
    trie_cluster_bounds = stats->cluster_bounds;
    trie_cluster_pruned = stats->cluster_pruned;
    trie_cluster_records_pruned = stats->cluster_records_pruned;
    trie_cluster_symbolic_pruned = stats->cluster_symbolic_pruned;
    trie_cluster_symbolic_records_pruned = stats->cluster_symbolic_records_pruned;
    trie_cluster_raw_ball_bounds = stats->cluster_raw_ball_bounds;
    trie_cluster_raw_ball_pruned = stats->cluster_raw_ball_pruned;
    trie_cluster_raw_ball_records_pruned = stats->cluster_raw_ball_records_pruned;
    trie_radial_candidates = stats->radial_candidates;
    trie_radial_pruned = stats->radial_pruned;
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

/* Root partitioning already computes exact MBRs for all populated children.
 * Merge them instead of making a second full scan over the root's records. */
static void trie_node_merge_child_mbbs(symbolic_trie_node *node,
                                       symbolic_trie_node *const *children,
                                       int fanout, int dimensions) {
    node->mbb_valid = 0;
    for (int bucket = 0; bucket < fanout; ++bucket) {
        const symbolic_trie_node *child = children[bucket];
        if (child == NULL || !child->mbb_valid) continue;
        if (!node->mbb_valid) {
            memcpy(node->min_word, child->min_word, (size_t) dimensions);
            memcpy(node->max_word, child->max_word, (size_t) dimensions);
            node->mbb_valid = 1;
            continue;
        }
        for (int d = 0; d < dimensions; ++d) {
            if (child->min_word[d] < node->min_word[d]) node->min_word[d] = child->min_word[d];
            if (child->max_word[d] > node->max_word[d]) node->max_word[d] = child->max_word[d];
        }
    }
}

#define TRIE_LEAF_IVF_MIN_SIZE 4096
#define TRIE_LEAF_IVF_TRAIN_SIZE 2048
#define TRIE_LEAF_IVF_MAX_ITERATIONS 10

static double trie_cluster_distance(const sax_type *word, const float *centroid,
                                    const float *centers, int dimensions, int alphabet) {
    double distance = 0.0;
    for (int d = 0; d < dimensions; ++d) {
        const double delta = (double) centers[(size_t) d * alphabet + word[d]] - centroid[d];
        distance += delta * delta;
    }
    return distance;
}

static int trie_cluster_leaf(symbolic_trie_node *node, int dimensions, int alphabet,
                             int cluster_count, const float *centers, int ts_length,
                             int radial_enabled) {
    const int size = node->size;
    if (size < TRIE_LEAF_IVF_MIN_SIZE || cluster_count >= size || centers == NULL) return 1;
    const int train_size = size < TRIE_LEAF_IVF_TRAIN_SIZE ? size : TRIE_LEAF_IVF_TRAIN_SIZE;
    float *centroids = calloc((size_t) cluster_count * dimensions, sizeof(*centroids));
    float *sums = calloc((size_t) cluster_count * dimensions, sizeof(*sums));
    int *assignments = malloc((size_t) size * sizeof(*assignments));
    int *counts = calloc((size_t) cluster_count, sizeof(*counts));
    if (centroids == NULL || sums == NULL || assignments == NULL || counts == NULL) {
        free(centroids); free(sums); free(assignments); free(counts); return 0;
    }
    for (int cluster = 0; cluster < cluster_count; ++cluster) {
        const int record = (int) (((long long) cluster * size) / cluster_count);
        for (int d = 0; d < dimensions; ++d)
            centroids[(size_t) cluster * dimensions + d] =
                centers[(size_t) d * alphabet + node->words[record][d]];
    }
    for (int i = 0; i < size; ++i) assignments[i] = -1;

    for (int iteration = 0; iteration < TRIE_LEAF_IVF_MAX_ITERATIONS; ++iteration) {
        memset(sums, 0, (size_t) cluster_count * dimensions * sizeof(*sums));
        memset(counts, 0, (size_t) cluster_count * sizeof(*counts));
        int changed = 0;
        for (int sample = 0; sample < train_size; ++sample) {
            const int record = (int) (((long long) sample * size) / train_size);
            const sax_type *word = node->words[record];
            int best = 0;
            double best_distance = trie_cluster_distance(word, centroids, centers, dimensions, alphabet);
            for (int cluster = 1; cluster < cluster_count; ++cluster) {
                const double distance = trie_cluster_distance(word,
                    centroids + (size_t) cluster * dimensions, centers, dimensions, alphabet);
                if (distance < best_distance) { best_distance = distance; best = cluster; }
            }
            if (assignments[record] != best) { assignments[record] = best; changed = 1; }
            ++counts[best];
            for (int d = 0; d < dimensions; ++d)
                sums[(size_t) best * dimensions + d] += centers[(size_t) d * alphabet + word[d]];
        }
        for (int cluster = 0; cluster < cluster_count; ++cluster) {
            if (counts[cluster] == 0) {
                const int record = (int) (((long long) cluster * size) / cluster_count);
                for (int d = 0; d < dimensions; ++d)
                    centroids[(size_t) cluster * dimensions + d] =
                        centers[(size_t) d * alphabet + node->words[record][d]];
            } else {
                for (int d = 0; d < dimensions; ++d)
                    centroids[(size_t) cluster * dimensions + d] =
                        sums[(size_t) cluster * dimensions + d] / counts[cluster];
            }
        }
        if (!changed) break;
    }

    memset(counts, 0, (size_t) cluster_count * sizeof(*counts));
    for (int record = 0; record < size; ++record) {
        const sax_type *word = node->words[record];
        int best = 0;
        double best_distance = trie_cluster_distance(word, centroids, centers, dimensions, alphabet);
        for (int cluster = 1; cluster < cluster_count; ++cluster) {
            const double distance = trie_cluster_distance(word,
                centroids + (size_t) cluster * dimensions, centers, dimensions, alphabet);
            if (distance < best_distance) { best_distance = distance; best = cluster; }
        }
        assignments[record] = best;
        ++counts[best];
    }
    int active = 0;
    for (int cluster = 0; cluster < cluster_count; ++cluster) if (counts[cluster] != 0) ++active;
    if (active < 2) { free(centroids); free(sums); free(assignments); free(counts); return 1; }

    trie_leaf_cluster *clusters = calloc((size_t) active, sizeof(*clusters));
    sax_type *min_words = malloc((size_t) active * dimensions);
    sax_type *max_words = calloc((size_t) active * dimensions, sizeof(*max_words));
    sax_type **ordered_words = malloc((size_t) size * sizeof(*ordered_words));
    file_position_type *ordered_positions = malloc((size_t) size * sizeof(*ordered_positions));
    ts_type *raw_centroids = calloc((size_t) active * ts_length, sizeof(*raw_centroids));
    double *raw_sums = calloc((size_t) active * ts_length, sizeof(*raw_sums));
    float *raw_radii = calloc((size_t) active, sizeof(*raw_radii));
    float *record_radii = radial_enabled ? malloc((size_t) size * sizeof(*record_radii)) : NULL;
    int *offsets = calloc((size_t) cluster_count, sizeof(*offsets));
    int *cursors = calloc((size_t) cluster_count, sizeof(*cursors));
    int *remap = malloc((size_t) cluster_count * sizeof(*remap));
    if (clusters == NULL || min_words == NULL || max_words == NULL || ordered_words == NULL ||
        ordered_positions == NULL || raw_centroids == NULL || raw_sums == NULL || raw_radii == NULL ||
        offsets == NULL || cursors == NULL || remap == NULL) {
        free(clusters); free(min_words); free(max_words); free(ordered_words); free(ordered_positions);
        free(raw_centroids); free(raw_sums); free(raw_radii); free(record_radii);
        free(offsets); free(cursors); free(remap); free(centroids); free(sums); free(assignments); free(counts);
        return 0;
    }
    int offset = 0, compact = 0;
    for (int cluster = 0; cluster < cluster_count; ++cluster) {
        offsets[cluster] = offset;
        cursors[cluster] = offset;
        if (counts[cluster] == 0) { remap[cluster] = -1; continue; }
        remap[cluster] = compact;
        clusters[compact].offset = offset;
        clusters[compact].size = counts[cluster];
        memset(min_words + (size_t) compact * dimensions, UCHAR_MAX, (size_t) dimensions);
        offset += counts[cluster]; ++compact;
    }
    for (int record = 0; record < size; ++record) {
        const int cluster = assignments[record];
        const int position = cursors[cluster]++;
        ordered_words[position] = node->words[record];
        ordered_positions[position] = node->positions[record];
        const int group = remap[cluster];
        sax_type *minimum = min_words + (size_t) group * dimensions;
        sax_type *maximum = max_words + (size_t) group * dimensions;
        for (int d = 0; d < dimensions; ++d) {
            if (node->words[record][d] < minimum[d]) minimum[d] = node->words[record][d];
            if (node->words[record][d] > maximum[d]) maximum[d] = node->words[record][d];
        }
        const ts_type *series = rawfile + node->positions[record];
        double *sum = raw_sums + (size_t) group * ts_length;
        for (int t = 0; t < ts_length; ++t) sum[t] += series[t];
    }
    for (int cluster = 0; cluster < active; ++cluster) {
        ts_type *centroid = raw_centroids + (size_t) cluster * ts_length;
        const double *sum = raw_sums + (size_t) cluster * ts_length;
        const int count = clusters[cluster].size;
        for (int t = 0; t < ts_length; ++t) centroid[t] = (ts_type) (sum[t] / count);
    }
    for (int cluster = 0; cluster < active; ++cluster) {
        const ts_type *centroid = raw_centroids + (size_t) cluster * ts_length;
        double maximum_distance_squared = 0.0;
        const int first = clusters[cluster].offset;
        const int last = first + clusters[cluster].size;
        for (int record = first; record < last; ++record) {
            const ts_type *series = rawfile + ordered_positions[record];
            double distance_squared = 0.0;
            for (int t = 0; t < ts_length; ++t) {
                const double delta = (double) series[t] - centroid[t];
                distance_squared += delta * delta;
            }
            if (record_radii != NULL)
                record_radii[record] = (float) sqrt(distance_squared);
            if (distance_squared > maximum_distance_squared)
                maximum_distance_squared = distance_squared;
        }
        /* Round upward: an undersized stored radius could invalidate the
         * lower bound after the double-to-float conversion. */
        raw_radii[cluster] = nextafterf((float) sqrt(maximum_distance_squared), INFINITY);
    }
    memcpy(node->words, ordered_words, (size_t) size * sizeof(*node->words));
    memcpy(node->positions, ordered_positions, (size_t) size * sizeof(*node->positions));
    node->clusters = clusters;
    node->cluster_min_words = min_words;
    node->cluster_max_words = max_words;
    node->cluster_raw_centroids = raw_centroids;
    node->cluster_raw_radii = raw_radii;
    node->record_raw_radii = record_radii;
    node->cluster_count = active;
    free(ordered_words); free(ordered_positions); free(offsets); free(cursors); free(remap);
    free(raw_sums);
    free(centroids); free(sums); free(assignments); free(counts);
    return 1;
}

typedef struct {
    symbolic_trie_node **items;
    size_t size;
    size_t capacity;
} trie_cluster_leaf_list;

/* The tree is already immutable at this point. First collect eligible leaves
 * on one thread, then let workers cluster independent leaves without sharing
 * per-leaf scratch state or modifying common tree counters. */
static int trie_collect_cluster_leaves(symbolic_trie_node *node,
                                       trie_cluster_leaf_list *leaves) {
    if (node == NULL) return 1;
    if (!node->leaf) {
        for (int i = 0; i < node->split_fanout; ++i)
            if (!trie_collect_cluster_leaves(node->children[i], leaves)) return 0;
        return 1;
    }
    if (node->size < TRIE_LEAF_IVF_MIN_SIZE) return 1;
    if (leaves->size == leaves->capacity) {
        size_t capacity = leaves->capacity == 0 ? 128 : leaves->capacity * 2;
        if (capacity > SIZE_MAX / sizeof(*leaves->items)) return 0;
        symbolic_trie_node **items = realloc(leaves->items, capacity * sizeof(*items));
        if (items == NULL) return 0;
        leaves->items = items;
        leaves->capacity = capacity;
    }
    leaves->items[leaves->size++] = node;
    return 1;
}

static int trie_cluster_leaves_parallel(struct symbolic_trie_index *trie,
                                        const isax_index *index, const float *centers,
                                        int *active_workers, size_t *eligible_leaves) {
    trie_cluster_leaf_list leaves = {0};
    if (!trie_collect_cluster_leaves(trie->root, &leaves)) {
        free(leaves.items);
        return 0;
    }
    if (eligible_leaves != NULL) *eligible_leaves = leaves.size;

    const int requested_workers = maxquerythread > 0 ? maxquerythread : 1;
    int workers = 1;
    unsigned long clustered_leaves = 0;
    unsigned long long cluster_count = 0;
    unsigned long radial_leaves = 0;
    unsigned long long radial_records = 0;
    int failed = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(requested_workers) reduction(+:clustered_leaves, cluster_count, radial_leaves, radial_records) reduction(|:failed)
    {
#pragma omp single
        workers = omp_get_num_threads();
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < leaves.size; ++i) {
            symbolic_trie_node *node = leaves.items[i];
            if (!trie_cluster_leaf(node, trie->dimensions, index->settings->sax_alphabet_cardinality,
                                   index->settings->trie_leaf_ivf, centers,
                                   index->settings->timeseries_size,
                                   index->settings->trie_leaf_ivf_radial_bound)) {
                failed = 1;
                continue;
            }
            if (node->cluster_count != 0) {
                ++clustered_leaves;
                cluster_count += (unsigned long long) node->cluster_count;
                if (node->record_raw_radii != NULL) {
                    ++radial_leaves;
                    radial_records += (unsigned long long) node->size;
                }
            }
        }
    }
#else
    for (size_t i = 0; i < leaves.size; ++i) {
        symbolic_trie_node *node = leaves.items[i];
        if (!trie_cluster_leaf(node, trie->dimensions, index->settings->sax_alphabet_cardinality,
                               index->settings->trie_leaf_ivf, centers,
                               index->settings->timeseries_size,
                               index->settings->trie_leaf_ivf_radial_bound)) {
            failed = 1;
            continue;
        }
        if (node->cluster_count != 0) {
            ++clustered_leaves;
            cluster_count += (unsigned long long) node->cluster_count;
            if (node->record_raw_radii != NULL) {
                ++radial_leaves;
                radial_records += (unsigned long long) node->size;
            }
        }
    }
#endif
    free(leaves.items);
    trie->clustered_leaves = clustered_leaves;
    trie->cluster_count = cluster_count;
    trie->radial_leaves = radial_leaves;
    trie->radial_records = radial_records;
    if (active_workers != NULL) *active_workers = workers;
    return !failed;
}

static float *trie_build_bin_centers(const isax_index *index, int dimensions) {
    const int alphabet = index->settings->sax_alphabet_cardinality;
    if (index->bins == NULL || alphabet < 2) return NULL;
    float *centers = malloc((size_t) dimensions * alphabet * sizeof(*centers));
    if (centers == NULL) return NULL;
    for (int d = 0; d < dimensions; ++d) {
        const float *bins = index->bins[d];
        if (bins == NULL) { free(centers); return NULL; }
        for (int symbol = 0; symbol < alphabet; ++symbol) {
            float lower, upper;
            if (symbol == 0) {
                const float width = alphabet > 2 && bins[1] > bins[0] ? bins[1] - bins[0] : 0.0f;
                lower = bins[0] - 0.5f * width; upper = bins[0];
            } else if (symbol == alphabet - 1) {
                const float width = bins[alphabet - 2] > bins[alphabet - 3]
                                        ? bins[alphabet - 2] - bins[alphabet - 3] : 0.0f;
                lower = bins[alphabet - 2]; upper = bins[alphabet - 2] + 0.5f * width;
            } else { lower = bins[symbol - 1]; upper = bins[symbol]; }
            centers[(size_t) d * alphabet + symbol] = isfinite(lower) && isfinite(upper)
                                                           ? 0.5f * (lower + upper) : (float) symbol;
        }
    }
    return centers;
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

static float trie_exact_distance(const ts_type *query, const ts_type *record,
                                 int length, float bsf, trie_query_stats *stats) {
    unsigned long long start = 0;
    if (profile_query_phases && stats != NULL) start = trie_monotonic_microseconds();
    float result = ts_ed((ts_type *) query, (ts_type *) record, length, bsf);
    if (start != 0) stats->exact_distance_microseconds += trie_monotonic_microseconds() - start;
    return result;
}

/* For every member x of a cluster, ||q - x|| is at least
 * max(0, ||q - c|| - r), where r encloses every member around c.  The
 * cluster is still formed in symbolic space; this separate raw-space ball is
 * what makes the additional bound exact. */
static messi_distance_interval trie_cluster_query_radius(const symbolic_trie_node *node,
                                                         int cluster, const ts_type *query,
                                                         int ts_length) {
    const ts_type *centroid = node->cluster_raw_centroids + (size_t) cluster * ts_length;
    double squared = 0.0;
    for (int t = 0; t < ts_length; ++t) {
        const double delta = (double) query[t] - centroid[t];
        squared += delta * delta;
    }
    return messi_distance_interval_from_squared(squared, ts_length);
}

static float trie_cluster_raw_ball_lower_bound(const symbolic_trie_node *node,
                                               int cluster,
                                               messi_distance_interval query_radius) {
    if (node->cluster_raw_centroids == NULL || node->cluster_raw_radii == NULL)
        return 0.0f;
    /* Downward rounding keeps finite-precision arithmetic from turning this
     * mathematical lower bound into an overestimate.  Radii were rounded up
     * when built. */
    const double lower = query_radius.lower - (double) node->cluster_raw_radii[cluster];
    if (lower <= 0.0) return 0.0f;
    return nextafterf((float) (lower * lower), 0.0f);
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
                             const struct symbolic_trie_index *trie,
                             const isax_index *index, const double *variances) {
    int best = -1;
    double best_score = -1.0;
    for (int d = 0; d < dimensions; ++d) {
        if (trie_dimension_is_used(node, d)) continue;
        const int fanout = trie_dimension_fanout(trie, index, d);
        if (fanout <= 1) continue;
        int counts[TRIE_MAX_FANOUT] = {0};
        for (int i = 0; i < node->size; ++i) ++counts[trie_bucket(node->words[i], d, fanout)];
        double entropy = 0.0;
        int occupied = 0;
        for (int b = 0; b < fanout; ++b) if (counts[b] != 0) {
            double p = (double) counts[b] / node->size;
            entropy -= p * log(p); ++occupied;
        }
        if (occupied < 2) continue;
        double score = entropy / log((double) fanout);
        if (variances != NULL) score *= variances[d];
        if (score > best_score) { best_score = score; best = d; }
    }
    return best;
}

static int trie_split_leaf(struct symbolic_trie_index *trie, isax_index *index, symbolic_trie_node *node) {
    int dimension = trie_choose_split(node, trie->split_dimensions, trie, index,
                                      index->settings->symbolic_variances);
    if (dimension < 0) { node->split_exhausted = 1; ++trie->split_exhausted_leaves; return 0; }
    trie_dimension_mask used = trie_dimension_mark_used(node->used_dimensions, dimension);
    const int fanout = trie_dimension_fanout(trie, index, dimension);
    symbolic_trie_node *children[TRIE_MAX_FANOUT] = {0};
    for (int i = 0; i < node->size; ++i) {
        int bucket = trie_bucket(node->words[i], dimension, fanout);
        if (children[bucket] == NULL && (children[bucket] = trie_node_create(trie->dimensions, node, used)) == NULL) return 0;
        if (!trie_leaf_append(children[bucket], node->words[i], node->positions[i], trie->dimensions)) return 0;
        node->words[i] = NULL;
    }
    free(node->words); free(node->positions);
    node->words = NULL; node->positions = NULL; node->size = node->capacity = 0;
    node->leaf = 0; node->split_dimension = dimension; node->split_fanout = fanout;
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
    const int split_dimensions = trie->split_dimensions;
    const int workers = maxquerythread > 0 ? maxquerythread : 1;
    unsigned long *histograms = calloc((size_t) workers * split_dimensions * TRIE_MAX_FANOUT,
                                       sizeof(*histograms));
    if (histograms == NULL) return -1;

#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
    {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        unsigned long *local = histograms + (size_t) worker * split_dimensions * TRIE_MAX_FANOUT;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (long sample = 0; sample < sample_size; ++sample) {
            long record = (long) (((unsigned long long) sample * node->size) / sample_size);
            const sax_type *word = node->words[record];
            for (int d = 0; d < split_dimensions; ++d) {
                const int fanout = trie_dimension_fanout(trie, index, d);
                if (fanout > 1) ++local[d * TRIE_MAX_FANOUT + trie_bucket(word, d, fanout)];
            }
        }
    }

    int dimension = -1;
    double best_score = -1.0;
    for (int d = 0; d < split_dimensions; ++d) {
        if (trie_dimension_is_used(node, d)) continue;
        const int fanout = trie_dimension_fanout(trie, index, d);
        if (fanout <= 1) continue;
        unsigned long counts[TRIE_MAX_FANOUT] = {0};
        for (int worker = 0; worker < workers; ++worker) {
            const unsigned long *local = histograms + (size_t) worker * split_dimensions * TRIE_MAX_FANOUT + d * TRIE_MAX_FANOUT;
            for (int bucket = 0; bucket < fanout; ++bucket) counts[bucket] += local[bucket];
        }
        double entropy = 0.0;
        int occupied = 0;
        for (int bucket = 0; bucket < fanout; ++bucket) if (counts[bucket] != 0) {
            double p = (double) counts[bucket] / sample_size;
            entropy -= p * log(p);
            ++occupied;
        }
        if (occupied < 2) continue;
        double score = entropy / log((double) fanout);
        if (index->settings->symbolic_variances != NULL) score *= index->settings->symbolic_variances[d];
        if (score > best_score) { best_score = score; dimension = d; }
    }
    free(histograms);
    if (dimension < 0) {
        node->split_exhausted = 1;
        ++trie->split_exhausted_leaves;
        return 0;
    }

    const int fanout = trie_dimension_fanout(trie, index, dimension);
    unsigned long *counts = calloc((size_t) workers * fanout, sizeof(*counts));
    unsigned long *offsets = calloc((size_t) workers * fanout, sizeof(*offsets));
    sax_type *local_min = malloc((size_t) workers * fanout * dimensions);
    sax_type *local_max = calloc((size_t) workers * fanout * dimensions, sizeof(*local_max));
    if (counts == NULL || offsets == NULL || local_min == NULL || local_max == NULL) {
        free(counts); free(offsets); free(local_min); free(local_max);
        return -1;
    }
    memset(local_min, UCHAR_MAX, (size_t) workers * fanout * dimensions);

#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
    {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        unsigned long *local = counts + (size_t) worker * fanout;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int i = 0; i < node->size; ++i)
            ++local[trie_bucket(node->words[i], dimension, fanout)];
    }

    unsigned long bucket_sizes[TRIE_MAX_FANOUT] = {0};
    for (int bucket = 0; bucket < fanout; ++bucket) {
        unsigned long offset = 0;
        for (int worker = 0; worker < workers; ++worker) {
            offsets[(size_t) worker * fanout + bucket] = offset;
            offset += counts[(size_t) worker * fanout + bucket];
        }
        bucket_sizes[bucket] = offset;
    }

    symbolic_trie_node *children[TRIE_MAX_FANOUT] = {0};
    const trie_dimension_mask used = trie_dimension_mark_used(node->used_dimensions, dimension);
    for (int bucket = 0; bucket < fanout; ++bucket) if (bucket_sizes[bucket] != 0) {
        children[bucket] = trie_node_create(dimensions, node, used);
        if (children[bucket] == NULL ||
            (children[bucket]->words = malloc(sizeof(*children[bucket]->words) * bucket_sizes[bucket])) == NULL ||
            (children[bucket]->positions = malloc(sizeof(*children[bucket]->positions) * bucket_sizes[bucket])) == NULL) {
            for (int b = 0; b < fanout; ++b) if (children[b] != NULL) {
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
        unsigned long *local_offsets = offsets + (size_t) worker * fanout;
        sax_type *min_word = local_min + (size_t) worker * fanout * dimensions;
        sax_type *max_word = local_max + (size_t) worker * fanout * dimensions;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int i = 0; i < node->size; ++i) {
            sax_type *word = node->words[i];
            int bucket = trie_bucket(word, dimension, fanout);
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

    for (int bucket = 0; bucket < fanout; ++bucket) if (children[bucket] != NULL) {
        for (int d = 0; d < dimensions; ++d) {
            sax_type min = UCHAR_MAX, max = 0;
            for (int worker = 0; worker < workers; ++worker) {
                const sax_type local_low = local_min[((size_t) worker * fanout + bucket) * dimensions + d];
                const sax_type local_high = local_max[((size_t) worker * fanout + bucket) * dimensions + d];
                if (local_low < min) min = local_low;
                if (local_high > max) max = local_high;
            }
            children[bucket]->min_word[d] = min;
            children[bucket]->max_word[d] = max;
        }
        children[bucket]->mbb_valid = 1;
    }
    trie_node_merge_child_mbbs(node, children, fanout, dimensions);
    free(counts); free(offsets); free(local_min); free(local_max);
    free(node->words); free(node->positions);
    node->words = NULL; node->positions = NULL; node->size = node->capacity = 0;
    node->leaf = 0; node->split_dimension = dimension; node->split_fanout = fanout;
    memcpy(node->children, children, sizeof(children));
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

#if HAVE_CBLAS
/* Project one bounded bulk-build block.  The BLAS implementation owns the
 * matrix-multiply threads; OpenMP is used only for independent work around
 * that call. */
static enum response trie_build_pca_block(isax_index *index,
                                          struct symbolic_trie_index *trie,
                                          long first_record, unsigned int records,
                                          int apply_znorm) {
    const int dimensions = trie->dimensions;
    const int ts_length = index->settings->timeseries_size;
    const int function_type = index->settings->function_type;
    if (records == 0 || (function_type != 5 && function_type != 6) ||
        dimensions <= 0 || ts_length <= 0 || index->pca_dim <= 0 ||
        (size_t) records > SIZE_MAX / (size_t) dimensions) return FAILURE;

    ts_type *projection = malloc(sizeof(*projection) * (size_t) records * dimensions);
    if (projection == NULL) return FAILURE;
    const ts_type *input = rawfile + (size_t) first_record * ts_length;
    ts_type *fft_input = NULL;
    int failed = 0;

    if (function_type == 5) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(maxquerythread)
#endif
        for (unsigned int row = 0; row < records; ++row) {
            if (apply_znorm)
                znorm(rawfile + (size_t) (first_record + row) * ts_length, ts_length);
        }
    } else {
        if ((size_t) records > SIZE_MAX / (size_t) index->pca_dim) {
            free(projection);
            return FAILURE;
        }
        fft_input = malloc(sizeof(*fft_input) * (size_t) records * index->pca_dim);
        if (fft_input == NULL) {
            free(projection);
            return FAILURE;
        }
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#endif
        {
            fftw_workspace fftw = {0};
            pthread_mutex_lock(&trie_fftw_plan_lock);
            fftw_workspace_init(&fftw, ts_length);
            pthread_mutex_unlock(&trie_fftw_plan_lock);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (unsigned int row = 0; row < records; ++row) {
                if (failed) continue;
                ts_type *ts = rawfile + (size_t) (first_record + row) * ts_length;
                if (apply_znorm) znorm(ts, ts_length);
                memcpy(fftw.ts, ts, sizeof(*ts) * ts_length);
                if (fft_full_real_from_ts(index, &fftw) != SUCCESS) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                    failed = 1;
                } else {
                    memcpy(fft_input + (size_t) row * index->pca_dim, fftw.transform,
                           sizeof(*fft_input) * index->pca_dim);
                }
            }
            pthread_mutex_lock(&trie_fftw_plan_lock);
            fftw_workspace_destroy(&fftw);
            pthread_mutex_unlock(&trie_fftw_plan_lock);
        }
        input = fft_input;
    }

    if (!failed && pca_project_batch(index, input, records, projection,
                                     maxquerythread > 0 ? maxquerythread : 1) != SUCCESS) failed = 1;
    if (!failed) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(maxquerythread)
#endif
        for (unsigned int row = 0; row < records; ++row) {
            const long record = first_record + row;
            sax_type *word = trie->word_arena + (size_t) record * trie->dimensions;
            const ts_type *values = projection + (size_t) row * dimensions;
            if (function_type == 5) spartan_from_pca(index, values, word);
            else sfa_from_fft(index, values, word);
            trie->root->words[record] = word;
            trie->root->positions[record] = (file_position_type) ((size_t) record * ts_length);
        }
    }

    free(fft_input);
    free(projection);
    return failed ? FAILURE : SUCCESS;
}
#endif

static int trie_insert(struct symbolic_trie_index *trie, isax_index *index, sax_type *word,
                       file_position_type position) {
    symbolic_trie_node *node = trie->root;
    while (!node->leaf) {
        trie_node_update_mbb(node, word, trie->dimensions);
        int bucket = trie_bucket(word, node->split_dimension, node->split_fanout);
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
    for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL) {
#ifdef _OPENMP
#pragma omp task firstprivate(i) if (depth < 6)
#endif
        trie_build_subtree(trie, index, node->children[i], depth + 1);
    }
}

enum response symbolic_trie_build(isax_index *index, const char *path, long ts_num,
                                  int filetype_int, int apply_znorm) {
    if (index == NULL || index->settings == NULL || ts_num <= 0) return FAILURE;
    double build_start = messi_monotonic_seconds();
    FILE *file = fopen(path, "rb");
    if (file == NULL) return FAILURE;
    if (fseek(file, (long) index->settings->dataset_header_bytes, SEEK_SET) != 0) {
        fclose(file);
        return FAILURE;
    }
    messi_build_progress build_progress;
    messi_build_progress_init(&build_progress);
    const size_t values = (size_t) ts_num * index->settings->timeseries_size;
    rawfile = malloc(values * sizeof(*rawfile));
    if (rawfile == NULL) { fclose(file); messi_build_progress_abort(&build_progress); return FAILURE; }
    if (filetype_int) {
        /* Convert in bounded chunks so integer input does not require a
         * second dataset-sized allocation beside rawfile. */
        const size_t chunk_values = 1U << 20;
        file_type *input = malloc(chunk_values * sizeof(*input));
        if (input == NULL) { fclose(file); free(rawfile); rawfile = NULL; messi_build_progress_abort(&build_progress); return FAILURE; }
        for (size_t offset = 0; offset < values; offset += chunk_values) {
            size_t count = values - offset < chunk_values ? values - offset : chunk_values;
            if (fread(input, sizeof(*input), count, file) != count) {
                free(input); fclose(file); free(rawfile); rawfile = NULL; messi_build_progress_abort(&build_progress); return FAILURE;
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(maxquerythread)
#endif
            for (size_t i = 0; i < count; ++i)
                rawfile[offset + i] = file_value_to_ts(input[i], filetype_int);
            messi_build_progress_update(&build_progress,
                10.0 * (double) (offset + count) / (double) values);
        }
        free(input);
    } else if (fread(rawfile, sizeof(*rawfile), values, file) != values) {
        fclose(file); free(rawfile); rawfile = NULL; messi_build_progress_abort(&build_progress); return FAILURE;
    } else {
        messi_build_progress_update(&build_progress, 10.0);
    }
    fclose(file);
    double read_end = messi_monotonic_seconds();
    struct symbolic_trie_index *trie = calloc(1, sizeof(*trie));
    if (trie == NULL ||
        (trie->root = trie_node_create(index->settings->n_segments, NULL,
                                       (trie_dimension_mask) { 0, 0 })) == NULL) {
        free(trie); free(rawfile); rawfile = NULL; messi_build_progress_abort(&build_progress); return FAILURE;
    }
    trie->dimensions = index->settings->n_segments;
    trie->split_dimensions = index->settings->trie_split_dimensions > 0
                                 ? index->settings->trie_split_dimensions
                                 : trie->dimensions;
    trie->dynamic_alphabet = index->settings->trie_dynamic_alphabet &&
                             index->settings->root_bit_cardinalities != NULL;
    trie->fanout = trie->dynamic_alphabet
                       ? (1 << index->settings->trie_max_bits)
                       : index->settings->trie_fanout;
    trie->bound_dimensions = index->settings->trie_bound_dimensions > 0
                                 ? index->settings->trie_bound_dimensions
                                 : trie->dimensions;
    /* Materialize all full words in a single arena.  Subsequent splits move
     * only pointers and positions, never raw series or word ownership. */
    trie->root->words = calloc((size_t) ts_num, sizeof(*trie->root->words));
    trie->root->positions = malloc((size_t) ts_num * sizeof(*trie->root->positions));
    if ((size_t) ts_num > SIZE_MAX / (size_t) trie->dimensions) {
        trie_node_destroy(trie->root); free(trie); free(rawfile); rawfile = NULL; messi_build_progress_abort(&build_progress); return FAILURE;
    }
    trie->word_arena = malloc((size_t) ts_num * (size_t) trie->dimensions * sizeof(*trie->word_arena));
    trie->root->capacity = trie->root->size = (int) ts_num;
    int failed = trie->root->words == NULL || trie->root->positions == NULL || trie->word_arena == NULL;
    unsigned long transform_completed = 0;
#if HAVE_CBLAS
    if ((index->settings->function_type == 5 || index->settings->function_type == 6) && !failed) {
        int reported_transform_percent = 10;
        for (long first = 0; first < ts_num && !failed;
             first += TRIE_PCA_PROJECTION_BLOCK_RECORDS) {
            const unsigned int records = (unsigned int) ((ts_num - first) <
                (long) TRIE_PCA_PROJECTION_BLOCK_RECORDS
                    ? (ts_num - first) : TRIE_PCA_PROJECTION_BLOCK_RECORDS);
            if (trie_build_pca_block(index, trie, first, records, apply_znorm) != SUCCESS) {
                failed = 1;
            } else {
                double percent = 10.0 +
                    45.0 * (double) (first + records) / (double) ts_num;
                int whole_percent = (int) percent;
                if (whole_percent > reported_transform_percent) {
                    messi_build_progress_update(&build_progress, percent);
                    reported_transform_percent = whole_percent;
                }
            }
        }
    } else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#endif
    {
        fftw_workspace fftw = {0};
        ts_type *transform = malloc(sizeof(*transform) * (size_t) trie->dimensions);
        unsigned long transform_pending = 0;
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
            sax_type *word = trie->word_arena + (size_t) i * trie->dimensions;
            if (transform == NULL || trie_word_from_ts(index, ts, word, transform, &fftw) != SUCCESS) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                failed = 1;
            } else {
                trie->root->words[i] = word;
                trie->root->positions[i] = (file_position_type) ((size_t) i * index->settings->timeseries_size);
            }
            if (++transform_pending == 1024 || i + 1 == ts_num) {
                unsigned long completed = __sync_add_and_fetch(&transform_completed,
                                                                transform_pending);
                messi_build_progress_update(&build_progress,
                    10.0 + 45.0 * (double) completed / (double) ts_num);
                transform_pending = 0;
            }
        }
        free(transform);
        pthread_mutex_lock(&trie_fftw_plan_lock);
        fftw_workspace_destroy(&fftw);
        pthread_mutex_unlock(&trie_fftw_plan_lock);
    }
    }
    if (failed) { trie_node_destroy(trie->root); free(trie->word_arena); free(trie); free(rawfile); rawfile = NULL; messi_build_progress_abort(&build_progress); return FAILURE; }
    double transform_end = messi_monotonic_seconds();
    messi_build_progress_update(&build_progress, 55.0);

    /* Split the initially enormous root before entering the task region.  The
     * sampled splitter uses OpenMP worksharing itself; nested worksharing from
     * inside the later omp single region would otherwise be serialized. */
    if (trie->root->size > index->settings->max_leaf_size) {
        int root_split = trie_split_root_sampled_parallel(trie, index, trie->root);
        if (root_split < 0) {
            trie_node_destroy(trie->root); free(trie->word_arena); free(trie); free(rawfile); rawfile = NULL;
            messi_build_progress_abort(&build_progress);
            return FAILURE;
        }
        if (root_split == 0) {
            for (long i = 0; i < ts_num; ++i)
                trie_node_update_mbb(trie->root, trie->root->words[i], trie->dimensions);
        }
    } else {
        for (long i = 0; i < ts_num; ++i)
            trie_node_update_mbb(trie->root, trie->root->words[i], trie->dimensions);
    }
    messi_build_progress_update(&build_progress, 65.0);
#ifdef _OPENMP
#pragma omp parallel num_threads(maxquerythread)
#pragma omp single
#endif
    {
        if (trie->root->leaf) {
            trie_build_subtree(trie, index, trie->root, 0);
        } else {
            for (int i = 0; i < trie->root->split_fanout; ++i) if (trie->root->children[i] != NULL) {
#ifdef _OPENMP
#pragma omp task firstprivate(i)
#endif
                trie_build_subtree(trie, index, trie->root->children[i], 1);
            }
        }
    }
    messi_build_progress_update(&build_progress, 90.0);
    if (index->settings->trie_leaf_ivf != 0) {
        const double clustering_start = messi_monotonic_seconds();
        float *centers = trie_build_bin_centers(index, trie->dimensions);
        int clustering_workers = 1;
        size_t eligible_leaves = 0;
        if (centers == NULL || !trie_cluster_leaves_parallel(trie, index, centers,
                                                              &clustering_workers, &eligible_leaves)) {
            fprintf(stderr, "warning: unable to build some trie leaf IVF directories; affected leaves use plain scans.\n");
        }
        free(centers);
        fprintf(stderr, ">>> trie leaf IVF timing\n");
        fprintf(stderr, "    eligible leaves  : %zu\n", eligible_leaves);
        fprintf(stderr, "    workers          : %d\n", clustering_workers);
        fprintf(stderr, "    clustered leaves : %lu\n", trie->clustered_leaves);
        fprintf(stderr, "    clusters         : %llu\n", trie->cluster_count);
        if (index->settings->trie_leaf_ivf_radial_bound) {
            char radial_count[32], radial_bytes[32];
            fprintf(stderr, "    radial leaves    : %lu\n", trie->radial_leaves);
            fprintf(stderr, "    radial records   : %s\n",
                    trie_format_compact_count(trie->radial_records, radial_count));
            fprintf(stderr, "    radius storage   : %sB\n",
                    trie_format_compact_count(trie->radial_records * sizeof(float), radial_bytes));
            if (trie->radial_leaves != trie->clustered_leaves)
                fprintf(stderr,
                        "warning: record-radius allocation failed for %lu clustered leaves; those leaves use ordinary IVF scans.\n",
                        trie->clustered_leaves - trie->radial_leaves);
        }
        fprintf(stderr, "    total            : %.3f s\n", messi_monotonic_seconds() - clustering_start);
    }
    trie->node_count = trie_node_count(trie->root);
    unsigned long internal_nodes = 0;
    unsigned long long nonempty_children = 0;
    unsigned long long child_slots = 0;
    double largest_child_share = 0.0;
    trie_collect_split_diagnostics(trie->root, &internal_nodes, &nonempty_children,
                                   &child_slots,
                                   &largest_child_share);
    double split_end = messi_monotonic_seconds();
    index->trie = trie;
    index->total_records = ts_num;
    messi_build_progress_finish(&build_progress);
    fprintf(stderr, ">>> trie build timing\n");
    fprintf(stderr, "    read       : %.3f s\n", read_end - build_start);
    fprintf(stderr, "    transform  : %.3f s\n", transform_end - read_end);
    fprintf(stderr, "    split      : %.3f s\n", split_end - transform_end);
    fprintf(stderr, "    total      : %.3f s\n", split_end - build_start);
    if (internal_nodes != 0) {
        char internal_count[32];
        fprintf(stderr, ">>> trie split diagnostics\n");
        fprintf(stderr, "    internal nodes : %s\n",
                trie_format_compact_count((unsigned long long) internal_nodes,
                                          internal_count));
        fprintf(stderr, "    max fanout     : %d\n", trie->fanout);
        fprintf(stderr, "    avg non-empty  : %.2f\n",
                (double) nonempty_children / (double) internal_nodes);
        fprintf(stderr, "    occupancy      : %.1f%%\n",
                child_slots == 0 ? 0.0 :
                    100.0 * (double) nonempty_children / (double) child_slots);
        fprintf(stderr, "    largest child  : %.1f%%\n",
                100.0 * largest_child_share / (double) internal_nodes);
    }
    return SUCCESS;
}

static float trie_refine_streaming_record(isax_index *index, const symbolic_trie_node *node,
                                          int record, float mbr_suffix,
                                          const ts_type *query, const ts_type *transform,
                                          float bsf, trie_query_stats *stats,
                                          trie_query_scratch *scratch,
                                          file_position_type *best_position) {
    if (stats != NULL) ++stats->lower_bounds;
    const unsigned long long bound_start = profile_query_phases && stats != NULL
                                               ? trie_monotonic_microseconds() : 0;
    const float lower = trie_record_lower_bound(index->trie, index, transform,
                                                node->words[record], bsf,
                                                mbr_suffix, scratch);
    if (bound_start != 0)
        stats->record_bound_microseconds += trie_monotonic_microseconds() - bound_start;
    if (lower > bsf) {
        return bsf;
    }
    if (stats != NULL) ++stats->exact_distances;
    const float distance = trie_exact_distance(query, rawfile + node->positions[record],
                                               index->settings->timeseries_size, bsf, stats);
    if (distance < bsf) {
        if (best_position != NULL) *best_position = node->positions[record];
        return distance;
    }
    return bsf;
}

static inline int trie_radial_simd_lanes(const isax_index *index) {
    if (!index->settings->SIMD_flag) return 1;
#if defined(__AVX512F__)
    return 16;
#elif ADS_HAVE_AVX2
    return 8;
#else
    return 1;
#endif
}

/* Return one survivor bit per contiguous float radius.  Using a prune mask
 * makes unordered (NaN) radii survive, matching the conservative scalar
 * predicate. */
static inline unsigned int trie_radial_survivor_mask(const float *radii,
                                                      messi_radial_window window,
                                                      int lanes) {
#if defined(__AVX512F__)
    if (lanes == 16) {
        const __m512 values = _mm512_loadu_ps(radii);
        const __mmask16 below = _mm512_cmp_ps_mask(
            values, _mm512_set1_ps(window.lower), _CMP_LT_OQ);
        const __mmask16 above = _mm512_cmp_ps_mask(
            values, _mm512_set1_ps(window.upper), _CMP_GT_OQ);
        return (unsigned int) ((~(below | above)) & 0xffffU);
    }
#elif ADS_HAVE_AVX2
    if (lanes == 8) {
        const __m256 values = _mm256_loadu_ps(radii);
        const __m256 below = _mm256_cmp_ps(
            values, _mm256_set1_ps(window.lower), _CMP_LT_OQ);
        const __m256 above = _mm256_cmp_ps(
            values, _mm256_set1_ps(window.upper), _CMP_GT_OQ);
        const unsigned int pruned = (unsigned int) _mm256_movemask_ps(
            _mm256_or_ps(below, above));
        return (~pruned) & 0xffU;
    }
#endif
    return messi_radial_radius_in_window(radii[0], window) ? 1U : 0U;
}

static float trie_scan_leaf_range(isax_index *index, const symbolic_trie_node *node,
                                  int offset, int count, float mbr_suffix,
                                  const float *record_radii,
                                  messi_distance_interval query_radius,
                                  const ts_type *query, const ts_type *transform, float bsf,
                                  trie_query_stats *stats, trie_query_scratch *scratch,
                                  file_position_type *best_position) {
    const int end = offset + count;
    const int streaming = index->settings->trie_streaming_leaf_scan ||
                          !trie_scratch_reserve(scratch, count);
    int candidate_count = 0;
    const unsigned long lower_bounds_before = stats != NULL ? stats->lower_bounds : 0;
    const unsigned long long bound_start = !streaming && profile_query_phases && stats != NULL
                                               ? trie_monotonic_microseconds() : 0;
    const int radial_applied = record_radii != NULL &&
        (!index->settings->trie_leaf_ivf_radial_bound_auto ||
         trie_radial_auto_gate_allows(scratch));
    const int radial_sampling = radial_applied &&
        index->settings->trie_leaf_ivf_radial_bound_auto &&
        scratch != NULL && scratch->radial_auto_gate != NULL &&
        __atomic_load_n(&scratch->radial_auto_gate->decision,
                        __ATOMIC_ACQUIRE) == 0;
    const int radial_lanes = radial_applied ? trie_radial_simd_lanes(index) : 1;
    messi_radial_window radial_window = { -INFINITY, INFINITY };
    float radial_window_bsf = 0.0f;
    int radial_window_valid = 0;
    unsigned long sampled_survivors = 0;

    if (!radial_applied) {
        for (int record = offset; record < end; ++record) {
            if (streaming) {
                bsf = trie_refine_streaming_record(index, node, record, mbr_suffix,
                                                   query, transform, bsf, stats, scratch,
                                                   best_position);
                continue;
            }
            if (stats != NULL) ++stats->lower_bounds;
            const float lower = trie_record_lower_bound(index->trie, index, transform,
                                                        node->words[record], bsf,
                                                        mbr_suffix, scratch);
            if (lower <= bsf) {
                scratch->candidates[candidate_count].lower_bound = lower;
                scratch->candidates[candidate_count++].record_index = record;
            }
        }
    } else {
        if (stats != NULL) stats->radial_candidates += (unsigned long) count;
        int block = offset;
        while (block < end) {
            if (!radial_window_valid || radial_window_bsf != bsf) {
                radial_window = messi_radial_window_from_bsf(query_radius, bsf);
                radial_window_bsf = bsf;
                radial_window_valid = 1;
            }
            int lanes = 1;
            if (radial_lanes > 1 && block + radial_lanes <= end)
                lanes = radial_lanes;
            unsigned int survivors = trie_radial_survivor_mask(
                record_radii + block, radial_window, lanes);
            if (radial_sampling && stats == NULL)
                sampled_survivors += (unsigned long) __builtin_popcount(survivors);
            while (survivors != 0U) {
                const int lane = __builtin_ctz(survivors);
                const int record = block + lane;
                survivors &= survivors - 1U;
                if (streaming) {
                    bsf = trie_refine_streaming_record(index, node, record, mbr_suffix,
                                                       query, transform, bsf, stats, scratch,
                                                       best_position);
                    continue;
                }
                if (stats != NULL) ++stats->lower_bounds;
                const float lower = trie_record_lower_bound(index->trie, index, transform,
                                                            node->words[record], bsf,
                                                            mbr_suffix, scratch);
                if (lower <= bsf) {
                    scratch->candidates[candidate_count].lower_bound = lower;
                    scratch->candidates[candidate_count++].record_index = record;
                }
            }
            block += lanes;
        }
    }
    unsigned long radial_pruned = 0;
    if (radial_applied) {
        radial_pruned = stats != NULL
            ? (unsigned long) count - (stats->lower_bounds - lower_bounds_before)
            : (radial_sampling ? (unsigned long) count - sampled_survivors : 0);
        if (stats != NULL) stats->radial_pruned += radial_pruned;
        if (radial_sampling)
            trie_radial_auto_gate_observe(scratch->radial_auto_gate,
                                          (unsigned long) count, radial_pruned);
    }
    if (streaming) return bsf;
    if (bound_start != 0)
        stats->record_bound_microseconds += trie_monotonic_microseconds() - bound_start;

    unsigned long long heap_start = profile_query_phases && stats != NULL
                                        ? trie_monotonic_microseconds() : 0;
    trie_heap_build(scratch->candidates, candidate_count);
    if (heap_start != 0)
        stats->candidate_heap_microseconds += trie_monotonic_microseconds() - heap_start;
    while (candidate_count > 0 && scratch->candidates[0].lower_bound <= bsf) {
        heap_start = profile_query_phases && stats != NULL ? trie_monotonic_microseconds() : 0;
        const trie_leaf_candidate candidate = trie_heap_pop(scratch->candidates, &candidate_count);
        if (heap_start != 0)
            stats->candidate_heap_microseconds += trie_monotonic_microseconds() - heap_start;
        if (stats != NULL) ++stats->exact_distances;
        const float distance = trie_exact_distance(query,
            rawfile + node->positions[candidate.record_index],
            index->settings->timeseries_size, bsf, stats);
        if (distance < bsf) {
            bsf = distance;
            if (best_position != NULL) *best_position = node->positions[candidate.record_index];
        }
    }
    return bsf;
}

static float trie_scan_leaf_best_first(isax_index *index, const symbolic_trie_node *node,
                                       const ts_type *query, const ts_type *transform, float bsf,
                                       trie_query_stats *stats, trie_query_scratch *scratch,
                                       file_position_type *best_position) {
    if (node->cluster_count == 0) {
        const float suffix = trie_record_mbr_suffix(index->trie, index, transform,
                                                    node->min_word, node->max_word, scratch);
        return trie_scan_leaf_range(index, node, 0, node->size, suffix, NULL,
                                    (messi_distance_interval) { 0.0, 0.0 },
                                    query, transform, bsf, stats, scratch, best_position);
    }
    typedef struct {
        int cluster;
        float bound;
        messi_distance_interval query_radius;
    } trie_cluster_bound;
    trie_cluster_bound ordered[64];
    int count = 0;
    for (int cluster = 0; cluster < node->cluster_count; ++cluster) {
        const trie_leaf_cluster *group = &node->clusters[cluster];
        if (stats != NULL) ++stats->cluster_bounds;
        const float symbolic_bound = trie_lower_bound(index->trie, index, transform,
            node->cluster_min_words + (size_t) cluster * index->trie->dimensions,
            node->cluster_max_words + (size_t) cluster * index->trie->dimensions,
            index->trie->dimensions, bsf, stats);
        if (symbolic_bound >= bsf) {
            if (stats != NULL) {
                ++stats->cluster_symbolic_pruned;
                stats->cluster_symbolic_records_pruned += (unsigned long) group->size;
                ++stats->cluster_pruned;
                stats->cluster_records_pruned += (unsigned long) group->size;
            }
            continue;
        }
        messi_distance_interval query_radius = { 0.0, 0.0 };
        const int radial_may_run = node->record_raw_radii != NULL &&
            (!index->settings->trie_leaf_ivf_radial_bound_auto ||
             trie_radial_auto_gate_allows(scratch));
        const int needs_query_radius = index->settings->trie_leaf_ivf_raw_ball_bound ||
                                       radial_may_run;
        if (needs_query_radius)
            query_radius = trie_cluster_query_radius(node, cluster, query,
                                                     index->settings->timeseries_size);
        float raw_ball_bound = 0.0f;
        if (index->settings->trie_leaf_ivf_raw_ball_bound) {
            if (stats != NULL) ++stats->cluster_raw_ball_bounds;
            raw_ball_bound = trie_cluster_raw_ball_lower_bound(node, cluster, query_radius);
            if (raw_ball_bound >= bsf) {
                if (stats != NULL) {
                    ++stats->cluster_raw_ball_pruned;
                    stats->cluster_raw_ball_records_pruned += (unsigned long) group->size;
                    ++stats->cluster_pruned;
                    stats->cluster_records_pruned += (unsigned long) group->size;
                }
                continue;
            }
        }
        const float bound = symbolic_bound > raw_ball_bound ? symbolic_bound : raw_ball_bound;
        ordered[count++] = (trie_cluster_bound) { cluster, bound, query_radius };
    }
    for (int i = 1; i < count; ++i) {
        trie_cluster_bound current = ordered[i]; int j = i - 1;
        while (j >= 0 && ordered[j].bound > current.bound) { ordered[j + 1] = ordered[j]; --j; }
        ordered[j + 1] = current;
    }
    for (int i = 0; i < count; ++i) {
        const int cluster = ordered[i].cluster;
        const trie_leaf_cluster *group = &node->clusters[cluster];
        const sax_type *minimum = node->cluster_min_words + (size_t) cluster * index->trie->dimensions;
        const sax_type *maximum = node->cluster_max_words + (size_t) cluster * index->trie->dimensions;
        const float suffix = trie_record_mbr_suffix(index->trie, index, transform,
                                                    minimum, maximum, scratch);
        bsf = trie_scan_leaf_range(index, node, group->offset, group->size, suffix,
                                   node->record_raw_radii, ordered[i].query_radius,
                                   query, transform, bsf, stats, scratch, best_position);
    }
    return bsf;
}

static float trie_search_node(isax_index *index, const symbolic_trie_node *node, const ts_type *query,
                              const ts_type *transform, float bsf, trie_query_stats *stats,
                              const symbolic_trie_node *skip_leaf, trie_query_scratch *scratch,
                              file_position_type *best_position) {
    if (node == skip_leaf) return bsf;
    if (stats != NULL) ++stats->checked_nodes;
    float lower = trie_lower_bound(index->trie, index, transform, node->min_word, node->max_word,
                                   index->trie->dimensions, bsf, stats);
    if (lower > bsf) return bsf;
    if (!node->leaf) {
        typedef struct { const symbolic_trie_node *node; float bound; } trie_child_bound;
        trie_child_bound ordered[TRIE_MAX_FANOUT];
        int n = 0;
        for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL) {
            ordered[n].node = node->children[i];
            ordered[n++].bound = trie_lower_bound(index->trie, index, transform,
                node->children[i]->min_word, node->children[i]->max_word,
                index->trie->dimensions, bsf, stats);
        }
        for (int i = 1; i < n; ++i) {
            trie_child_bound current = ordered[i]; int j = i - 1;
            while (j >= 0 && ordered[j].bound > current.bound) { ordered[j + 1] = ordered[j]; --j; }
            ordered[j + 1] = current;
        }
        for (int i = 0; i < n; ++i) if (ordered[i].bound <= bsf)
            bsf = trie_search_node(index, ordered[i].node, query, transform, bsf, stats, skip_leaf, scratch,
                                   best_position);
        return bsf;
    }
    return trie_scan_leaf_best_first(index, node, query, transform, bsf, stats, scratch, best_position);
}

static const symbolic_trie_node *trie_seed_leaf(isax_index *index, const sax_type *query_word,
                                                 const ts_type *transform, trie_query_stats *stats) {
    const symbolic_trie_node *node = index->trie->root;
    while (!node->leaf) {
        int bucket = trie_bucket(query_word, node->split_dimension, node->split_fanout);
        const symbolic_trie_node *next = node->children[bucket];
        if (next == NULL) {
            float best_bound = FLT_MAX;
            for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL) {
                float bound = trie_lower_bound(index->trie, index, transform,
                                                node->children[i]->min_word,
                                                node->children[i]->max_word,
                                                index->trie->dimensions, best_bound, stats);
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

/* The single-query engine deliberately mirrors the iSAX MESSI engine: first
 * distribute independent top-level work, then process low-bound leaves from
 * several locked priority queues.  A trie has one physical root, so it first
 * expands that root into a small frontier; thereafter its workers have the
 * same work-unit granularity as iSAX's many first-buffer-layer roots. */
typedef struct {
    const symbolic_trie_node *node;
    float lower_bound;
} trie_leaf_work;

typedef struct {
    trie_leaf_work *items;
    int size;
    int capacity;
    pthread_mutex_t lock;
} trie_leaf_queue;

static void trie_leaf_queue_sift_down(trie_leaf_work *items, int size, int position) {
    while (1) {
        int left = 2 * position + 1, right = left + 1, smallest = position;
        if (left < size && items[left].lower_bound < items[smallest].lower_bound) smallest = left;
        if (right < size && items[right].lower_bound < items[smallest].lower_bound) smallest = right;
        if (smallest == position) return;
        trie_leaf_work tmp = items[position]; items[position] = items[smallest]; items[smallest] = tmp;
        position = smallest;
    }
}

static int trie_leaf_queue_push(trie_leaf_queue *queue, trie_leaf_work item,
                                trie_query_stats *stats) {
    unsigned long long start = profile_query_phases && stats != NULL
                                   ? trie_monotonic_microseconds() : 0;
    int ok = 1;
    pthread_mutex_lock(&queue->lock);
    if (queue->size == queue->capacity) {
        int capacity = queue->capacity == 0 ? 128 : queue->capacity * 2;
        trie_leaf_work *items = realloc(queue->items, (size_t) capacity * sizeof(*items));
        if (items == NULL) ok = 0;
        else { queue->items = items; queue->capacity = capacity; }
    }
    if (ok) {
        int position = queue->size++;
        queue->items[position] = item;
        while (position > 0) {
            int parent = (position - 1) / 2;
            if (queue->items[parent].lower_bound <= queue->items[position].lower_bound) break;
            trie_leaf_work tmp = queue->items[parent]; queue->items[parent] = queue->items[position];
            queue->items[position] = tmp; position = parent;
        }
    }
    pthread_mutex_unlock(&queue->lock);
    if (start != 0) stats->queue_microseconds += trie_monotonic_microseconds() - start;
    return ok;
}

static int trie_leaf_queue_pop(trie_leaf_queue *queue, trie_leaf_work *item,
                               trie_query_stats *stats) {
    unsigned long long start = profile_query_phases && stats != NULL
                                   ? trie_monotonic_microseconds() : 0;
    pthread_mutex_lock(&queue->lock);
    if (queue->size == 0) {
        pthread_mutex_unlock(&queue->lock);
        if (start != 0) stats->queue_microseconds += trie_monotonic_microseconds() - start;
        return 0;
    }
    *item = queue->items[0];
    queue->items[0] = queue->items[--queue->size];
    trie_leaf_queue_sift_down(queue->items, queue->size, 0);
    pthread_mutex_unlock(&queue->lock);
    if (start != 0) stats->queue_microseconds += trie_monotonic_microseconds() - start;
    return 1;
}

static symbolic_trie_node **trie_parallel_frontier(symbolic_trie_node *root, int target, int *count) {
    int capacity = target > 8 ? target * 2 : 16;
    symbolic_trie_node **frontier = malloc((size_t) capacity * sizeof(*frontier));
    if (frontier == NULL) return NULL;
    *count = 1; frontier[0] = root;
    for (int current = 0; current < *count && *count < target; ) {
        symbolic_trie_node *node = frontier[current];
        if (node->leaf) { ++current; continue; }
        int children = 0;
        for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL) ++children;
        if (children == 0) { ++current; continue; }
        if (*count - 1 + children > capacity) {
            int new_capacity = capacity * 2;
            while (new_capacity < *count - 1 + children) new_capacity *= 2;
            symbolic_trie_node **expanded = realloc(frontier, (size_t) new_capacity * sizeof(*frontier));
            if (expanded == NULL) { free(frontier); return NULL; }
            frontier = expanded; capacity = new_capacity;
        }
        frontier[current] = frontier[--*count];
        for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL)
            frontier[(*count)++] = node->children[i];
    }
    return frontier;
}

static void trie_parallel_collect_leaves(isax_index *index, const symbolic_trie_node *node,
                                         const ts_type *transform, float bsf,
                                         const symbolic_trie_node *skip_leaf,
                                         trie_leaf_queue *queues, int queue_count,
                                         int *next_queue, trie_query_stats *stats, int *failed) {
    if (node == skip_leaf || __atomic_load_n(failed, __ATOMIC_RELAXED)) return;
    ++stats->checked_nodes;
    float lower = trie_lower_bound(index->trie, index, transform, node->min_word, node->max_word,
                                   index->trie->dimensions, bsf, stats);
    if (lower >= bsf) return;
    if (node->leaf) {
        int queue_number = __sync_fetch_and_add(next_queue, 1) % queue_count;
        if (!trie_leaf_queue_push(&queues[queue_number], (trie_leaf_work) {node, lower}, stats))
            __atomic_store_n(failed, 1, __ATOMIC_RELAXED);
        return;
    }
    for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL)
        trie_parallel_collect_leaves(index, node->children[i], transform, bsf, skip_leaf,
                                     queues, queue_count, next_queue, stats, failed);
}

static float trie_parallel_exact_search(isax_index *index, const ts_type *query,
                                        const ts_type *transform, float bsf,
                                        const symbolic_trie_node *skip_leaf,
                                        trie_query_stats *result_stats,
                                        file_position_type *result_position,
                                        trie_radial_auto_gate *radial_auto_gate) {
    if (result_stats != NULL) memset(result_stats, 0, sizeof(*result_stats));
#ifndef _OPENMP
    trie_query_scratch scratch = {0};
    scratch.radial_auto_gate = radial_auto_gate;
    float distance = trie_search_node(index, index->trie->root, query, transform, bsf,
                                      result_stats, skip_leaf, &scratch, result_position);
    free(scratch.candidates);
    return distance;
#else
    const int workers = maxquerythread > 0 ? maxquerythread : 1;
    /* Match iSAX: --queue-number controls the number of shared leaf-work
     * queues.  The CLI defaults it to the active worker count. */
    const int queue_count = N_PQUEUE > 0 ? N_PQUEUE : workers;
    int frontier_count = 0;
    symbolic_trie_node **frontier = trie_parallel_frontier(index->trie->root, workers * 4, &frontier_count);
    trie_leaf_queue *queues = calloc((size_t) queue_count, sizeof(*queues));
    trie_query_stats *worker_stats = calloc((size_t) workers, sizeof(*worker_stats));
    trie_query_scratch *worker_scratch = calloc((size_t) workers, sizeof(*worker_scratch));
    if (frontier == NULL || queues == NULL || worker_stats == NULL || worker_scratch == NULL) {
        free(frontier); free(queues); free(worker_stats); free(worker_scratch);
        trie_query_scratch scratch = {0};
        scratch.radial_auto_gate = radial_auto_gate;
        float distance = trie_search_node(index, index->trie->root, query, transform, bsf,
                                          result_stats, skip_leaf, &scratch, result_position);
        free(scratch.candidates); return distance;
    }
    for (int i = 0; i < workers; ++i) {
        worker_scratch[i].radial_auto_gate = radial_auto_gate;
        trie_prepare_record_lb_table(index->trie, index, transform, &worker_scratch[i]);
    }
    for (int i = 0; i < queue_count; ++i) pthread_mutex_init(&queues[i].lock, NULL);
    pthread_rwlock_t bsf_lock = PTHREAD_RWLOCK_INITIALIZER;
    float shared_bsf = bsf;
    file_position_type shared_position = result_position == NULL ? QUERY_RESULT_NO_POSITION : *result_position;
    int next_queue = 0, failed = 0;
#pragma omp parallel num_threads(workers)
    {
        int worker = omp_get_thread_num();
        trie_query_stats *stats = &worker_stats[worker];
        unsigned long long worker_start = trie_monotonic_microseconds();
        float traversal_bsf;
        pthread_rwlock_rdlock(&bsf_lock); traversal_bsf = shared_bsf; pthread_rwlock_unlock(&bsf_lock);
        const unsigned long long frontier_start = trie_monotonic_microseconds();
        const unsigned long long mbr_before_frontier = stats->mbr_bound_microseconds;
        const unsigned long long queue_before_frontier = stats->queue_microseconds;
#pragma omp for schedule(dynamic, 1) nowait
        for (int i = 0; i < frontier_count; ++i)
            trie_parallel_collect_leaves(index, frontier[i], transform, traversal_bsf, skip_leaf,
                                         queues, queue_count, &next_queue, stats, &failed);
        if (profile_query_phases) {
            const unsigned long long elapsed = trie_monotonic_microseconds() - frontier_start;
            const unsigned long long direct = (stats->mbr_bound_microseconds - mbr_before_frontier) +
                                              (stats->queue_microseconds - queue_before_frontier);
            stats->frontier_microseconds += elapsed > direct ? elapsed - direct : 0;
        }
#pragma omp barrier
        for (int offset = 0; offset < queue_count && !__atomic_load_n(&failed, __ATOMIC_RELAXED); ++offset) {
            int queue_number = (worker + offset) % queue_count;
            trie_leaf_work work;
            while (trie_leaf_queue_pop(&queues[queue_number], &work, stats)) {
                float current_bsf;
                pthread_rwlock_rdlock(&bsf_lock); current_bsf = shared_bsf; pthread_rwlock_unlock(&bsf_lock);
                if (work.lower_bound > current_bsf) continue;
                file_position_type candidate_position;
                pthread_rwlock_rdlock(&bsf_lock);
                candidate_position = shared_position;
                pthread_rwlock_unlock(&bsf_lock);
                float distance = trie_scan_leaf_best_first(index, work.node, query, transform,
                                                            current_bsf, stats, &worker_scratch[worker],
                                                            &candidate_position);
                if (distance < current_bsf) {
                    pthread_rwlock_wrlock(&bsf_lock);
                    if (distance < shared_bsf) {
                        shared_bsf = distance;
                        shared_position = candidate_position;
                    }
                    pthread_rwlock_unlock(&bsf_lock);
                }
            }
        }
        unsigned long long worker_elapsed = trie_monotonic_microseconds() - worker_start;
        if (profile_query_phases) {
            const unsigned long long bounded = stats->mbr_bound_microseconds +
                                               stats->record_bound_microseconds + stats->exact_distance_microseconds +
                                               stats->frontier_microseconds + stats->queue_microseconds +
                                               stats->candidate_heap_microseconds;
            stats->synchronization_microseconds = worker_elapsed > bounded ? worker_elapsed - bounded : 0;
            stats->traversal_microseconds = stats->frontier_microseconds + stats->queue_microseconds +
                                            stats->candidate_heap_microseconds + stats->synchronization_microseconds;
        }
    }
    if (result_stats != NULL) for (int i = 0; i < workers; ++i) {
        result_stats->checked_nodes += worker_stats[i].checked_nodes;
        result_stats->lower_bounds += worker_stats[i].lower_bounds;
        result_stats->exact_distances += worker_stats[i].exact_distances;
        result_stats->mbr_bound_microseconds += worker_stats[i].mbr_bound_microseconds;
        result_stats->record_bound_microseconds += worker_stats[i].record_bound_microseconds;
        result_stats->exact_distance_microseconds += worker_stats[i].exact_distance_microseconds;
        result_stats->traversal_microseconds += worker_stats[i].traversal_microseconds;
        result_stats->frontier_microseconds += worker_stats[i].frontier_microseconds;
        result_stats->queue_microseconds += worker_stats[i].queue_microseconds;
        result_stats->candidate_heap_microseconds += worker_stats[i].candidate_heap_microseconds;
        result_stats->synchronization_microseconds += worker_stats[i].synchronization_microseconds;
        result_stats->cluster_bounds += worker_stats[i].cluster_bounds;
        result_stats->cluster_pruned += worker_stats[i].cluster_pruned;
        result_stats->cluster_records_pruned += worker_stats[i].cluster_records_pruned;
        result_stats->cluster_symbolic_pruned += worker_stats[i].cluster_symbolic_pruned;
        result_stats->cluster_symbolic_records_pruned += worker_stats[i].cluster_symbolic_records_pruned;
        result_stats->cluster_raw_ball_bounds += worker_stats[i].cluster_raw_ball_bounds;
        result_stats->cluster_raw_ball_pruned += worker_stats[i].cluster_raw_ball_pruned;
        result_stats->cluster_raw_ball_records_pruned += worker_stats[i].cluster_raw_ball_records_pruned;
        result_stats->radial_candidates += worker_stats[i].radial_candidates;
        result_stats->radial_pruned += worker_stats[i].radial_pruned;
    }
    for (int i = 0; i < workers; ++i) free(worker_scratch[i].candidates);
    for (int i = 0; i < queue_count; ++i) { pthread_mutex_destroy(&queues[i].lock); free(queues[i].items); }
    pthread_rwlock_destroy(&bsf_lock);
    free(frontier); free(queues); free(worker_stats); free(worker_scratch);
    if (!__atomic_load_n(&failed, __ATOMIC_RELAXED)) {
        if (result_position != NULL) *result_position = shared_position;
        return shared_bsf;
    }
    /* Queue allocation failed mid-search.  Re-run serially rather than ever
     * returning a partial result. */
    trie_query_scratch scratch = {0};
    scratch.radial_auto_gate = radial_auto_gate;
    float distance = trie_search_node(index, index->trie->root, query, transform, shared_bsf,
                                      result_stats, skip_leaf, &scratch, result_position);
    free(scratch.candidates); return distance;
#endif
}

enum response symbolic_trie_query_file(isax_index *index, const char *path, int query_count,
                                       int filetype_int, int apply_znorm, float minimum_distance) {
    if (index == NULL || path == NULL || query_count < 0) return FAILURE;
    FILE *file = fopen(path, "rb");
    if (file == NULL) return FAILURE;
    if (fseek(file, (long) index->settings->query_header_bytes, SEEK_SET) != 0) {
        fclose(file);
        return FAILURE;
    }
    int ts_length = index->settings->timeseries_size;
    ts_type *query = malloc(sizeof(*query) * (size_t) ts_length);
    file_type *query_int = filetype_int ? malloc(sizeof(*query_int) * (size_t) ts_length) : NULL;
    ts_type *transform = malloc(sizeof(*transform) * (size_t) index->settings->n_segments);
    sax_type *word = malloc((size_t) index->settings->n_segments);
    trie_query_scratch scratch = {0};
    fftw_workspace fftw = {0}; fftw_workspace_init(&fftw, ts_length);
    if (query == NULL || (filetype_int && query_int == NULL) || transform == NULL || word == NULL) {
        fclose(file); free(query); free(query_int); free(transform); free(word); fftw_workspace_destroy(&fftw); return FAILURE;
    }
    unsigned long long cumulative_microseconds = 0;
    for (int i = 0; i < query_count; ++i) {
        int read_ok;
        if (filetype_int) {
            read_ok = fread(query_int, sizeof(*query_int), (size_t) ts_length, file) == (size_t) ts_length;
            for (int j = 0; read_ok && j < ts_length; ++j)
                query[j] = file_value_to_ts(query_int[j], filetype_int);
        } else {
            read_ok = fread(query, sizeof(*query), (size_t) ts_length, file) == (size_t) ts_length;
        }
        if (!read_ok || (apply_znorm && (znorm(query, ts_length), 0))) {
            fclose(file); free(query); free(query_int); free(transform); free(word); free(scratch.candidates);
            fftw_workspace_destroy(&fftw); return FAILURE;
        }
        unsigned long long query_start = trie_monotonic_microseconds();
        if (trie_word_from_ts(index, query, word, transform, &fftw) != SUCCESS) {
            fclose(file); free(query); free(query_int); free(transform); free(word); free(scratch.candidates);
            fftw_workspace_destroy(&fftw); return FAILURE;
        }
        trie_radial_auto_gate radial_auto_gate;
        const int radial_auto_initialized =
            trie_radial_auto_gate_init(index, &radial_auto_gate);
        scratch.radial_auto_gate = radial_auto_initialized ? &radial_auto_gate : NULL;
        trie_prepare_record_lb_table(index->trie, index, transform, &scratch);
        /* The transform drives lower bounds; the word selects the seed path. */
        trie_query_stats stats = {0};
        unsigned long long search_start = trie_monotonic_microseconds();
        const symbolic_trie_node *seed_leaf = trie_seed_leaf(index, word, transform, &stats);
        float bsf = seed_leaf == NULL ? FLT_MAX :
                trie_search_node(index, seed_leaf, query, transform, FLT_MAX, &stats, NULL, &scratch, NULL);
        stats.approximate_distance = bsf;
        float distance;
        if (maxquerythread > 1) {
            if (profile_query_phases) {
                const unsigned long long seed_elapsed = trie_monotonic_microseconds() - search_start;
                const unsigned long long bounded = stats.mbr_bound_microseconds +
                                                   stats.record_bound_microseconds + stats.exact_distance_microseconds +
                                                   stats.candidate_heap_microseconds;
                stats.frontier_microseconds += seed_elapsed > bounded ? seed_elapsed - bounded : 0;
                stats.traversal_microseconds = stats.frontier_microseconds + stats.queue_microseconds +
                                               stats.candidate_heap_microseconds + stats.synchronization_microseconds;
            }
            if (minimum_distance < bsf) bsf = minimum_distance;
            trie_query_stats parallel_stats = {0};
            distance = trie_parallel_exact_search(index, query, transform, bsf, seed_leaf,
                                                  &parallel_stats, NULL,
                                                  scratch.radial_auto_gate);
            stats.checked_nodes += parallel_stats.checked_nodes;
            stats.lower_bounds += parallel_stats.lower_bounds;
            stats.exact_distances += parallel_stats.exact_distances;
            stats.mbr_bound_microseconds += parallel_stats.mbr_bound_microseconds;
            stats.record_bound_microseconds += parallel_stats.record_bound_microseconds;
            stats.exact_distance_microseconds += parallel_stats.exact_distance_microseconds;
            stats.traversal_microseconds += parallel_stats.traversal_microseconds;
            stats.frontier_microseconds += parallel_stats.frontier_microseconds;
            stats.queue_microseconds += parallel_stats.queue_microseconds;
            stats.candidate_heap_microseconds += parallel_stats.candidate_heap_microseconds;
            stats.synchronization_microseconds += parallel_stats.synchronization_microseconds;
            stats.cluster_bounds += parallel_stats.cluster_bounds;
            stats.cluster_pruned += parallel_stats.cluster_pruned;
            stats.cluster_records_pruned += parallel_stats.cluster_records_pruned;
            stats.cluster_symbolic_pruned += parallel_stats.cluster_symbolic_pruned;
            stats.cluster_symbolic_records_pruned += parallel_stats.cluster_symbolic_records_pruned;
            stats.cluster_raw_ball_bounds += parallel_stats.cluster_raw_ball_bounds;
            stats.cluster_raw_ball_pruned += parallel_stats.cluster_raw_ball_pruned;
            stats.cluster_raw_ball_records_pruned += parallel_stats.cluster_raw_ball_records_pruned;
            stats.radial_candidates += parallel_stats.radial_candidates;
            stats.radial_pruned += parallel_stats.radial_pruned;
        } else {
            distance = trie_search_node(index, index->trie->root, query, transform,
                                        minimum_distance < bsf ? minimum_distance : bsf, &stats, seed_leaf, &scratch, NULL);
        }
        if (profile_query_phases && maxquerythread <= 1) {
            const unsigned long long search_elapsed = trie_monotonic_microseconds() - search_start;
            const unsigned long long bounded = stats.mbr_bound_microseconds +
                                               stats.record_bound_microseconds + stats.exact_distance_microseconds +
                                               stats.candidate_heap_microseconds;
            stats.frontier_microseconds += search_elapsed > bounded ? search_elapsed - bounded : 0;
            stats.traversal_microseconds = stats.frontier_microseconds + stats.queue_microseconds +
                                           stats.candidate_heap_microseconds + stats.synchronization_microseconds;
        }
        trie_radial_auto_gate_destroy(&radial_auto_gate, radial_auto_initialized);
        scratch.radial_auto_gate = NULL;
        stats.total_microseconds = trie_monotonic_microseconds() - query_start;
        cumulative_microseconds += stats.total_microseconds;
        trie_save_query_stats(index->trie, &stats, distance, cumulative_microseconds);
        if (SHOULD_REPORT_QUERY(i, query_count)) {
            PRINT_STATS_HEADER();
            trie_print_query_stats(i, index->trie, &stats, distance, cumulative_microseconds);
        }
    }
    fclose(file); free(query); free(query_int); free(transform); free(word); free(scratch.candidates);
    fftw_workspace_destroy(&fftw); return SUCCESS;
}

enum response symbolic_trie_query_file_batch(isax_index *index, const char *path, int query_count,
                                             int filetype_int, int apply_znorm, float minimum_distance) {
    if (index == NULL || path == NULL || query_count < 0) return FAILURE;
    FILE *file = fopen(path, "rb");
    int length = index->settings->timeseries_size, dimensions = index->settings->n_segments;
    ts_type *queries = malloc((size_t) query_count * length * sizeof(*queries));
    float *distances = malloc((size_t) query_count * sizeof(*distances));
    trie_query_stats *stats = calloc((size_t) query_count, sizeof(*stats));
    const size_t query_values = (size_t) query_count * length;
    if (file == NULL || queries == NULL || distances == NULL || stats == NULL) {
        if (file) fclose(file); free(queries); free(distances); free(stats); return FAILURE;
    }
    if (fseek(file, (long) index->settings->query_header_bytes, SEEK_SET) != 0) {
        fclose(file); free(queries); free(distances); free(stats); return FAILURE;
    }
    int read_ok = 1;
    if (filetype_int) {
        file_type *input = malloc(query_values * sizeof(*input));
        if (input == NULL || fread(input, sizeof(*input), query_values, file) != query_values) {
            read_ok = 0;
        } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(maxquerythread)
#endif
            for (size_t i = 0; i < query_values; ++i)
                queries[i] = file_value_to_ts(input[i], filetype_int);
        }
        free(input);
    } else if (fread(queries, sizeof(*queries), query_values, file) != query_values) {
        read_ok = 0;
    }
    fclose(file);
    if (!read_ok) { free(queries); free(distances); free(stats); return FAILURE; }
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
            unsigned long long query_start = trie_monotonic_microseconds();
            if (apply_znorm) znorm(query, length);
            if (transform == NULL || word == NULL || trie_word_from_ts(index, query, word, transform, &fftw) != SUCCESS) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                failed = 1;
            } else {
                trie_radial_auto_gate radial_auto_gate;
                const int radial_auto_initialized =
                    trie_radial_auto_gate_init(index, &radial_auto_gate);
                scratch.radial_auto_gate = radial_auto_initialized ? &radial_auto_gate : NULL;
                trie_prepare_record_lb_table(index->trie, index, transform, &scratch);
                /* One worker owns one query: avoid nested per-query tasks. */
                unsigned long long search_start = trie_monotonic_microseconds();
                const symbolic_trie_node *seed_leaf = trie_seed_leaf(index, word, transform, &stats[i]);
                float bsf = seed_leaf == NULL ? FLT_MAX :
                        trie_search_node(index, seed_leaf, query, transform, FLT_MAX, &stats[i], NULL, &scratch, NULL);
                stats[i].approximate_distance = bsf;
                if (minimum_distance < bsf) bsf = minimum_distance;
                distances[i] = trie_search_node(index, index->trie->root, query, transform, bsf,
                                                &stats[i], seed_leaf, &scratch, NULL);
                if (profile_query_phases) {
                    const unsigned long long search_elapsed = trie_monotonic_microseconds() - search_start;
                    const unsigned long long bounded = stats[i].mbr_bound_microseconds +
                                                       stats[i].record_bound_microseconds + stats[i].exact_distance_microseconds +
                                                       stats[i].candidate_heap_microseconds;
                    stats[i].frontier_microseconds += search_elapsed > bounded ? search_elapsed - bounded : 0;
                    stats[i].traversal_microseconds = stats[i].frontier_microseconds + stats[i].queue_microseconds +
                                                      stats[i].candidate_heap_microseconds + stats[i].synchronization_microseconds;
                }
                trie_radial_auto_gate_destroy(&radial_auto_gate, radial_auto_initialized);
                scratch.radial_auto_gate = NULL;
                stats[i].total_microseconds = trie_monotonic_microseconds() - query_start;
            }
        }
        free(transform); free(word); free(scratch.candidates);
        pthread_mutex_lock(&trie_fftw_plan_lock);
        fftw_workspace_destroy(&fftw);
        pthread_mutex_unlock(&trie_fftw_plan_lock);
    }
    if (!failed) {
        unsigned long long cumulative_microseconds = 0;
        for (int i = 0; i < query_count; ++i) {
            cumulative_microseconds += stats[i].total_microseconds;
            trie_save_query_stats(index->trie, &stats[i], distances[i], cumulative_microseconds);
            if (SHOULD_REPORT_QUERY(i, query_count)) {
                PRINT_STATS_HEADER();
                trie_print_query_stats(i, index->trie, &stats[i], distances[i], cumulative_microseconds);
            }
        }
        fflush(stderr);
    }
    free(queries); free(distances); free(stats); return failed ? FAILURE : SUCCESS;
}

query_result symbolic_trie_exact_search(isax_index *index, const ts_type *query,
                                        const ts_type *transform, float bsf) {
    query_result result = { FLT_MAX, NULL, 0, QUERY_RESULT_NO_POSITION };
    if (index == NULL || index->trie == NULL) return result;
    trie_radial_auto_gate radial_auto_gate;
    const int radial_auto_initialized = trie_radial_auto_gate_init(index, &radial_auto_gate);
    if (maxquerythread <= 1) {
        trie_query_scratch scratch = {0};
        scratch.radial_auto_gate = radial_auto_initialized ? &radial_auto_gate : NULL;
        trie_prepare_record_lb_table(index->trie, index, transform, &scratch);
        result.distance = trie_search_node(index, index->trie->root, query, transform, bsf,
                                           NULL, NULL, &scratch, &result.record_position);
        free(scratch.candidates);
    }
    else result.distance = trie_parallel_exact_search(index, query, transform, bsf, NULL, NULL,
                                                       &result.record_position,
                                                       radial_auto_initialized ? &radial_auto_gate : NULL);
    trie_radial_auto_gate_destroy(&radial_auto_gate, radial_auto_initialized);
    return result;
}

void symbolic_trie_destroy(isax_index *index) {
    if (index != NULL && index->trie != NULL) {
        trie_node_destroy(index->trie->root);
        free(index->trie->word_arena);
        free(index->trie); index->trie = NULL;
    }
}
