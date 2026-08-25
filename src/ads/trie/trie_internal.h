#ifndef MESSI_SYMBOLIC_TRIE_INTERNAL_H
#define MESSI_SYMBOLIC_TRIE_INTERNAL_H

/* Private implementation contract shared by trie.c and trie_topk.c.
 * Public users must include ads/trie/trie.h instead. */

#include "ads/trie/trie.h"
#include "ads/calc_utils.h"

#include <stdint.h>
#include <pthread.h>

#define TRIE_MAX_FANOUT 256

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
    trie_dimension_mask used_dimensions;
    unsigned char leaf;
    unsigned char mbb_valid;
    unsigned char split_exhausted;
    trie_leaf_cluster *clusters;
    sax_type *cluster_min_words;
    sax_type *cluster_max_words;
    int cluster_count;
} symbolic_trie_node;

struct symbolic_trie_index {
    symbolic_trie_node *root;
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
};

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
    float approximate_distance;
} trie_query_stats;

typedef struct {
    float lower_bound;
    int record_index;
} trie_leaf_candidate;

typedef struct {
    trie_leaf_candidate *candidates;
    int capacity;
    float record_lb_table[MESSI_RECORD_LB_MAX_DIMENSIONS][256];
    int record_lb_table_ready;
} trie_query_scratch;

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

float trie_lower_bound(const struct symbolic_trie_index *trie, isax_index *index,
                       const ts_type *transform, sax_type *sax_min, sax_type *sax_max,
                       int dimensions, float bsf, trie_query_stats *stats);
float trie_record_lower_bound(const struct symbolic_trie_index *trie, isax_index *index,
                              const ts_type *transform, sax_type *word, float bsf,
                              float mbr_suffix, const trie_query_scratch *scratch);
float trie_record_mbr_suffix(const struct symbolic_trie_index *trie, isax_index *index,
                             const ts_type *transform, const sax_type *min_word,
                             const sax_type *max_word, const trie_query_scratch *scratch);
void trie_prepare_record_lb_table(const struct symbolic_trie_index *trie,
                                  const isax_index *index, const ts_type *transform,
                                  trie_query_scratch *scratch);
int trie_scratch_reserve(trie_query_scratch *scratch, int capacity);
float trie_exact_distance(const ts_type *query, const ts_type *record, int length,
                          float bsf, trie_query_stats *stats);
void trie_heap_build(trie_leaf_candidate *heap, int size);
trie_leaf_candidate trie_heap_pop(trie_leaf_candidate *heap, int *size);
enum response trie_word_from_transform(isax_index *index, const ts_type *transform,
                                       sax_type *word);
const symbolic_trie_node *trie_seed_leaf(isax_index *index, const sax_type *query_word,
                                         const ts_type *transform, trie_query_stats *stats);
int trie_leaf_queue_pop(trie_leaf_queue *queue, trie_leaf_work *item, trie_query_stats *stats);
symbolic_trie_node **trie_parallel_frontier(symbolic_trie_node *root, int target, int *count);
void trie_parallel_collect_leaves(isax_index *index, const symbolic_trie_node *node,
                                  const ts_type *transform, float bsf,
                                  const symbolic_trie_node *skip_leaf,
                                  trie_leaf_queue *queues, int queue_count, int *next_queue,
                                  trie_query_stats *stats, int *failed);

symbolic_trie_topk_result trie_topk_exact_search_with_stats(
    isax_index *index, const ts_type *query, const ts_type *transform,
    size_t k, float minimum_distance, trie_query_stats *stats);

#endif
