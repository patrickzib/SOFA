#include "config.h"
#include "globals.h"
#include "trie_internal.h"
#include "ads/inmemory_index_engine.h"
#include "ads/parallel_query_engine.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* The 1-NN engine intentionally remains in trie.c with its scalar BSF hot
 * path.  This file owns the bounded max-heap required by exact top-k. */
typedef struct { float distance; file_position_type position; } trie_topk_entry;
typedef struct { trie_topk_entry *entries; size_t size; size_t capacity; } trie_topk_heap;

static int trie_topk_entry_worse(const trie_topk_entry *left, const trie_topk_entry *right) {
    return left->distance > right->distance ||
           (left->distance == right->distance && left->position > right->position);
}
static void trie_topk_sift_down(trie_topk_heap *heap, size_t root) {
    while (1) { size_t child = root * 2 + 1;
        if (child >= heap->size) return;
        if (child + 1 < heap->size && trie_topk_entry_worse(&heap->entries[child + 1], &heap->entries[child])) ++child;
        if (!trie_topk_entry_worse(&heap->entries[child], &heap->entries[root])) return;
        trie_topk_entry temporary = heap->entries[root]; heap->entries[root] = heap->entries[child];
        heap->entries[child] = temporary; root = child;
    }
}
static void trie_topk_sift_up(trie_topk_heap *heap, size_t child) {
    while (child != 0) { size_t parent = (child - 1) / 2;
        if (!trie_topk_entry_worse(&heap->entries[child], &heap->entries[parent])) return;
        trie_topk_entry temporary = heap->entries[parent]; heap->entries[parent] = heap->entries[child];
        heap->entries[child] = temporary; child = parent;
    }
}
static float trie_topk_threshold(const trie_topk_heap *heap) {
    return heap->size < heap->capacity ? FLT_MAX : heap->entries[0].distance;
}
static int trie_topk_insert(trie_topk_heap *heap, float distance, file_position_type position) {
    trie_topk_entry entry = { distance, position };
    for (size_t i = 0; i < heap->size; ++i) {
        if (heap->entries[i].position != position) continue;
        if (!trie_topk_entry_worse(&heap->entries[i], &entry)) return 0;
        heap->entries[i] = entry;
        for (size_t root = heap->size / 2; root-- > 0; ) trie_topk_sift_down(heap, root);
        return 1;
    }
    if (heap->size < heap->capacity) { heap->entries[heap->size] = entry; trie_topk_sift_up(heap, heap->size++); return 1; }
    if (!trie_topk_entry_worse(&heap->entries[0], &entry)) return 0;
    heap->entries[0] = entry; trie_topk_sift_down(heap, 0); return 1;
}
static int trie_topk_compare_ascending(const void *left, const void *right) {
    const trie_topk_entry *a = left, *b = right;
    if (a->distance < b->distance) return -1;
    if (a->distance > b->distance) return 1;
    return a->position < b->position ? -1 : a->position > b->position;
}
static float trie_topk_shared_threshold(const trie_topk_heap *local, trie_topk_heap *shared,
                                        pthread_mutex_t *shared_lock) {
    float threshold = trie_topk_threshold(local);
    if (shared != NULL && shared_lock != NULL) {
        pthread_mutex_lock(shared_lock); float global_threshold = trie_topk_threshold(shared); pthread_mutex_unlock(shared_lock);
        if (global_threshold < threshold) threshold = global_threshold;
    }
    return threshold;
}
static void trie_topk_publish(const trie_topk_heap *local, trie_topk_heap *shared, pthread_mutex_t *shared_lock) {
    if (shared == NULL || shared_lock == NULL || local->size == 0) return;
    pthread_mutex_lock(shared_lock);
    for (size_t i = 0; i < local->size; ++i) trie_topk_insert(shared, local->entries[i].distance, local->entries[i].position);
    pthread_mutex_unlock(shared_lock);
}

static int trie_topk_scan_leaf_range(isax_index *index, const symbolic_trie_node *node,
                                     int offset, int count, float mbr_suffix,
                                     const ts_type *query, const ts_type *transform,
                                     trie_topk_heap *results, trie_topk_heap *shared,
                                     pthread_mutex_t *shared_lock, trie_query_stats *stats,
                                     trie_query_scratch *scratch) {
    const int end = offset + count; int changed = 0;
    const float published_threshold = trie_topk_shared_threshold(results, shared, shared_lock);
    if (!trie_scratch_reserve(scratch, count)) {
        for (int i = offset; i < end; ++i) {
            if (stats != NULL) ++stats->lower_bounds;
            float threshold = trie_topk_threshold(results); if (published_threshold < threshold) threshold = published_threshold;
            float lower = trie_record_lower_bound(index->trie, index, transform, node->words[i], threshold, mbr_suffix, scratch);
            if (lower > threshold) continue;
            if (stats != NULL) ++stats->exact_distances;
            float d = trie_exact_distance(query, rawfile + node->positions[i], index->settings->timeseries_size,
                                          nextafterf(threshold, FLT_MAX), stats);
            if (d <= threshold && trie_topk_insert(results, d, node->positions[i])) changed = 1;
        }
        return changed;
    }
    int candidate_count = 0;
    for (int i = offset; i < end; ++i) {
        if (stats != NULL) ++stats->lower_bounds;
        float threshold = trie_topk_threshold(results); if (published_threshold < threshold) threshold = published_threshold;
        float lower = trie_record_lower_bound(index->trie, index, transform, node->words[i], threshold, mbr_suffix, scratch);
        if (lower <= threshold) scratch->candidates[candidate_count++] = (trie_leaf_candidate) { lower, i };
    }
    trie_heap_build(scratch->candidates, candidate_count);
    while (candidate_count > 0) {
        float threshold = trie_topk_threshold(results); if (published_threshold < threshold) threshold = published_threshold;
        if (scratch->candidates[0].lower_bound > threshold) break;
        trie_leaf_candidate candidate = trie_heap_pop(scratch->candidates, &candidate_count);
        if (stats != NULL) ++stats->exact_distances;
        file_position_type position = node->positions[candidate.record_index];
        float d = trie_exact_distance(query, rawfile + position, index->settings->timeseries_size,
                                      nextafterf(threshold, FLT_MAX), stats);
        if (d <= threshold && trie_topk_insert(results, d, position)) changed = 1;
    }
    return changed;
}

static int trie_topk_scan_leaf(isax_index *index, const symbolic_trie_node *node,
                               const ts_type *query, const ts_type *transform,
                               trie_topk_heap *results, trie_topk_heap *shared,
                               pthread_mutex_t *shared_lock, trie_query_stats *stats,
                               trie_query_scratch *scratch) {
    int changed = 0;
    if (node->cluster_count == 0) {
        float suffix = trie_record_mbr_suffix(index->trie, index, transform, node->min_word, node->max_word, scratch);
        return trie_topk_scan_leaf_range(index, node, 0, node->size, suffix, query, transform,
                                         results, shared, shared_lock, stats, scratch);
    }
    typedef struct { int cluster; float bound; } trie_topk_cluster_bound;
    trie_topk_cluster_bound ordered[64]; int count = 0;
    for (int cluster = 0; cluster < node->cluster_count; ++cluster) {
        float threshold = trie_topk_shared_threshold(results, shared, shared_lock);
        const trie_leaf_cluster *group = &node->clusters[cluster];
        if (stats != NULL) ++stats->cluster_bounds;
        float bound = trie_lower_bound(index->trie, index, transform,
            node->cluster_min_words + (size_t) cluster * index->trie->dimensions,
            node->cluster_max_words + (size_t) cluster * index->trie->dimensions,
            index->trie->dimensions, threshold, stats);
        if (bound > threshold) { if (stats != NULL) { ++stats->cluster_pruned; stats->cluster_records_pruned += group->size; } }
        else ordered[count++] = (trie_topk_cluster_bound) { cluster, bound };
    }
    for (int i = 1; i < count; ++i) { trie_topk_cluster_bound value = ordered[i]; int j = i - 1;
        while (j >= 0 && ordered[j].bound > value.bound) { ordered[j + 1] = ordered[j]; --j; } ordered[j + 1] = value; }
    for (int i = 0; i < count; ++i) {
        int cluster = ordered[i].cluster; const trie_leaf_cluster *group = &node->clusters[cluster];
        const sax_type *minimum = node->cluster_min_words + (size_t) cluster * index->trie->dimensions;
        const sax_type *maximum = node->cluster_max_words + (size_t) cluster * index->trie->dimensions;
        float suffix = trie_record_mbr_suffix(index->trie, index, transform, minimum, maximum, scratch);
        changed |= trie_topk_scan_leaf_range(index, node, group->offset, group->size, suffix, query, transform,
                                             results, shared, shared_lock, stats, scratch);
    }
    return changed;
}

static void trie_topk_search_node(isax_index *index, const symbolic_trie_node *node,
                                  const ts_type *query, const ts_type *transform,
                                  trie_topk_heap *results, trie_query_stats *stats,
                                  trie_query_scratch *scratch) {
    if (stats != NULL) ++stats->checked_nodes;
    float threshold = trie_topk_threshold(results);
    float lower = trie_lower_bound(index->trie, index, transform, node->min_word, node->max_word,
                                   index->trie->dimensions, threshold, stats);
    if (lower > threshold) return;
    if (node->leaf) { trie_topk_scan_leaf(index, node, query, transform, results, NULL, NULL, stats, scratch); return; }
    typedef struct { const symbolic_trie_node *node; float bound; } trie_topk_child_bound;
    trie_topk_child_bound ordered[TRIE_MAX_FANOUT]; int count = 0;
    for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL) {
        threshold = trie_topk_threshold(results);
        float bound = trie_lower_bound(index->trie, index, transform, node->children[i]->min_word,
                                       node->children[i]->max_word, index->trie->dimensions, threshold, stats);
        if (bound <= threshold) ordered[count++] = (trie_topk_child_bound) { node->children[i], bound };
    }
    for (int i = 1; i < count; ++i) { trie_topk_child_bound value = ordered[i]; int j = i - 1;
        while (j >= 0 && ordered[j].bound > value.bound) { ordered[j + 1] = ordered[j]; --j; } ordered[j + 1] = value; }
    for (int i = 0; i < count; ++i)
        if (ordered[i].bound <= trie_topk_threshold(results)) trie_topk_search_node(index, ordered[i].node, query, transform, results, stats, scratch);
}
static const symbolic_trie_node *trie_topk_first_leaf(const symbolic_trie_node *node) {
    while (node != NULL && !node->leaf) { const symbolic_trie_node *next = NULL;
        for (int i = 0; i < node->split_fanout; ++i) if (node->children[i] != NULL) { next = node->children[i]; break; }
        node = next;
    }
    return node;
}

static int trie_topk_parallel_search(isax_index *index, const ts_type *query, const ts_type *transform,
                                     trie_topk_heap *results, trie_query_stats *result_stats) {
    if (result_stats != NULL) memset(result_stats, 0, sizeof(*result_stats));
#ifndef _OPENMP
    trie_query_scratch scratch = {0}; trie_prepare_record_lb_table(index->trie, index, transform, &scratch);
    trie_topk_search_node(index, index->trie->root, query, transform, results, result_stats, &scratch);
    free(scratch.candidates); return results->size == results->capacity;
#else
    const int workers = maxquerythread > 0 ? maxquerythread : 1;
    if (workers <= 1) { trie_query_scratch scratch = {0}; trie_prepare_record_lb_table(index->trie, index, transform, &scratch);
        trie_topk_search_node(index, index->trie->root, query, transform, results, result_stats, &scratch);
        free(scratch.candidates); return results->size == results->capacity; }
    const int queue_count = N_PQUEUE > 0 ? N_PQUEUE : workers; int frontier_count = 0;
    symbolic_trie_node **frontier = trie_parallel_frontier(index->trie->root, workers * 4, &frontier_count);
    trie_leaf_queue *queues = calloc((size_t) queue_count, sizeof(*queues));
    trie_query_scratch *scratch = calloc((size_t) workers, sizeof(*scratch));
    trie_topk_heap *local = calloc((size_t) workers, sizeof(*local));
    trie_query_stats *worker_stats = calloc((size_t) workers, sizeof(*worker_stats));
    if (frontier == NULL || queues == NULL || scratch == NULL || local == NULL || worker_stats == NULL) goto serial_fallback;
    sax_type *seed_word = malloc((size_t) index->settings->n_segments);
    const symbolic_trie_node *seed = seed_word != NULL && trie_word_from_transform(index, transform, seed_word) == SUCCESS
        ? trie_seed_leaf(index, seed_word, transform, result_stats) : trie_topk_first_leaf(index->trie->root);
    free(seed_word);
    trie_query_scratch seed_scratch = {0}; trie_prepare_record_lb_table(index->trie, index, transform, &seed_scratch);
    if (seed != NULL) trie_topk_scan_leaf(index, seed, query, transform, results, NULL, NULL, result_stats, &seed_scratch);
    free(seed_scratch.candidates);
    pthread_mutex_t result_lock = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < queue_count; ++i) pthread_mutex_init(&queues[i].lock, NULL);
    for (int i = 0; i < workers; ++i) { local[i].entries = calloc(results->capacity, sizeof(*local[i].entries)); local[i].capacity = results->capacity;
        if (local[i].entries == NULL) { for (int j = 0; j < i; ++j) free(local[j].entries); goto queue_fallback; }
        trie_prepare_record_lb_table(index->trie, index, transform, &scratch[i]); }
    { int next_queue = 0, failed = 0; const float initial_threshold = trie_topk_threshold(results);
#pragma omp parallel num_threads(workers)
      { int worker = omp_get_thread_num();
#pragma omp for schedule(dynamic, 1) nowait
        for (int i = 0; i < frontier_count; ++i) trie_parallel_collect_leaves(index, frontier[i], transform, initial_threshold, seed, queues, queue_count, &next_queue, &worker_stats[worker], &failed);
#pragma omp barrier
        for (int offset = 0; offset < queue_count && !__atomic_load_n(&failed, __ATOMIC_RELAXED); ++offset) { trie_leaf_work work; int queue = (worker + offset) % queue_count;
          while (trie_leaf_queue_pop(&queues[queue], &work, &worker_stats[worker])) { float threshold = trie_topk_shared_threshold(&local[worker], results, &result_lock);
            if (work.lower_bound > threshold) continue;
            if (trie_topk_scan_leaf(index, work.node, query, transform, &local[worker], results, &result_lock, &worker_stats[worker], &scratch[worker])) trie_topk_publish(&local[worker], results, &result_lock);
          }
        }
      }
      for (int i = 0; i < workers; ++i) { trie_topk_publish(&local[i], results, &result_lock);
        if (result_stats != NULL) { result_stats->checked_nodes += worker_stats[i].checked_nodes; result_stats->lower_bounds += worker_stats[i].lower_bounds; result_stats->exact_distances += worker_stats[i].exact_distances; result_stats->cluster_bounds += worker_stats[i].cluster_bounds; result_stats->cluster_pruned += worker_stats[i].cluster_pruned; result_stats->cluster_records_pruned += worker_stats[i].cluster_records_pruned; }
        free(local[i].entries); free(scratch[i].candidates); }
      for (int i = 0; i < queue_count; ++i) { pthread_mutex_destroy(&queues[i].lock); free(queues[i].items); }
      pthread_mutex_destroy(&result_lock); free(frontier); free(queues); free(scratch); free(local); free(worker_stats);
      if (!__atomic_load_n(&failed, __ATOMIC_RELAXED)) return results->size == results->capacity;
      frontier = NULL; queues = NULL; scratch = NULL; local = NULL; worker_stats = NULL;
    }
    results->size = 0;
serial_fallback:
    free(frontier); free(queues); free(scratch); free(local); free(worker_stats);
    { trie_query_scratch serial_scratch = {0}; trie_prepare_record_lb_table(index->trie, index, transform, &serial_scratch);
      trie_topk_search_node(index, index->trie->root, query, transform, results, result_stats, &serial_scratch);
      free(serial_scratch.candidates); return results->size == results->capacity; }
queue_fallback:
    for (int i = 0; i < queue_count; ++i) pthread_mutex_destroy(&queues[i].lock);
    pthread_mutex_destroy(&result_lock); free(frontier); free(queues); free(scratch); free(local); free(worker_stats);
    { trie_query_scratch serial_scratch = {0}; trie_prepare_record_lb_table(index->trie, index, transform, &serial_scratch);
      trie_topk_search_node(index, index->trie->root, query, transform, results, result_stats, &serial_scratch);
      free(serial_scratch.candidates); return results->size == results->capacity; }
#endif
}

symbolic_trie_topk_result trie_topk_exact_search_with_stats(isax_index *index, const ts_type *query,
    const ts_type *transform, size_t k, float minimum_distance, trie_query_stats *stats) {
    symbolic_trie_topk_result output = {0, NULL, NULL};
    if (index == NULL || index->trie == NULL || query == NULL || transform == NULL || k == 0 || k > (size_t) index->total_records) return output;
    trie_topk_heap heap = { calloc(k, sizeof(*heap.entries)), 0, k }; if (heap.entries == NULL) return output;
    (void) minimum_distance;
    if (!trie_topk_parallel_search(index, query, transform, &heap, stats)) { free(heap.entries); return output; }
    qsort(heap.entries, heap.size, sizeof(*heap.entries), trie_topk_compare_ascending);
    output.distances = malloc(heap.size * sizeof(*output.distances)); output.positions = malloc(heap.size * sizeof(*output.positions));
    if (output.distances == NULL || output.positions == NULL) { symbolic_trie_topk_result_destroy(&output); free(heap.entries); return output; }
    for (size_t i = 0; i < heap.size; ++i) { output.distances[i] = heap.entries[i].distance; output.positions[i] = heap.entries[i].position / index->settings->timeseries_size; }
    output.count = heap.size; free(heap.entries); return output;
}
symbolic_trie_topk_result symbolic_trie_exact_topk_search(isax_index *index, const ts_type *query,
    const ts_type *transform, size_t k, float minimum_distance) {
    return trie_topk_exact_search_with_stats(index, query, transform, k, minimum_distance, NULL);
}
void symbolic_trie_topk_result_destroy(symbolic_trie_topk_result *result) {
    if (result == NULL) return;
    free(result->distances); free(result->positions); result->distances = NULL; result->positions = NULL; result->count = 0;
}
