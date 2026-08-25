#ifdef VALUES
#include <values.h>
#endif

#include <float.h>
#include "config.h"
#include "../../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <sys/wait.h>
#include <time.h>
#include <sys/time.h>

#include "omp.h"
#include "ads/isax_query_engine.h"
#include "ads/inmemory_index_engine.h"
#include "ads/inmemory_query_engine.h"
#include "ads/parallel_inmemory_query_engine.h"
#include "ads/parallel_index_engine.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/pqueue.h"
#include "ads/sax/sax.h"
#include "ads/sfa/sfa.h"
#include "ads/sfa/dft.h"
#include "ads/calc_utils.h"
#include "ads/isax_node_split.h"
#include "ads/inmemory_topk_engine.h"
#include "ads/pthread_barrier.h"

#define NTHREADS 4

int N_PQUEUE = 1;

static void *exact_search_old_worker_inmemory(void *rfdata);
static void *mindistance_worker_inmemory(void *essdata);
static void *readworker_inmemory(void *read_pointer);
static void *topk_readworker_inmemory(void *read_pointer);
static void *exact_search_worker_inmemory_hybridpqueue(void *rfdata);
static void insert_tree_node_m_hybridpqueue_time(float *paa, isax_node *node,
                                                  isax_index *index, float bsf, pqueue_t **pq,
                                                  pthread_mutex_t *lock_queue, int *tnumber,
                                                  unsigned long int *time_lb);


static float calculate_node_distance_inmemory_m(isax_index *index, isax_node *node,
                                                ts_type *query, ts_type *paa, float bsf) {
    float distmin;
    COUNT_CHECKED_NODE()

    // If node has buffered data
    if (node->buffer != NULL) {
        int i;

        if (node->buffer->full_buffer_size > 0) {
#pragma omp parallel for num_threads(maxquerythread) reduction(min : bsf)
            for (i = 0; i < node->buffer->full_buffer_size; i++) {
                distmin = messi_minidist_raw(index, paa, node->buffer->full_sax_buffer[i],
                                             index->settings->max_sax_cardinalities, bsf);

                if (distmin < bsf) {
                    float dist = ts_ed(query, node->buffer->full_ts_buffer[i],
                                   index->settings->timeseries_size, bsf);
                    if (dist < bsf) {
                        bsf = dist;
                    }
                }
            }
        }

        if (node->buffer->tmp_full_buffer_size > 0) {
#pragma omp parallel for num_threads(maxquerythread) reduction(min : bsf)
            for (i = 0; i < node->buffer->tmp_full_buffer_size; i++) {
                distmin = messi_minidist_raw(index, paa, node->buffer->tmp_full_sax_buffer[i],
                                             index->settings->max_sax_cardinalities, bsf);

                if (distmin < bsf) {
                    float dist = ts_ed(query, node->buffer->tmp_full_ts_buffer[i],
                                   index->settings->timeseries_size, bsf);
                    if (dist < bsf) {
                        bsf = dist;
                    }
                }
            }
        }

        if (node->buffer->partial_buffer_size > 0) {
#pragma omp parallel for num_threads(maxquerythread) reduction(min : bsf)
            for (i = 0; i < node->buffer->partial_buffer_size; i++) {
                distmin = messi_minidist_raw(index, paa, node->buffer->partial_sax_buffer[i],
                                             index->settings->max_sax_cardinalities, bsf);

                if (distmin < bsf) {
                    float dist = ts_ed(query, &(rawfile[*node->buffer->partial_position_buffer[i]]),
                                       index->settings->timeseries_size, bsf);
                    if (dist < bsf) {
                        bsf = dist;
                    }
                }
            }
        }
    }

    //////////////////////////////////////

    return bsf;
}


static query_result refine_answer_inmemory_m(ts_type *ts, ts_type *paa, isax_index *index,
                                             query_result approximate_bsf_result,
                                             float minimum_distance, int limit) {
    query_result bsf_result = approximate_bsf_result;

    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;

    int j = 0;
    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size,
                               cmp_pri, get_pri, set_pri, get_pos, set_pos);

    // Insert all root nodes in heap.
    isax_node *current_root_node = index->first_node;

    while (current_root_node != NULL) {
        query_result *mindist_result = malloc(sizeof(query_result));

        mindist_result->distance = messi_minidist(index, paa, current_root_node->isax_values,
                                                  current_root_node->isax_cardinalities, minimum_distance);
        mindist_result->node = current_root_node;
        pqueue_insert(pq, mindist_result);
        current_root_node = current_root_node->next;
    }
    query_result *n;
    int checks = 0;
    while ((n = pqueue_pop(pq))) {
        // The best node has a worse mindist, so search is finished!
        if (n->distance >= bsf_result.distance || n->distance > minimum_distance) {
            pqueue_insert(pq, n);
            break;
        } else {
            // If it is a leaf, check its real distance.
            if (n->node->is_leaf) {
                // *** ADAPTIVE SPLITTING ***
                if (!n->node->has_full_data_file &&
                    (n->node->leaf_size > index->settings->min_leaf_size)) {
                    free(n);
                    continue;
                }
                // *** EXTRA BOUNDING ***
                if (tight_bound) {
                    j++;
                    float mindistance = calculate_minimum_distance_inmemory(index, n->node, ts, paa);

                    if (mindistance >= bsf_result.distance) {
                        free(n);
                        continue;
                    }
                }
                // *** REAL DISTANCE ***
                checks++;
                float distance = calculate_node_distance_inmemory_m(index, n->node, ts, paa, bsf_result.distance);
                if (distance < bsf_result.distance) {
                    bsf_result.distance = distance;
                    bsf_result.node = n->node;
                }
                if (checks > limit) {
                    pqueue_insert(pq, n);
                    break;
                }
            } else {
                // If it is an intermediate node calculate mindist for children
                // and push them in the queue
                if (n->node->left_child->isax_cardinalities != NULL) {
                    if (n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file &&
                        aggressive_check) {
                        float distance = calculate_node_distance_inmemory(index, n->node->left_child, ts, paa,
                                                                          bsf_result.distance);
                        if (distance < bsf_result.distance) {
                            bsf_result.distance = distance;
                            bsf_result.node = n->node->left_child;
                        }
                    } else {
                        query_result *mindist_result = malloc(sizeof(query_result));

                        mindist_result->distance = messi_minidist(index, paa,
                                                                  n->node->left_child->isax_values,
                                                                  n->node->left_child->isax_cardinalities,
                                                                  minimum_distance);
                        mindist_result->node = n->node->left_child;
                        pqueue_insert(pq, mindist_result);
                    }
                }
                if (n->node->right_child->isax_cardinalities != NULL) {
                    if (n->node->right_child->is_leaf && !n->node->left_child->has_partial_data_file &&
                        aggressive_check) {
                        float distance = calculate_node_distance_inmemory(index, n->node->right_child, ts, paa,
                                                                          bsf_result.distance);
                        if (distance < bsf_result.distance) {
                            bsf_result.distance = distance;
                            bsf_result.node = n->node->right_child;
                        }
                    } else {
                        query_result *mindist_result = malloc(sizeof(query_result));

                        mindist_result->distance = messi_minidist(index, paa,
                                                                  n->node->right_child->isax_values,
                                                                  n->node->right_child->isax_cardinalities,
                                                                  minimum_distance);
                        mindist_result->node = n->node->right_child;
                        pqueue_insert(pq, mindist_result);
                    }
                }
            }

            // Free the node currently popped.
            free(n);
        }
    }
    // Free the nodes that where not popped.
    while ((n = pqueue_pop(pq))) {
        free(n);
    }
    // Free the priority queue.
    pqueue_free(pq);
    return bsf_result;
}


query_result exact_search_parads_inmemory(ts_type *ts, ts_type *paa, isax_index *index,
                                          float minimum_distance, int min_checked_leaves) {
    query_result approximate_result = approximate_search_inmemory_pRecBuf(ts, paa, index, 1);
    query_result bsf_result = approximate_result;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    int i;
    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }
    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer_inmemory_m(ts, paa, index, approximate_result, minimum_distance,
                                                      min_checked_leaves);
    }

    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size,
                               cmp_pri, get_pri, set_pri, get_pos, set_pos);


    query_result *do_not_remove = &approximate_result;

    SET_APPROXIMATE(approximate_result.distance);

    RESET_BYTES_ACCESSED

    if (approximate_result.node != NULL) {
        // Insert approximate result in heap.
        pqueue_insert(pq, &approximate_result);
    }

    // Insert all root nodes in heap.
    isax_node *current_root_node = index->first_node;


    pthread_t threadid[maxquerythread];

    refind_answer_fonction_data rfdata;
    pthread_mutex_t lock_queue = PTHREAD_MUTEX_INITIALIZER, lock_current_root_node = PTHREAD_MUTEX_INITIALIZER;
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;
    pthread_barrier_t lock_barrier;
    pthread_barrier_init(&lock_barrier, NULL, 1);
    rfdata.paa = paa;
    rfdata.ts = ts;

    rfdata.lock_queue = &lock_queue;
    rfdata.lock_current_root_node = &lock_current_root_node;
    rfdata.lock_bsf = &lock_bsf;

    rfdata.index = index;
    rfdata.minimum_distance = minimum_distance;

    rfdata.pq = pq;
    // Insert all root nodes in heap.
    rfdata.current_root_node = index->first_node;
    rfdata.bsf_result = &bsf_result;

    query_result *n;

    rfdata.lock_barrier = &lock_barrier;
    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, exact_search_old_worker_inmemory, (void *) &(rfdata));
    }
    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }


    // Free the nodes that where not popped.
    while ((n = pqueue_pop(pq))) {
        if (n != do_not_remove)
            free(n);
    }
    // Free the priority queue.
    pthread_barrier_destroy(&lock_barrier);
    pqueue_free(pq);

    bsf_result = *(rfdata.bsf_result);
    return *(rfdata.bsf_result);
}

static void *exact_search_old_worker_inmemory(void *rfdata) {
    isax_node *current_root_node;
    query_result *n;
    isax_index *index = ((refind_answer_fonction_data *) rfdata)->index;
    ts_type *paa = ((refind_answer_fonction_data *) rfdata)->paa;
    ts_type *ts = ((refind_answer_fonction_data *) rfdata)->ts;
    pqueue_t *pq = ((refind_answer_fonction_data *) rfdata)->pq;
    query_result *do_not_remove = ((refind_answer_fonction_data *) rfdata)->bsf_result;
    float minimum_distance = ((refind_answer_fonction_data *) rfdata)->minimum_distance;
    float bsfdisntance;
    int limit = ((refind_answer_fonction_data *) rfdata)->limit;
    int checks = 0;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    query_result *bsf_result = (((refind_answer_fonction_data *) rfdata)->bsf_result);
    while (1) {
        pthread_mutex_lock(((refind_answer_fonction_data *) rfdata)->lock_current_root_node);
        current_root_node = ((refind_answer_fonction_data *) rfdata)->current_root_node;
        if (current_root_node != NULL) {
            ((refind_answer_fonction_data *) rfdata)->current_root_node = ((refind_answer_fonction_data *) rfdata)->current_root_node->next;
            pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_current_root_node);
        } else {
            pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_current_root_node);
            break;
        }


        query_result *mindist_result = malloc(sizeof(query_result));
        mindist_result->distance = minidist_paa_to_isax(paa, current_root_node->isax_values, current_root_node->isax_cardinalities, index->settings, 0);

        mindist_result->node = current_root_node;
        pthread_mutex_lock(((refind_answer_fonction_data *) rfdata)->lock_queue);

        pqueue_insert(pq, mindist_result);

        pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_queue);
    }
    pthread_barrier_wait(((refind_answer_fonction_data *) rfdata)->lock_barrier);

    while (1) {


        pthread_mutex_lock(((refind_answer_fonction_data *) rfdata)->lock_queue);
        n = pqueue_pop(pq);
        pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_queue);
        if (n == NULL)
            break;
        pthread_rwlock_rdlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
        bsfdisntance = bsf_result->distance;
        pthread_rwlock_unlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
        // The best node has a worse mindist, so search is finished!
        if (n->distance >= bsfdisntance || n->distance > minimum_distance) {
            pthread_mutex_lock(((refind_answer_fonction_data *) rfdata)->lock_queue);
            pqueue_insert(pq, n);
            pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_queue);
            break;
        } else {
            // If it is a leaf, check its real distance.
            if (n->node->is_leaf) {
                // *** ADAPTIVE SPLITTING ***
                if (!n->node->has_full_data_file &&
                    (n->node->leaf_size > index->settings->min_leaf_size)) {
                    split_node(index, n->node, 1, 1);
                    pthread_mutex_lock(((refind_answer_fonction_data *) rfdata)->lock_queue);
                    pqueue_insert(pq, n);
                    pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_queue);
                    continue;
                }
                // *** EXTRA BOUNDING ***
                if (tight_bound) {
                    float mindistance = calculate_minimum_distance_inmemory(index, n->node, ts, paa);
                    if (mindistance >= bsfdisntance) {
                        if (n != do_not_remove)
                            free(n);
                        continue;
                    }
                }
                // *** REAL DISTANCE ***
                checks++;
                float distance = calculate_node_distance_inmemory(index, n->node, ts, paa, bsfdisntance);
                if (distance < bsfdisntance) {
                    pthread_rwlock_wrlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
                    if (distance < bsf_result->distance) {
                        bsf_result->distance = distance;
                        bsf_result->node = n->node;
                    }

                    pthread_rwlock_unlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
                }
            } else {
                // If it is an intermediate node calculate mindist for children
                // and push them in the queue
                if (n->node->left_child->isax_cardinalities != NULL) {
                    if (n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file &&
                        aggressive_check) {
                        float distance = calculate_node_distance_inmemory(index, n->node, ts, paa, bsfdisntance);
                        if (distance < bsfdisntance) {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
                            if (distance < bsf_result->distance) {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->left_child;
                            }
                            pthread_rwlock_unlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
                        }
                    } else {
                        query_result *mindist_result = malloc(sizeof(query_result));
                        mindist_result->distance = minidist_paa_to_isax(paa, n->node->left_child->isax_values, n->node->left_child->isax_cardinalities, index->settings, 0);
                        mindist_result->node = n->node->left_child;
                        pthread_mutex_lock(((refind_answer_fonction_data *) rfdata)->lock_queue);
                        pqueue_insert(pq, mindist_result);
                        pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_queue);
                    }
                }
                if (n->node->right_child->isax_cardinalities != NULL) {
                    if (n->node->right_child->is_leaf && !n->node->left_child->has_partial_data_file &&
                        aggressive_check) {
                        float distance = calculate_node_distance_inmemory(index, n->node, ts, paa, bsfdisntance);
                        if (distance < bsfdisntance) {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
                            if (distance < bsf_result->distance) {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->right_child;
                            }
                            pthread_rwlock_unlock(((refind_answer_fonction_data *) rfdata)->lock_bsf);
                        }
                    } else {
                        query_result *mindist_result = malloc(sizeof(query_result));
                        mindist_result->distance = minidist_paa_to_isax(paa, n->node->right_child->isax_values, n->node->right_child->isax_cardinalities, index->settings, 0);
                        mindist_result->node = n->node->right_child;
                        pthread_mutex_lock(((refind_answer_fonction_data *) rfdata)->lock_queue);
                        pqueue_insert(pq, mindist_result);
                        pthread_mutex_unlock(((refind_answer_fonction_data *) rfdata)->lock_queue);
                    }
                }
            }

            // Free the node currently popped.
            if (n != do_not_remove)
                free(n);
        }
    }
    return NULL;
}

query_result exact_search_serial_ParIS_inmemory(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance,
                                                int min_checked_leaves) {

    RESET_BYTES_ACCESSED


    pthread_t threadid[maxquerythread];
    COUNT_INPUT_TIME_START
    query_result approximate_result = approximate_search_inmemory_pRecBuf(ts, paa, index, 1);
    int sum_of_lab = 0;
    float bsf_distance;
    unsigned long j;

    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }

    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer_inmemory(ts, paa, index, approximate_result, minimum_distance,
                                                    min_checked_leaves);
    }

    COUNT_INPUT_TIME_END
    int i;

#ifdef AUTO_TUNE
    float *mindists = malloc(sizeof(float) * index->sax_cache_size);
#endif

    SET_APPROXIMATE(approximate_result.distance);
    ParIS_LDCW_data *essdata = malloc(sizeof(ParIS_LDCW_data) * (maxquerythread));
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;
    COUNT_CAL_TIME_START
    for (i = 0; i < (maxquerythread - 1); i++) {
        essdata[i].index = index;
        essdata[i].lock_bsf = &lock_bsf;
        essdata[i].start_number = i * (index->sax_cache_size / maxquerythread);
        essdata[i].stop_number = (i + 1) * (index->sax_cache_size / maxquerythread);
        essdata[i].paa = paa;
        essdata[i].ts = ts;
        essdata[i].bsfdistance = approximate_result.distance;
        essdata[i].sum_of_lab = 0;
    }
    essdata[maxquerythread - 1].index = index;
    essdata[maxquerythread - 1].lock_bsf = &lock_bsf;
    essdata[maxquerythread - 1].start_number = (maxquerythread - 1) * (index->sax_cache_size / maxquerythread);
    essdata[maxquerythread - 1].stop_number = index->sax_cache_size;
    essdata[maxquerythread - 1].paa = paa;
    essdata[maxquerythread - 1].ts = ts;
    essdata[maxquerythread - 1].bsfdistance = approximate_result.distance;
    essdata[maxquerythread - 1].sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, mindistance_worker_inmemory, (void *) &(essdata[i]));
    }
    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
        sum_of_lab += essdata[i].sum_of_lab;
    }
    COUNT_CAL_TIME_END
    unsigned long *label_number = malloc(sizeof(unsigned long) * (sum_of_lab));
    float *minidisvector = malloc(sizeof(float) * (sum_of_lab));

    sum_of_lab = 0;
    COUNT_OUTPUT_TIME_START
    for (i = 0; i < maxquerythread; i++) {
        memcpy(&(label_number[sum_of_lab]), essdata[i].label_number, sizeof(unsigned long) * essdata[i].sum_of_lab);
        memcpy(&(minidisvector[sum_of_lab]), essdata[i].minidisvector, sizeof(float) * essdata[i].sum_of_lab);
        free(essdata[i].label_number);
        free(essdata[i].minidisvector);
        sum_of_lab += essdata[i].sum_of_lab;

    }
    pthread_t readthread[maxquerythread * maxreadthread];
    ParIS_read_worker_data readpointer;
    unsigned long readcounter = 0;
    float bsfdistance = (approximate_result.distance);


    readpointer.ts = ts;
    readpointer.index = index;
    readpointer.counter = &readcounter;
    readpointer.load_point = label_number;
    readpointer.lock_bsf = &lock_bsf;
    readpointer.bsf2 = &bsfdistance;
    readpointer.minidisvector = minidisvector;
    readpointer.sum_of_lab = sum_of_lab;
    readpointer.rawfile = rawfile;

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_create(&(readthread[i]), NULL, readworker_inmemory, (void *) &(readpointer));
    }

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_join(readthread[i], NULL);
    }
    COUNT_OUTPUT_TIME_END
    approximate_result.distance = bsfdistance;
    free(essdata);
    free(minidisvector);
    free(label_number);
    return approximate_result;
}

pqueue_bsf exact_topk_serial_ParIS_inmemory(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance,
                                            int min_checked_leaves, int k) {

    RESET_BYTES_ACCESSED

    pqueue_bsf *pq_bsf = pqueue_bsf_init(k);
    pthread_t threadid[maxquerythread];
    COUNT_INPUT_TIME_START
    approximate_topk_inmemory(ts, paa, index, pq_bsf);

    int sum_of_lab = 0;

    unsigned long j;

    // Early termination...


    if (pq_bsf->knn[k - 1] == FLT_MAX || min_checked_leaves > 1) {
        refine_topk_answer_inmemory(ts, paa, index, pq_bsf, minimum_distance, min_checked_leaves);
    }

    COUNT_INPUT_TIME_END
    int i;

#ifdef AUTO_TUNE
    float *mindists = malloc(sizeof(float) * index->sax_cache_size);
#endif

    SET_APPROXIMATE(pq_bsf->knn[k - 1]);
    ParIS_LDCW_data *essdata = malloc(sizeof(ParIS_LDCW_data) * (maxquerythread));
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;
    COUNT_CAL_TIME_START
    for (i = 0; i < (maxquerythread - 1); i++) {
        essdata[i].index = index;
        essdata[i].lock_bsf = &lock_bsf;
        essdata[i].start_number = i * (index->sax_cache_size / maxquerythread);
        essdata[i].stop_number = (i + 1) * (index->sax_cache_size / maxquerythread);
        essdata[i].paa = paa;
        essdata[i].ts = ts;
        essdata[i].bsfdistance = pq_bsf->knn[k - 1];
        essdata[i].sum_of_lab = 0;
    }
    essdata[maxquerythread - 1].index = index;
    essdata[maxquerythread - 1].lock_bsf = &lock_bsf;
    essdata[maxquerythread - 1].start_number = (maxquerythread - 1) * (index->sax_cache_size / maxquerythread);
    essdata[maxquerythread - 1].stop_number = index->sax_cache_size;
    essdata[maxquerythread - 1].paa = paa;
    essdata[maxquerythread - 1].ts = ts;
    essdata[maxquerythread - 1].bsfdistance = pq_bsf->knn[k - 1];
    essdata[maxquerythread - 1].sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, mindistance_worker_inmemory, (void *) &(essdata[i]));
    }
    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
        sum_of_lab += essdata[i].sum_of_lab;
    }
    COUNT_CAL_TIME_END
    unsigned long *label_number = malloc(sizeof(unsigned long) * (sum_of_lab));
    float *minidisvector = malloc(sizeof(float) * (sum_of_lab));

    sum_of_lab = 0;
    COUNT_OUTPUT_TIME_START
    for (i = 0; i < maxquerythread; i++) {
        memcpy(&(label_number[sum_of_lab]), essdata[i].label_number, sizeof(unsigned long) * essdata[i].sum_of_lab);
        memcpy(&(minidisvector[sum_of_lab]), essdata[i].minidisvector, sizeof(float) * essdata[i].sum_of_lab);
        free(essdata[i].label_number);
        free(essdata[i].minidisvector);
        sum_of_lab += essdata[i].sum_of_lab;

    }
    pthread_t readthread[maxquerythread * maxreadthread];
    ParIS_read_worker_data readpointer;
    unsigned long readcounter = 0;

    readpointer.ts = ts;
    readpointer.index = index;
    readpointer.counter = &readcounter;
    readpointer.load_point = label_number;
    readpointer.lock_bsf = &lock_bsf;
    readpointer.minidisvector = minidisvector;
    readpointer.sum_of_lab = sum_of_lab;
    readpointer.rawfile = rawfile;
    readpointer.pq_bsf = pq_bsf;


    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_create(&(readthread[i]), NULL, topk_readworker_inmemory, (void *) &(readpointer));
    }

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_join(readthread[i], NULL);
    }
    COUNT_OUTPUT_TIME_END
    free(essdata);
    free(minidisvector);
    free(label_number);
    return *pq_bsf;
}


static void *mindistance_worker_inmemory(void *essdata) {

    unsigned long i;
    float bsfdistance, mindist;
    isax_index *index = ((ParIS_LDCW_data *) essdata)->index;
    unsigned long start_number = ((ParIS_LDCW_data *) essdata)->start_number;
    unsigned long stop_number = ((ParIS_LDCW_data *) essdata)->stop_number;
    ts_type *paa = ((ParIS_LDCW_data *) essdata)->paa;
    ts_type *ts = ((ParIS_LDCW_data *) essdata)->ts;
    ((ParIS_LDCW_data *) essdata)->label_number = malloc(sizeof(unsigned long) * 10000);
    ((ParIS_LDCW_data *) essdata)->minidisvector = malloc(sizeof(float) * 10000);

    int max_number = 10000;

    for (i = start_number; i < stop_number; i++) {

        sax_type *sax = &index->sax_cache[i * index->settings->n_segments];

        mindist = minidist_paa_to_isax_rawa_SIMD(paa, sax, index->settings->max_sax_cardinalities, index->settings);
        if (mindist <= ((ParIS_LDCW_data *) essdata)->bsfdistance) {
            if (((ParIS_LDCW_data *) essdata)->sum_of_lab >= max_number) {
                max_number = (max_number + 10000);
                ((ParIS_LDCW_data *) essdata)->label_number = (unsigned long *) realloc(
                        ((ParIS_LDCW_data *) essdata)->label_number, sizeof(unsigned long) * max_number);
                ((ParIS_LDCW_data *) essdata)->minidisvector = (float *) realloc(
                        ((ParIS_LDCW_data *) essdata)->minidisvector, sizeof(float) * max_number);

            }
            ((ParIS_LDCW_data *) essdata)->label_number[((ParIS_LDCW_data *) essdata)->sum_of_lab] = i;
            ((ParIS_LDCW_data *) essdata)->minidisvector[((ParIS_LDCW_data *) essdata)->sum_of_lab] = mindist;
            ((ParIS_LDCW_data *) essdata)->sum_of_lab++;
        }
    }
    return NULL;
}

static void *readworker_inmemory(void *read_pointer) {
    isax_index *index = ((ParIS_read_worker_data *) read_pointer)->index;
    float bsf, dist;
    ts_type *ts_buffer;
    ts_type *ts = ((ParIS_read_worker_data *) read_pointer)->ts;
    unsigned long t = 0, p;
    unsigned long sum_of_lab = ((ParIS_read_worker_data *) read_pointer)->sum_of_lab;
    float *minidisvector = ((ParIS_read_worker_data *) read_pointer)->minidisvector;
    bsf = *(((ParIS_read_worker_data *) read_pointer)->bsf2);

    while (1) {
        pthread_rwlock_rdlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
        bsf = *(((ParIS_read_worker_data *) read_pointer)->bsf2);
        pthread_rwlock_unlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
        t = __sync_fetch_and_add(((ParIS_read_worker_data *) read_pointer)->counter, 1);
        if (t >= sum_of_lab) {
            break;
        }
        p = ((ParIS_read_worker_data *) read_pointer)->load_point[t];
        if (minidisvector[t] < bsf) {
            ts_buffer = &((ParIS_read_worker_data *) read_pointer)->rawfile[p * index->settings->ts_byte_size /
                                                                            sizeof(ts_type)];
            dist = ts_ed(ts, ts_buffer, index->settings->timeseries_size, bsf);

            if (dist < bsf) {
                pthread_rwlock_wrlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
                if (dist < *(((ParIS_read_worker_data *) read_pointer)->bsf2)) {
                    *(((ParIS_read_worker_data *) read_pointer)->bsf2) = dist;
                }
                pthread_rwlock_unlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
            }
        }
    }
    return NULL;
}

static void *topk_readworker_inmemory(void *read_pointer) {
    isax_index *index = ((ParIS_read_worker_data *) read_pointer)->index;
    float bsf, dist;
    pqueue_bsf *pq_bsf = ((ParIS_read_worker_data *) read_pointer)->pq_bsf;
    ts_type *ts_buffer;
    ts_type *ts = ((ParIS_read_worker_data *) read_pointer)->ts;
    unsigned long t = 0, p;
    unsigned long sum_of_lab = ((ParIS_read_worker_data *) read_pointer)->sum_of_lab;
    float *minidisvector = ((ParIS_read_worker_data *) read_pointer)->minidisvector;
    bsf = FLT_MAX;

    while (1) {
        pthread_rwlock_rdlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
        bsf = pq_bsf->knn[pq_bsf->k - 1];
        pthread_rwlock_unlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
        t = __sync_fetch_and_add(((ParIS_read_worker_data *) read_pointer)->counter, 1);
        if (t >= sum_of_lab) {
            break;
        }
        p = ((ParIS_read_worker_data *) read_pointer)->load_point[t];
        if (minidisvector[t] < bsf) {
            ts_buffer = &((ParIS_read_worker_data *) read_pointer)->rawfile[p * index->settings->ts_byte_size /
                                                                            sizeof(ts_type)];
            dist = ts_ed(ts, ts_buffer, index->settings->timeseries_size, bsf);

            if (dist <= bsf) {
                pthread_rwlock_wrlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
                pqueue_bsf_insert(pq_bsf, dist, p, NULL);
                pthread_rwlock_unlock(((ParIS_read_worker_data *) read_pointer)->lock_bsf);
            }
        }
    }
    return NULL;
}

query_result exact_search_MESSI(ts_type *ts, ts_type *paa, isax_index *index, node_list *nodelist,
                                float minimum_distance, int min_checked_leaves, int kn) {
    query_result approximate_result = approximate_search_inmemory_pRecBuf(ts, paa, index, kn);

    query_result bsf_result = approximate_result;
    int node_counter = 0;
    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }
    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer_inmemory_m(
                ts, paa, index, approximate_result, minimum_distance,min_checked_leaves);
    }
    pqueue_t **allpq = malloc(sizeof(pqueue_t *) * N_PQUEUE);

    pthread_mutex_t ququelock[N_PQUEUE];
    int queuelabel[N_PQUEUE];

    SET_APPROXIMATE(approximate_result.distance);

    pthread_t threadid[maxquerythread];
    MESSI_workerdata workerdata[maxquerythread];
    pthread_mutex_t lock_queue = PTHREAD_MUTEX_INITIALIZER, lock_current_root_node = PTHREAD_MUTEX_INITIALIZER;
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;
    pthread_barrier_t lock_barrier;
    pthread_barrier_init(&lock_barrier, NULL, maxquerythread);
    if (profile_query_phases) {
        TOTAL_TREE_TRAVERSAL_TIME = 0;
        TOTAL_LB_DIST_CALC_TIME = 0;
        TOTAL_MBR_DIST_CALC_TIME = 0;
        TOTAL_RECORD_LB_DIST_CALC_TIME = 0;
        TOTAL_REAL_DIST_CALC_TIME = 0;
    }

    for (int i = 0; i < N_PQUEUE; i++) {
        allpq[i] = pqueue_init(index->settings->root_nodes_size / N_PQUEUE,
                               cmp_pri, get_pri, set_pri, get_pos, set_pos);
        pthread_mutex_init(&ququelock[i], NULL);
        queuelabel[i] = 1;
    }

    if (nodelist == NULL) {
        fprintf(stderr, "exact_search_MESSI called with NULL nodelist.\n");
        return approximate_result;
    }

    for (int i = 0; i < maxquerythread; i++) {
        workerdata[i].paa = paa;
        workerdata[i].ts = ts;
        workerdata[i].lock_queue = &lock_queue;
        workerdata[i].lock_current_root_node = &lock_current_root_node;
        workerdata[i].lock_bsf = &lock_bsf;
        workerdata[i].nodelist = nodelist->nlist;
        workerdata[i].amountnode = nodelist->node_amount;
        workerdata[i].index = index;
        workerdata[i].minimum_distance = minimum_distance;
        workerdata[i].node_counter = &node_counter;
        workerdata[i].pq = allpq[i];
        workerdata[i].bsf_result = &bsf_result;
        workerdata[i].lock_barrier = &lock_barrier;
        workerdata[i].alllock = ququelock;
        workerdata[i].allqueuelabel = queuelabel;
        workerdata[i].allpq = allpq;
        workerdata[i].startqueuenumber = i % N_PQUEUE;
    }

    for (int i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, exact_search_worker_inmemory_hybridpqueue, (void *) &(workerdata[i]));
    }
    for (int i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }
    if (profile_query_phases)
        total_tree_pass_time = (double) TOTAL_TREE_TRAVERSAL_TIME;

    pthread_barrier_destroy(&lock_barrier);

    for (int i = 0; i < N_PQUEUE; i++) {
        pqueue_free(allpq[i]);
    }
    free(allpq);
    bsf_result = bsf_result;

    return bsf_result;

}


static float calculate_node_distance2_inmemory_profiled(isax_index *index, isax_node *node,
                                                         ts_type *query, ts_type *paa, float bsf,
                                                         unsigned long int *record_bound_time,
                                                         unsigned long int *exact_distance_time) {
    if (node->buffer == NULL) return bsf;
    COUNT_CHECKED_NODE()
    int exact_count = 0;
    for (int i = 0; i < node->buffer->partial_buffer_size; ++i) {
        struct timeval start, end;
        gettimeofday(&start, NULL);
        float lower = messi_minidist_raw(index, paa, node->buffer->partial_sax_buffer[i],
                                         index->settings->max_sax_cardinalities, bsf);
        gettimeofday(&end, NULL);
        *record_bound_time += (unsigned long int)
            ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec - start.tv_usec);
        if (lower < bsf) {
            gettimeofday(&start, NULL);
            float distance = ts_ed(query, rawfile + *node->buffer->partial_position_buffer[i],
                                   index->settings->timeseries_size, bsf);
            gettimeofday(&end, NULL);
            *exact_distance_time += (unsigned long int)
                ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec - start.tv_usec);
            ++exact_count;
            if (distance < bsf) bsf = distance;
        }
    }
    __sync_fetch_and_add(&LBDcalculationnumber, node->buffer->partial_buffer_size);
    __sync_fetch_and_add(&RDcalculationnumber, exact_count);
    return bsf;
}
static void *exact_search_worker_inmemory_hybridpqueue(void *rfdata) {
    isax_node *current_root_node;
    query_result *n;
    isax_index *index = ((MESSI_workerdata *) rfdata)->index;
    ts_type *paa = ((MESSI_workerdata *) rfdata)->paa;
    ts_type *ts = ((MESSI_workerdata *) rfdata)->ts;
    float minimum_distance = ((MESSI_workerdata *) rfdata)->minimum_distance;
    bool finished = true;
    int current_root_node_number;
    query_result *bsf_result = (((MESSI_workerdata *) rfdata)->bsf_result);
    float bsfdisntance = bsf_result->distance;
    int tnumber = rand() % N_PQUEUE;
    int startqueuenumber = ((MESSI_workerdata *) rfdata)->startqueuenumber;

    struct timeval current_time;
    struct timeval pq_insert_time_start;
    struct timeval pq_remove_time_start;

    unsigned long int total_pq_insert_time=0;
    unsigned long int total_pq_remove_time=0;
    unsigned long int total_mbr_dist_calc_time=0;
    unsigned long int total_record_lb_dist_calc_time=0;
    unsigned long int total_real_dist_calc_time=0;
    unsigned long int total_tree_traversal_time=0;
    struct timeval worker_phase_start;
    if (profile_query_phases) gettimeofday(&worker_phase_start, NULL);


    gettimeofday(&pq_insert_time_start, NULL);

    while (1) {
        current_root_node_number = __sync_fetch_and_add(((MESSI_workerdata *) rfdata)->node_counter, 1);
        if (current_root_node_number >= ((MESSI_workerdata *) rfdata)->amountnode)
            break;
        current_root_node = ((MESSI_workerdata *) rfdata)->nodelist[current_root_node_number];

        if (profile_query_phases) {
            insert_tree_node_m_hybridpqueue_time(
                paa, current_root_node, index, bsfdisntance,
                ((MESSI_workerdata *) rfdata)->allpq, ((MESSI_workerdata *) rfdata)->alllock,
                &tnumber, &total_mbr_dist_calc_time);
        } else {
            insert_tree_node_m_hybridpqueue(
                paa, current_root_node, index, bsfdisntance,
                ((MESSI_workerdata *) rfdata)->allpq, ((MESSI_workerdata *) rfdata)->alllock,
                &tnumber);
        }
    }

    pthread_barrier_wait(((MESSI_workerdata *) rfdata)->lock_barrier);

    gettimeofday(&current_time, NULL);
    total_pq_insert_time += ((current_time.tv_sec*1000000 + (current_time.tv_usec)) -
        (pq_insert_time_start.tv_sec*1000000 + (pq_insert_time_start.tv_usec)));

    while (1) {
        gettimeofday(&pq_remove_time_start, NULL);
        pthread_mutex_lock(&(((MESSI_workerdata *) rfdata)->alllock[startqueuenumber]));
        n = pqueue_pop(((MESSI_workerdata *) rfdata)->allpq[startqueuenumber]);
        pthread_mutex_unlock(&(((MESSI_workerdata *) rfdata)->alllock[startqueuenumber]));

        gettimeofday(&current_time, NULL);
        total_pq_remove_time += ((current_time.tv_sec*1000000 + (current_time.tv_usec)) -
            (pq_remove_time_start.tv_sec*1000000 + (pq_remove_time_start.tv_usec)));

        if (n == NULL) {
            break;
        }

        bsfdisntance = bsf_result->distance;
        if (n->distance > bsfdisntance || n->distance > minimum_distance) {
            break;
        }

        // If it is a leaf, check its real distance.
        if (n->node->is_leaf) {
            float distance = profile_query_phases
                             ? calculate_node_distance2_inmemory_profiled(index, n->node, ts, paa,
                                                                          bsfdisntance, &total_record_lb_dist_calc_time,
                                                                          &total_real_dist_calc_time)
                             : calculate_node_distance2_inmemory(index, n->node, ts, paa, bsfdisntance);

            if (distance < bsfdisntance) {
                pthread_rwlock_wrlock(((MESSI_workerdata *) rfdata)->lock_bsf);
                if (distance < bsf_result->distance) {
                    bsf_result->distance = distance;
                    bsf_result->node = n->node;
                    bsfdisntance = distance;
                }
                pthread_rwlock_unlock(((MESSI_workerdata *) rfdata)->lock_bsf);
            }
        }

        free(n);
    }

    if ((((MESSI_workerdata *) rfdata)->allqueuelabel[startqueuenumber]) == 1) {
        gettimeofday(&pq_remove_time_start, NULL);

        (((MESSI_workerdata *) rfdata)->allqueuelabel[startqueuenumber]) = 0;
        pthread_mutex_lock(&(((MESSI_workerdata *) rfdata)->alllock[startqueuenumber]));
        while ((n = pqueue_pop(((MESSI_workerdata *) rfdata)->allpq[startqueuenumber]))) {
            free(n);
        }
        pthread_mutex_unlock(&(((MESSI_workerdata *) rfdata)->alllock[startqueuenumber]));

        gettimeofday(&current_time, NULL);
        total_pq_remove_time += ((current_time.tv_sec*1000000 + (current_time.tv_usec)) - (pq_remove_time_start.tv_sec*1000000 + (pq_remove_time_start.tv_usec)));
    }

    while (1) {
        finished = true;
        for (int i = 0; i < N_PQUEUE; i++) {
            if ((((MESSI_workerdata *) rfdata)->allqueuelabel[i]) == 1) {
                finished = false;
                while (1) {
                    gettimeofday(&pq_remove_time_start, NULL);

                    pthread_mutex_lock(&(((MESSI_workerdata *) rfdata)->alllock[i]));
                    n = pqueue_pop(((MESSI_workerdata *) rfdata)->allpq[i]);
                    pthread_mutex_unlock(&(((MESSI_workerdata *) rfdata)->alllock[i]));

                    gettimeofday(&current_time, NULL);
                    total_pq_remove_time += ((current_time.tv_sec*1000000 + (current_time.tv_usec))
                        - (pq_remove_time_start.tv_sec*1000000 + (pq_remove_time_start.tv_usec)));

                    if (n == NULL) {
                        ((MESSI_workerdata *) rfdata)->allqueuelabel[i] = 0;
                        break;
                    }

                    bsfdisntance = bsf_result->distance;
                    if (n->distance > bsfdisntance || n->distance > minimum_distance) {
                        break;
                    }

                    // If it is a leaf, check its real distance.
                    if (n->node->is_leaf) {
                        float distance = profile_query_phases
                                         ? calculate_node_distance2_inmemory_profiled(index, n->node, ts, paa,
                                                                                      bsf_result->distance,
                                                                                      &total_record_lb_dist_calc_time,
                                                                                      &total_real_dist_calc_time)
                                         : calculate_node_distance2_inmemory(index, n->node, ts, paa,
                                                                            bsf_result->distance);

                        if (distance < bsfdisntance) {
                            pthread_rwlock_wrlock(((MESSI_workerdata *) rfdata)->lock_bsf);
                            if (distance < bsf_result->distance) {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node;
                                bsfdisntance = distance;
                            }
                            pthread_rwlock_unlock(((MESSI_workerdata *) rfdata)->lock_bsf);
                        }
                    }

                    free(n);
                }

            }
        }
        if (finished) {
            break;
        }
    }

    __sync_fetch_and_add(&TOTAL_PQ_INSERT_TIME,(int)total_pq_insert_time);
    __sync_fetch_and_add(&TOTAL_PQ_REMOVE_TIME,(int)total_pq_remove_time);
    __sync_fetch_and_add(&TOTAL_MBR_DIST_CALC_TIME,(int)total_mbr_dist_calc_time);
    __sync_fetch_and_add(&TOTAL_RECORD_LB_DIST_CALC_TIME,(int)total_record_lb_dist_calc_time);
    __sync_fetch_and_add(&TOTAL_LB_DIST_CALC_TIME,
                         (int)(total_mbr_dist_calc_time + total_record_lb_dist_calc_time));
    __sync_fetch_and_add(&TOTAL_REAL_DIST_CALC_TIME,(int)total_real_dist_calc_time);
    if (profile_query_phases) {
        struct timeval worker_phase_end;
        gettimeofday(&worker_phase_end, NULL);
        const unsigned long int worker_elapsed = (unsigned long int)
            ((worker_phase_end.tv_sec - worker_phase_start.tv_sec) * 1000000L +
             worker_phase_end.tv_usec - worker_phase_start.tv_usec);
        const unsigned long int direct = total_mbr_dist_calc_time +
                                         total_record_lb_dist_calc_time + total_real_dist_calc_time;
        total_tree_traversal_time = worker_elapsed > direct ? worker_elapsed - direct : 0;
    }
    __sync_fetch_and_add(&TOTAL_TREE_TRAVERSAL_TIME,(int)total_tree_traversal_time);
    return NULL;
}

void insert_tree_node_m_hybridpqueue(float *paa, isax_node *node, isax_index *index, float bsf, pqueue_t **pq,
                                     pthread_mutex_t *lock_queue, int *tnumber) {

    float distance = messi_minidist_range_raw(
        index, paa,
        node->mbb_sax_min, node->mbb_sax_max,
        index->settings->max_sax_cardinalities, bsf);

    if (distance < bsf) {
        if (node->is_leaf) {
            query_result *mindist_result = malloc(sizeof(query_result));
            mindist_result->node = node;
            mindist_result->distance = distance;
            pthread_mutex_lock(&lock_queue[*tnumber]);
            pqueue_insert(pq[*tnumber], mindist_result);
            pthread_mutex_unlock(&lock_queue[*tnumber]);
            *tnumber = (*tnumber + 1) % N_PQUEUE;
            added_tree_node++;
        } else {
            if (node->left_child->isax_cardinalities != NULL) {
                insert_tree_node_m_hybridpqueue(paa, node->left_child, index, bsf, pq, lock_queue, tnumber);
            }
            if (node->right_child->isax_cardinalities != NULL) {
                insert_tree_node_m_hybridpqueue(paa, node->right_child, index, bsf, pq, lock_queue, tnumber);
            }
        }
    }
}

static void insert_tree_node_m_hybridpqueue_time(float *paa, isax_node *node,
                                                  isax_index *index, float bsf, pqueue_t **pq,
                                                  pthread_mutex_t *lock_queue, int *tnumber,
                                                  unsigned long int *time_lb) {
    struct timeval current_time;
    struct timeval lb_dist_time_start;

    gettimeofday(&lb_dist_time_start, NULL);
    float distance = messi_minidist_range_raw(
            index, paa,
            node->mbb_sax_min, node->mbb_sax_max,
            index->settings->max_sax_cardinalities, bsf);

    gettimeofday(&current_time, NULL);
    *time_lb += ((int) (current_time.tv_sec * 1000000 + (current_time.tv_usec)) -
                 (int) (lb_dist_time_start.tv_sec * 1000000 + (lb_dist_time_start.tv_usec)));

    if (distance < bsf) {
        if (node->is_leaf) {
            query_result *mindist_result = malloc(sizeof(query_result));
            mindist_result->node = node;
            mindist_result->distance = distance;
            pthread_mutex_lock(&lock_queue[*tnumber]);
            pqueue_insert(pq[*tnumber], mindist_result);
            pthread_mutex_unlock(&lock_queue[*tnumber]);
            *tnumber = (*tnumber + 1) % N_PQUEUE;
            added_tree_node++;
        } else {
            if (node->left_child->isax_cardinalities != NULL) {
                insert_tree_node_m_hybridpqueue_time(paa, node->left_child, index, bsf, pq, lock_queue, tnumber,
                                                     time_lb);
            }
            if (node->right_child->isax_cardinalities != NULL) {
                insert_tree_node_m_hybridpqueue_time(paa, node->right_child, index, bsf, pq, lock_queue, tnumber,
                                                     time_lb);
            }
        }
    }
}
