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

#include <time.h>
#include <sys/time.h>

#include "ads/isax_query_engine.h"
#include "ads/inmemory_query_engine.h"
#include "ads/parallel_query_engine.h"
#include "ads/parallel_inmemory_query_engine.h"
#include "ads/parallel_index_engine.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/pqueue.h"
#include "ads/sax/sax.h"
#include "ads/isax_node_split.h"
#include "ads/sfa/dft.h"
#include "ads/sfa/sfa.h"
#include "ads/spartan/spartan.h"
#include "ads/calc_utils.h"

float *MINDISTS;

void *compute_mindists_in(void *ptr) {
    struct args_in *arguments = (struct args_in *) ptr;
    unsigned long i;

    for (i = arguments->from; i < arguments->to; i++) {
        sax_type *sax = &(arguments->index->sax_cache[i * arguments->index->settings->n_segments]);
        MINDISTS[i] = messi_minidist_raw(arguments->index, arguments->paa, sax,
                                         arguments->index->settings->max_sax_cardinalities,
                                         FLT_MAX);
    }

    return NULL;
}

query_result approximate_search_inmemory_pRecBuf(ts_type *ts, ts_type *paa, isax_index *index, int kn) {
    query_result result;

    sax_type *sax = malloc(sizeof(sax_type) * index->settings->n_segments);


    //SFA
    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        sfa_from_fft(index, paa, sax);
    } else if (index->settings->function_type == 5) {
        spartan_from_pca(index, paa, sax);
    } else {
        sax_from_paa(paa, sax, index->settings);
    }

    root_mask_type root_mask = 0;
    root_mask = isax_root_mask_from_sax(index, sax, kn);

    COUNT_INIT_TIME_END

    COUNT_TREE_PASS_TIME_START

    if ((&((parallel_first_buffer_layer * )(index->fbl))->soft_buffers[(int) root_mask])->initialized) {
        isax_node *node = (&((parallel_first_buffer_layer * )(index->fbl))->soft_buffers[(int) root_mask])->node;
        // Traverse tree

        // Adaptive splitting

        while (!node->is_leaf) {
            // TODO: confirm whether this should be "sax_bit_cardinality - 1 - split_mask".
            int location = index->settings->sax_bit_cardinality -
                           node->split_data->split_mask[node->split_data->splitpoint];
            root_mask_type mask = index->settings->bit_masks[location];

            if (sax[node->split_data->splitpoint] & mask) {
                node = node->right_child;
            } else {
                node = node->left_child;
            }

            // Adaptive splitting
        }

        result.distance = calculate_node_distance_inmemory(index, node, ts, paa, FLT_MAX);
        result.node = node;
    } else {
        result.node = NULL;
        result.distance = FLT_MAX;
    }
    COUNT_TREE_PASS_TIME_END

    free(sax);

    return result;
}

float calculate_node_distance_inmemory(isax_index *index, isax_node *node, ts_type *query, ts_type *paa, float bsf) {
    float distmin;
    COUNT_CHECKED_NODE()

    // If node has buffered data
    if (node->buffer != NULL) {
        int i;

        RDcalculationnumber = RDcalculationnumber + node->buffer->full_buffer_size;
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

        RDcalculationnumber = RDcalculationnumber + node->buffer->tmp_full_buffer_size;
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

        RDcalculationnumber = RDcalculationnumber + node->buffer->partial_buffer_size;
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
    return bsf;
}

//debugging only!!!
float calculate_node_distance2_inmemory(isax_index *index, isax_node *node, ts_type *query, ts_type *paa, float bsf) {
    float distmin;
    COUNT_CHECKED_NODE()

    // If node has buffered data
    if (node->buffer != NULL) {
        int i;
        int real_dist_counter = 0;

        for (i = 0; i < node->buffer->partial_buffer_size; i++) {
            distmin = messi_minidist_raw(
                index, paa, node->buffer->partial_sax_buffer[i],
                index->settings->max_sax_cardinalities, bsf);

            if (distmin < bsf) {
                float dist = ts_ed(
                    query, &(rawfile[*node->buffer->partial_position_buffer[i]]),
                    index->settings->timeseries_size, bsf);

                real_dist_counter++;

                if (dist < bsf) {
                    bsf = dist;
                }
            }
        }

        __sync_fetch_and_add(&LBDcalculationnumber, node->buffer->partial_buffer_size);
        __sync_fetch_and_add(&RDcalculationnumber, real_dist_counter);
    }
    return bsf;
}

query_result refine_answer_inmemory(ts_type *ts, ts_type *paa, isax_index *index,
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

        mindist_result->distance = minidist_paa_to_isax(paa, current_root_node->isax_values,
                                                        current_root_node->isax_cardinalities,
                                                        index->settings, 0);
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
                float distance = calculate_node_distance_inmemory(index, n->node, ts, paa, bsf_result.distance);
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
                        mindist_result->distance = minidist_paa_to_isax(paa, n->node->left_child->isax_values,
                                                                        n->node->left_child->isax_cardinalities,
                                                                        index->settings, 0);
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
                        mindist_result->distance = minidist_paa_to_isax(paa, n->node->right_child->isax_values,
                                                                        n->node->right_child->isax_cardinalities,
                                                                        index->settings, 0);
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


float calculate_minimum_distance_inmemory(isax_index *index, isax_node *node, ts_type *raw_query, ts_type *query) {
    float bsfLeaf = minidist_paa_to_isax(query, node->isax_values,
                                         node->isax_cardinalities,
                                         index->settings, 0);
    float bsfRecord = FLT_MAX;

    if (!index->has_wedges) {
        int i = 0;


        if (node->buffer != NULL) {
            for (i = 0; i < node->buffer->partial_buffer_size; i++) {
                float mindist = minidist_paa_to_isax_raw_SIMD(query, node->buffer->partial_sax_buffer[i],
                                                              index->settings->max_sax_cardinalities,
                                                              index->settings);
                if (mindist < bsfRecord) {
                    bsfRecord = mindist;
                }
            }

            for (i = 0; i < node->buffer->tmp_partial_buffer_size; i++) {
                float mindist = minidist_paa_to_isax_raw_SIMD(query, node->buffer->tmp_partial_sax_buffer[i],
                                                              index->settings->max_sax_cardinalities,
                                                              index->settings);
                if (mindist < bsfRecord) {
                    bsfRecord = mindist;
                }
            }
        }
    } else {
        int i = 0;
        if (node->wedges[0] == FLT_MIN) {
            bsfRecord = FLT_MAX;
        } else {
            bsfRecord = 0;
            ts_type *min_wedge = &node->wedges[0];
            ts_type *max_wedge = &node->wedges[index->settings->timeseries_size];
            if (raw_query[i] > max_wedge[i]) {
                bsfRecord += (raw_query[i] - max_wedge[i]) * (raw_query[i] - max_wedge[i]);
            } else if (raw_query[i] < max_wedge[i] && raw_query[i] > min_wedge[i]) {
            } else {
                bsfRecord += (min_wedge[i] - raw_query[i]) * (min_wedge[i] - raw_query[i]);
            }
        }

    }
    float bsf = (bsfRecord == FLT_MAX) ? bsfLeaf : bsfRecord;


    return bsf;
}
