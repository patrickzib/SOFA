//
//  isax_query_engine.c
//  al_isax
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//

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
#include <stdbool.h>

#include "ads/isax_query_engine.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/pqueue.h"
#include "ads/sax/sax.h"
#include "ads/isax_node_split.h"

#define NTHREADS 4

static float *mindists;

struct mindist_worker_data {
    unsigned long from;
    unsigned long to;
    ts_type *paa;
    isax_index *index;
};


static void *compute_mindists(void *ptr) {
    struct mindist_worker_data *arguments = ptr;
    unsigned long i;

    for(i=arguments->from; i<arguments->to; i++) {
        sax_type *sax = &(arguments->index->sax_cache[i * arguments->index->settings->n_segments]);
        mindists[i] = minidist_paa_to_isax(arguments->paa, sax,
                                               arguments->index->settings->max_sax_cardinalities,
                                               arguments->index->settings, 1);
    }
    return NULL;
}


query_result exact_search_serial(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves) 
{
    
    RESET_BYTES_ACCESSED
    mindists=malloc(sizeof(float) * index->sax_cache_size);
    unsigned long j;
    for (j = 0; j < index->sax_cache_size; j++)
        mindists[j] = FLT_MAX;
    
    
    query_result approximate_result = approximate_search(ts, paa, index);
    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }
    
    if(approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }
    
    
    
    unsigned long i;
    COUNT_INPUT_TIME_START
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = malloc(index->settings->ts_byte_size);
    COUNT_INPUT_TIME_END
    
    SET_APPROXIMATE(approximate_result.distance);

    pthread_t thread[NTHREADS];
    struct mindist_worker_data arguments[NTHREADS];
    for(i=0; i<NTHREADS; i++) {
        arguments[i].from = i*(index->sax_cache_size / NTHREADS);
        if(i < (NTHREADS-1)) {
            arguments[i].to = (i+1)*(index->sax_cache_size / NTHREADS);
        }
        else {
            arguments[i].to = index->sax_cache_size;
        }
        arguments[i].paa = paa;
        arguments[i].index = index;
        pthread_create(&thread[i], NULL, compute_mindists, &arguments[i]);
    }
    
    for(i=0; i<NTHREADS;i++) {
        pthread_join(thread[i], NULL);
    }
    for(i=0; i<index->sax_cache_size; i++) {
        if(mindists[i] <= approximate_result.distance) {
            COUNT_INPUT_TIME_START
            fseek(raw_file, i * index->settings->ts_byte_size, SEEK_SET);
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
            COUNT_INPUT_TIME_END
            float dist = ts_euclidean_distance(ts, ts_buffer, index->settings->timeseries_size, approximate_result.distance);
            if(dist < approximate_result.distance) {

                approximate_result.distance = dist;

    #ifdef STORE_ANSWER
                memcpy(index->answer, ts_buffer, index->settings->timeseries_size * sizeof(ts_type));
    #endif
            }
            INCREASE_BYTES_ACCESSED(index->settings->ts_byte_size)
        }
    }
    free(ts_buffer);
    fclose(raw_file);
    free(mindists);
    
    return approximate_result;
}
pqueue_bsf exact_topk_serial(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves, int k) 
{
    
    RESET_BYTES_ACCESSED
    pqueue_bsf *pq_bsf= pqueue_bsf_init(k);
    mindists=malloc(sizeof(float) * index->sax_cache_size);
    unsigned long j;
    for (j = 0; j < index->sax_cache_size; j++)
        mindists[j] = FLT_MAX;
    
    
    approximate_topk(ts, paa, index,pq_bsf);
    
    
    // Early termination...
    if (pq_bsf->knn[k-1] == 0) {
        return *pq_bsf;
    }
    if(pq_bsf->knn[k-1] == FLT_MAX  || min_checked_leaves > 1) {
        refine_topk_answer(ts, paa, index, pq_bsf, minimum_distance, min_checked_leaves);
    }
    unsigned long i;
    COUNT_INPUT_TIME_START
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = malloc(index->settings->ts_byte_size);
    COUNT_INPUT_TIME_END
    
    SET_APPROXIMATE(pq_bsf->knn[k-1]);

    pthread_t thread[NTHREADS];
    struct mindist_worker_data arguments[NTHREADS];
    for(i=0; i<NTHREADS; i++) {
        arguments[i].from = i*(index->sax_cache_size / NTHREADS);
        if(i < (NTHREADS-1)) {
            arguments[i].to = (i+1)*(index->sax_cache_size / NTHREADS);
        }
        else {
            arguments[i].to = index->sax_cache_size;
        }
        arguments[i].paa = paa;
        arguments[i].index = index;
        pthread_create(&thread[i], NULL, compute_mindists, &arguments[i]);
    }
    
    for(i=0; i<NTHREADS;i++) {
        pthread_join(thread[i], NULL);
    }
    for(i=0; i<index->sax_cache_size; i++) {
        if(mindists[i] <= pq_bsf->knn[k-1]) {
            COUNT_INPUT_TIME_START
            fseek(raw_file, i * index->settings->ts_byte_size, SEEK_SET);
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
            COUNT_INPUT_TIME_END
            float dist = ts_euclidean_distance(ts, ts_buffer, index->settings->timeseries_size, pq_bsf->knn[k-1]);
                if (dist <= pq_bsf->knn[pq_bsf->k-1]) {
                #ifdef STORE_ANSWER
                memcpy(index->answer, ts_buffer, index->settings->timeseries_size * sizeof(ts_type));
                #endif
            
                pqueue_bsf_insert(pq_bsf,dist,i,NULL);
            }
            INCREASE_BYTES_ACCESSED(index->settings->ts_byte_size)
        }
    }
    free(ts_buffer);
    fclose(raw_file);
    free(mindists);
    return *pq_bsf;
}

query_result  approximate_search (ts_type *ts, ts_type *paa, isax_index *index) 
{
    query_result result;

    sax_type *sax = malloc(sizeof(sax_type) * index->settings->n_segments);
    sax_from_paa(paa, sax, index->settings);

    root_mask_type root_mask = 0;
    CREATE_MASK(root_mask, index, sax);

    if (index->fbl->soft_buffers[(int) root_mask].initialized) {
        isax_node *node = index->fbl->soft_buffers[(int) root_mask].node;
        if (node->is_leaf && !node->has_full_data_file &&
            (node->leaf_size > index->settings->min_leaf_size))
        {
            split_node(index, node, 0, 1);
        }

        while (!node->is_leaf) {
            int location = index->settings->sax_bit_cardinality - 1 -
            node->split_data->split_mask[node->split_data->splitpoint];
            root_mask_type mask = index->settings->bit_masks[location];

            if(sax[node->split_data->splitpoint] & mask)
            {
                node = node->right_child;
            }
            else
            {
                node = node->left_child;
            }

            if (node->is_leaf && !node->has_full_data_file &&
                (node->leaf_size > index->settings->min_leaf_size))
            {
                split_node(index, node, 0, 1);
            }
        }

        result.distance = calculate_node_distance(index, node, ts, FLT_MAX);


        result.node = node;
    }
    else {
        result.node = NULL;
        result.distance = FLT_MAX;
    }

    free(sax);

    return result;
}

void  approximate_topk (ts_type *ts, ts_type *paa, isax_index *index, pqueue_bsf *pq_bsf) 
{

    sax_type *sax = malloc(sizeof(sax_type) * index->settings->n_segments);
    sax_from_paa(paa, sax, index->settings);

    root_mask_type root_mask = 0;
    CREATE_MASK(root_mask, index, sax);

    if (index->fbl->soft_buffers[(int) root_mask].initialized) {
        isax_node *node = index->fbl->soft_buffers[(int) root_mask].node;
        if (node->is_leaf && !node->has_full_data_file &&
            (node->leaf_size > index->settings->min_leaf_size))
        {
            split_node(index, node, 0, 1);
        }

        while (!node->is_leaf) {
            int location = index->settings->sax_bit_cardinality - 1 -
            node->split_data->split_mask[node->split_data->splitpoint];
            root_mask_type mask = index->settings->bit_masks[location];

            if(sax[node->split_data->splitpoint] & mask)
            {
                node = node->right_child;
            }
            else
            {
                node = node->left_child;
            }

            if (node->is_leaf && !node->has_full_data_file &&
                (node->leaf_size > index->settings->min_leaf_size))
            {
                split_node(index, node, 0, 1);
            }
        }

        calculate_node_topk(index, node, ts, pq_bsf);
    }
    for (int i = 0; i < pq_bsf->k-1; ++i)
    {
        pq_bsf->knn[i]=pq_bsf->knn[pq_bsf->k-1];
    }
    free(sax);
}
query_result refine_answer (ts_type *ts, ts_type *paa, isax_index *index, 
              query_result approximate_bsf_result, 
                            float minimum_distance, int limit)  
{ 
  query_result bsf_result = approximate_bsf_result; 
 
   int tight_bound = index->settings->tight_bound; 
  int aggressive_check = index->settings->aggressive_check; 
 
    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size, 
                               cmp_pri, get_pri, set_pri, get_pos, set_pos); 
 
    // Insert all root nodes in heap. 
    isax_node *current_root_node = index->first_node; 
    while (current_root_node != NULL) { 
    query_result * mindist_result = malloc(sizeof(query_result)); 
    mindist_result->distance =  minidist_paa_to_isax(paa, current_root_node->isax_values, current_root_node->isax_cardinalities, index->settings, 0); 
    mindist_result->node = current_root_node; 
        pqueue_insert(pq, mindist_result); 
        current_root_node = current_root_node->next; 
    } 
    query_result * n; 
    int checks = 0; 
    while ((n = pqueue_pop(pq))) 
    { 
    // The best node has a worse mindist, so search is finished! 
        if (n->distance >= bsf_result.distance || n->distance > minimum_distance) { 
            pqueue_insert(pq, n); 
            break; 
        } 
        else { 
          // If it is a leaf, check its real distance. 
            if (n->node->is_leaf) { 
        // *** ADAPTIVE SPLITTING *** 
            if (!n->node->has_full_data_file && 
          (n->node->leaf_size > index->settings->min_leaf_size)) 
            { 
          // Split and push again in the queue 
                split_node(index, n->node, 0, 1); 
          pqueue_insert(pq, n); 
                continue; 
            } 
        // *** EXTRA BOUNDING *** 
        if(tight_bound) { 
          float mindistance = calculate_minimum_distance(index, n->node, ts, paa); 
          if(mindistance >= bsf_result.distance) 
          { 
            free(n); 
            continue; 
          } 
        } 
        // *** REAL DISTANCE *** 
        checks++; 
        float distance = calculate_node_distance(index, n->node, ts, bsf_result.distance); 
        if (distance < bsf_result.distance) 
        { 
          bsf_result.distance = distance; 
          bsf_result.node = n->node; 
        } 
        if(checks > limit) { 
          pqueue_insert(pq, n); 
          break; 
        } 
      } 
            else { 
              // If it is an intermediate node calculate mindist for children 
              // and push them in the queue 
                if (n->node->left_child->isax_cardinalities != NULL) { 
          if(n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){ 
            float distance = calculate_node_distance(index, n->node->left_child, ts, bsf_result.distance); 
            if (distance < bsf_result.distance) 
            { 
              bsf_result.distance = distance; 
              bsf_result.node = n->node->left_child; 
            } 
          } 
          else { 
                    query_result * mindist_result = malloc(sizeof(query_result)); 
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->left_child->isax_values, n->node->left_child->isax_cardinalities, index->settings, 0); 
          mindist_result->node = n->node->left_child; 
                    pqueue_insert(pq, mindist_result); 
          } 
                } 
                if (n->node->right_child->isax_cardinalities != NULL) { 
          if(n->node->right_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){ 
            float distance = calculate_node_distance(index, n->node->right_child, ts, bsf_result.distance); 
            if (distance < bsf_result.distance) 
            { 
              bsf_result.distance = distance; 
              bsf_result.node = n->node->right_child; 
            } 
          } 
          else { 
                    query_result * mindist_result = malloc(sizeof(query_result)); 
          mindist_result->distance =  minidist_paa_to_isax(paa, n->node->right_child->isax_values, n->node->right_child->isax_cardinalities, index->settings, 0); 
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
    while ((n = pqueue_pop(pq))) 
    { 
        free(n); 
    } 
    // Free the priority queue. 
    pqueue_free(pq); 
   return bsf_result; 
} 
void refine_topk_answer (ts_type *ts, ts_type *paa, isax_index *index, 
              pqueue_bsf *pq_bsf, 
                            float minimum_distance, int limit)  
{  
   int tight_bound = index->settings->tight_bound; 
  int aggressive_check = index->settings->aggressive_check; 
 
    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size, 
                               cmp_pri, get_pri, set_pri, get_pos, set_pos); 
 
    // Insert all root nodes in heap. 
    isax_node *current_root_node = index->first_node; 
    while (current_root_node != NULL) { 
    query_result * mindist_result = malloc(sizeof(query_result)); 
    mindist_result->distance =  minidist_paa_to_isax(paa, current_root_node->isax_values, current_root_node->isax_cardinalities, index->settings, 0); 
    mindist_result->node = current_root_node; 
        if (mindist_result->distance < pq_bsf->knn[pq_bsf->k-1])
        {
            pqueue_insert(pq, mindist_result);
        }
        current_root_node = current_root_node->next; 
    } 
    query_result * n; 
    int checks = 0; 
    while ((n = pqueue_pop(pq))) 
    { 
    // The best node has a worse mindist, so search is finished! 
        if (n->distance >= pq_bsf->knn[pq_bsf->k-1] || n->distance > minimum_distance) {
            pqueue_insert(pq, n);
            break;
        } 
        else { 
          // If it is a leaf, check its real distance. 
            if (n->node->is_leaf) { 
        // *** ADAPTIVE SPLITTING *** 
            if (!n->node->has_full_data_file && 
          (n->node->leaf_size > index->settings->min_leaf_size)) 
            { 
          // Split and push again in the queue 
                split_node(index, n->node, 0, 1); 
          pqueue_insert(pq, n); 
                continue; 
            } 
        // *** EXTRA BOUNDING *** 
        if(tight_bound) { 
          float mindistance = calculate_minimum_distance(index, n->node, ts, paa); 
                    if(mindistance >= pq_bsf->knn[pq_bsf->k-1])
                    {
                        free(n);
                        continue;
                    }
        } 
        // *** REAL DISTANCE *** 
        checks++; 
        calculate_node_topk(index, n->node, ts, pq_bsf);
 
                if(pq_bsf->knn[pq_bsf->k-1] < FLT_MAX) {
                    pqueue_insert(pq, n);
                    break;
                }
      } 
            else { 
              // If it is an intermediate node calculate mindist for children 
              // and push them in the queue 
                if (n->node->left_child->isax_cardinalities != NULL) { 
          if(n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){ 
            calculate_node_topk(index, n->node->left_child, ts, pq_bsf);

          } 
          else { 
                    query_result * mindist_result = malloc(sizeof(query_result)); 
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->left_child->isax_values, n->node->left_child->isax_cardinalities, index->settings, 0); 
          mindist_result->node = n->node->left_child;  
            if (mindist_result->distance < pq_bsf->knn[pq_bsf->k-1])
            {
                pqueue_insert(pq, mindist_result);
            }
          } 
                } 
                if (n->node->right_child->isax_cardinalities != NULL) { 
          if(n->node->right_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){ 
            calculate_node_topk(index, n->node->right_child, ts, pq_bsf);
          } 
          else { 
                    query_result * mindist_result = malloc(sizeof(query_result)); 
          mindist_result->distance =  minidist_paa_to_isax(paa, n->node->right_child->isax_values, n->node->right_child->isax_cardinalities, index->settings, 0); 
                    mindist_result->node = n->node->right_child; 
            if (mindist_result->distance < pq_bsf->knn[pq_bsf->k-1])
            {
                pqueue_insert(pq, mindist_result);
            } 
          } 
                } 
            } 
 
            // Free the node currently popped. 
           free(n); 
        } 
    } 
    // Free the nodes that where not popped. 
    while ((n = pqueue_pop(pq))) 
    { 
        free(n); 
    } 
    // Free the priority queue. 
    for (int i = 0; i < pq_bsf->k-1; ++i)
    {
        pq_bsf->knn[i]=pq_bsf->knn[pq_bsf->k-1];
    }
    pqueue_free(pq); 
} 
query_result  approximate_search_SIMD (ts_type *ts, ts_type *paa, isax_index *index) 
{
    query_result result;

    sax_type *sax = malloc(sizeof(sax_type) * index->settings->n_segments);
    sax_from_paa(paa, sax, index->settings);

    root_mask_type root_mask = 0;
    CREATE_MASK(root_mask, index, sax);

    if (index->fbl->soft_buffers[(int) root_mask].initialized) {
        isax_node *node = index->fbl->soft_buffers[(int) root_mask].node;
        if (node->is_leaf && !node->has_full_data_file &&
            (node->leaf_size > index->settings->min_leaf_size))
        {
            split_node(index, node, 0, 1);
        }

        while (!node->is_leaf) {
            int location = index->settings->sax_bit_cardinality - 1 -
            node->split_data->split_mask[node->split_data->splitpoint];
            root_mask_type mask = index->settings->bit_masks[location];

            if(sax[node->split_data->splitpoint] & mask)
            {
                node = node->right_child;
            }
            else
            {
                node = node->left_child;
            }

            if (node->is_leaf && !node->has_full_data_file &&
                (node->leaf_size > index->settings->min_leaf_size))
            {
                split_node(index, node, 0, 1);
            }
        }

        result.distance = calculate_node_distance(index, node, ts, FLT_MAX);
        result.node = node;
    }
    else {
        result.node = NULL;
        result.distance = FLT_MAX;
    }

    free(sax);
    return result;
}


query_result exact_search (ts_type *ts, ts_type *paa, isax_index *index,
                           float minimum_distance, int min_checked_leaves) 
{
    query_result approximate_result = approximate_search(ts, paa, index);
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;

    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }
    if(approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }
    query_result bsf_result = approximate_result;
    COUNT_QUEUE_TIME_START
    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size,
                               cmp_pri, get_pri, set_pri, get_pos, set_pos);
    COUNT_QUEUE_TIME_END


    query_result *do_not_remove = &approximate_result;

    SET_APPROXIMATE(approximate_result.distance);

    RESET_BYTES_ACCESSED

    if(approximate_result.node != NULL) {
        // Insert approximate result in heap.
        COUNT_QUEUE_TIME_START
        pqueue_insert(pq, &approximate_result);
        COUNT_QUEUE_TIME_END
    }

    // Insert all root nodes in heap.
    isax_node *current_root_node = index->first_node;
    while (current_root_node != NULL) {
        query_result * mindist_result = malloc(sizeof(query_result));
        mindist_result->distance =  minidist_paa_to_isax(paa, current_root_node->isax_values, current_root_node->isax_cardinalities, index->settings, 0);
        mindist_result->node = current_root_node;
        COUNT_QUEUE_TIME_START
        pqueue_insert(pq, mindist_result);
        COUNT_QUEUE_TIME_END
        current_root_node = current_root_node->next;
    }
    query_result * n;
    int checks = 0;
    while ((n = pqueue_pop(pq)))
    {
        // The best node has a worse mindist, so search is finished!
        if (n->distance >= bsf_result.distance || n->distance > minimum_distance) {
            COUNT_QUEUE_TIME_START
            pqueue_insert(pq, n);
            COUNT_QUEUE_TIME_END
            break;
        }
        else {
            // If it is a leaf, check its real distance.
            if (n->node->is_leaf) {
                // *** ADAPTIVE SPLITTING ***
                if (!n->node->has_full_data_file &&
                    (n->node->leaf_size > index->settings->min_leaf_size))
                {
                    // Split and push again in the queue
                    split_node(index, n->node, 0, 1);
                    COUNT_QUEUE_TIME_START
                    pqueue_insert(pq, n);
                    COUNT_QUEUE_TIME_END
                    continue;
                }
                // *** EXTRA BOUNDING ***
                if(tight_bound) {
                    COUNT_CAL_TIME_START
                    float mindistance = calculate_minimum_distance(index, n->node, ts, paa);
                    COUNT_CAL_TIME_END
                    if(mindistance >= bsf_result.distance)
                    {
                        if(n != do_not_remove)
                            free(n);
                        continue;
                    }
                }
                // *** REAL DISTANCE ***
                checks++;

                COUNT_CAL_TIME_START
                float distance = calculate_node_distance(index, n->node, ts, bsf_result.distance);
                COUNT_CAL_TIME_END
                if (distance < bsf_result.distance)
                {
                    bsf_result.distance = distance;
                    bsf_result.node = n->node;
                }
            }
            else {
                // If it is an intermediate node calculate mindist for children
                // and push them in the queue
                if (n->node->left_child->isax_cardinalities != NULL) {
                    if(n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){
                        COUNT_CAL_TIME_START
                        float distance = calculate_node_distance(index, n->node->left_child, ts, bsf_result.distance);
                        COUNT_CAL_TIME_END
                        if (distance < bsf_result.distance)
                        {
                            bsf_result.distance = distance;
                            bsf_result.node = n->node->left_child;
                        }
                    }
                    else {
                    query_result * mindist_result = malloc(sizeof(query_result));
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->left_child->isax_values, n->node->left_child->isax_cardinalities, index->settings, 0);
                    mindist_result->node = n->node->left_child;
                    COUNT_QUEUE_TIME_START
                    pqueue_insert(pq, mindist_result);
                    COUNT_QUEUE_TIME_END
                    }
                }
                if (n->node->right_child->isax_cardinalities != NULL) {
                    if(n->node->right_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){
                        COUNT_CAL_TIME_START
                        float distance = calculate_node_distance(index, n->node->right_child, ts, bsf_result.distance);
                        COUNT_CAL_TIME_END
                        if (distance < bsf_result.distance)
                        {
                            bsf_result.distance = distance;
                            bsf_result.node = n->node->right_child;
                        }
                    }
                    else {
                    query_result * mindist_result = malloc(sizeof(query_result));
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->right_child->isax_values, n->node->right_child->isax_cardinalities, index->settings, 0);
                    mindist_result->node = n->node->right_child;
                    COUNT_QUEUE_TIME_START
                    pqueue_insert(pq, mindist_result);
                    COUNT_QUEUE_TIME_END
                    }
                }
            }

            // Free the node currently popped.
            if(n != do_not_remove)
                free(n);
        }
    }

    // Free the nodes that where not popped.
    while ((n = pqueue_pop(pq)))
    {
        if(n != do_not_remove)
            free(n);
    }
    // Free the priority queue.
    COUNT_QUEUE_TIME_START
    pqueue_free(pq);
    COUNT_QUEUE_TIME_END
   return bsf_result;
}

pqueue_bsf exact_topk (ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves, int k)
{
    pqueue_bsf *pq_bsf = pqueue_bsf_init(k);
    approximate_topk(ts, paa, index, pq_bsf);
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;

    // Early termination...
    if (pq_bsf->knn[k-1] == 0) {
        return *pq_bsf;
    }
    if(pq_bsf->knn[k-1] == FLT_MAX || min_checked_leaves > 1) {
        refine_topk_answer(ts, paa, index, pq_bsf, minimum_distance, min_checked_leaves);
    }

    COUNT_QUEUE_TIME_START
    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size, cmp_pri, get_pri, set_pri, get_pos, set_pos);
    COUNT_QUEUE_TIME_END

    SET_APPROXIMATE(pq_bsf->knn[k-1]);

    RESET_BYTES_ACCESSED

    // Insert all root nodes in heap.
    isax_node *current_root_node = index->first_node;
    while (current_root_node != NULL) {
        query_result * mindist_result = malloc(sizeof(query_result));
        mindist_result->distance =  minidist_paa_to_isax(paa, current_root_node->isax_values, current_root_node->isax_cardinalities, index->settings, 0);
        mindist_result->node = current_root_node;
        if (mindist_result->distance < pq_bsf->knn[k-1]) {
            COUNT_QUEUE_TIME_START
            pqueue_insert(pq, mindist_result);
            COUNT_QUEUE_TIME_END
        } else {
            free(mindist_result);
        }
        current_root_node = current_root_node->next;
    }
    query_result * n;
    int checks = 0;
    while ((n = pqueue_pop(pq)))
    {
        // The best node has a worse mindist, so search is finished!
        if (n->distance >= pq_bsf->knn[pq_bsf->k-1] || n->distance > minimum_distance) {
            COUNT_QUEUE_TIME_START
            pqueue_insert(pq, n);
            COUNT_QUEUE_TIME_END
            break;
        }
        else {
            // If it is a leaf, check its real distance.
            if (n->node->is_leaf) {
                // *** ADAPTIVE SPLITTING ***
                if (!n->node->has_full_data_file &&
                    (n->node->leaf_size > index->settings->min_leaf_size))
                {
                    // Split and push again in the queue
                    
                    split_node(index, n->node, 0, 1);
                    COUNT_QUEUE_TIME_START
                    pqueue_insert(pq, n);
                    COUNT_QUEUE_TIME_END
                    continue;
                }
                // *** EXTRA BOUNDING ***
                if(tight_bound) {
                    COUNT_CAL_TIME_START
                    float mindistance = calculate_minimum_distance(index, n->node, ts, paa);
                    COUNT_CAL_TIME_END
                    if(mindistance >= pq_bsf->knn[pq_bsf->k-1])
                    {
                        free(n);
                        continue;
                    }
                }
                // *** REAL DISTANCE ***
                checks++;
                calculate_node_topk(index, n->node, ts, pq_bsf);

            }
            else {
                // If it is an intermediate node calculate mindist for children
                // and push them in the queue
                if (n->node->left_child->isax_cardinalities != NULL) {
                    if(n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){
                        calculate_node_topk(index, n->node->left_child, ts, pq_bsf);
                    }
                    else {
                    query_result * mindist_result = malloc(sizeof(query_result));
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->left_child->isax_values, n->node->left_child->isax_cardinalities, index->settings, 0);
                    mindist_result->node = n->node->left_child;
                    COUNT_QUEUE_TIME_START
                    pqueue_insert(pq, mindist_result);
                    COUNT_QUEUE_TIME_END
                    }
                }
                if (n->node->right_child->isax_cardinalities != NULL) {
                    if(n->node->right_child->is_leaf && !n->node->right_child->has_partial_data_file && aggressive_check){
                        calculate_node_topk(index, n->node->right_child, ts, pq_bsf);
                    }
                    else {
                    query_result * mindist_result = malloc(sizeof(query_result));
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->right_child->isax_values, n->node->right_child->isax_cardinalities, index->settings, 0);
                    mindist_result->node = n->node->right_child;
                    COUNT_QUEUE_TIME_START
                    pqueue_insert(pq, mindist_result);
                    COUNT_QUEUE_TIME_END
                    }
                }
            }
        }
    }

    while ((n = pqueue_pop(pq)))
    {
        free(n);
    }
    // Free the priority queue.
    COUNT_QUEUE_TIME_START
    pqueue_free(pq);
    COUNT_QUEUE_TIME_END
    return *pq_bsf;
}
