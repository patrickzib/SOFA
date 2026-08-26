//
//  isax_query_engine.h
//  al_isax
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//

#ifndef al_isax_isax_query_engine_h
#define al_isax_isax_query_engine_h
#include "config.h"
#include "../../globals.h"
#include "isax_index.h"
#include "isax_node.h"
#include "pqueue.h"

typedef struct query_result {
    float distance;
    isax_node *node;
    size_t pqueue_position;
} query_result;


static inline int
cmp_pri(double next, double curr)
{
	return (next > curr);
}


static inline double
get_pri(void *a)
{
	return (double) ((query_result *) a)->distance;
}


static inline void
set_pri(void *a, double pri)
{
	((query_result *) a)->distance = (float)pri;
}


static inline size_t
get_pos(void *a)
{
	return ((query_result *) a)->pqueue_position;
}


static inline void
set_pos(void *a, size_t pos)
{
	((query_result *) a)->pqueue_position = pos;
}


/// VARIOUS QUERY TYPES
query_result exact_search_serial(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves);
pqueue_bsf exact_topk_serial(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves, int k);

query_result exact_search (ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves);
pqueue_bsf exact_topk (ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves, int k);

/// HELPE FUNCTIONS
query_result approximate_search (ts_type *ts, ts_type *paa, isax_index *index);

void  approximate_topk (ts_type *ts, ts_type *paa, isax_index *index, pqueue_bsf *pq_bsf);
query_result approximate_search_SIMD (ts_type *ts, ts_type *paa, isax_index *index);
query_result refine_answer (ts_type *ts, ts_type *paa, isax_index *index,
							query_result approximate_bsf_result,
							float minimum_distance, int limit);
void refine_topk_answer (ts_type *ts, ts_type *paa, isax_index *index, 
              pqueue_bsf *pq_bsf, 
                            float minimum_distance, int limit);


#endif
