
#ifndef al_parallel_inmemory_query_engine_h
#define al_parallel_inmemory_query_engine_h
#include "config.h"
#include "../../globals.h"
#include "sax/ts.h"
#include "sax/sax.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "isax_index.h"
#include "isax_query_engine.h"
#include "parallel_query_engine.h"
#include "isax_node.h"
#include "pqueue.h"
#include "isax_first_buffer_layer.h"
#include "ads/isax_node_split.h"
#include "inmemory_index_engine.h"

typedef struct localStack {
    isax_node **val; 
    int top;
    int bottom;
}localStack;

typedef struct MESSI_workerdata
{
	isax_node *current_root_node;
	ts_type *paa,*paaU,*paaL,*ts,*uo,*lo;
	pqueue_t *pq;
	isax_index *index;
	float minimum_distance;
	int limit;
	pthread_mutex_t *lock_current_root_node;
	pthread_mutex_t *lock_queue;
    pthread_barrier_t *lock_barrier;
	pthread_rwlock_t *lock_bsf;
	query_result *bsf_result;
	int *node_counter;
	isax_node **nodelist;
	int amountnode;
	localStack *localstk; 
	localStack *allstk;
	pthread_mutex_t *locallock,*alllock;
	int *queuelabel,*allqueuelabel;
	pqueue_t **allpq;
	int startqueuenumber;
	int warpWind;
	pqueue_bsf *pq_bsf;
}MESSI_workerdata;

query_result exact_search_parads_inmemory (ts_type *ts, ts_type *paa, isax_index *index,
                           float minimum_distance, int min_checked_leaves);
query_result exact_search_serial_ParIS_inmemory(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves);
pqueue_bsf exact_topk_serial_ParIS_inmemory(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves,int k);

query_result exact_search_MESSI (ts_type *ts, ts_type *paa, isax_index *index, node_list *nodelist,
                                 float minimum_distance, int min_checked_leaves, int kn);
void insert_tree_node_m_hybridpqueue(float *paa,isax_node *node,isax_index *index,float bsf,pqueue_t **pq,pthread_mutex_t *lock_queue,int *tnumber);
extern int N_PQUEUE;

#endif
