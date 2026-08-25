
#ifndef al_inmemory_index_engine_h
#define al_inmemory_index_engine_h
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

void index_creation_pRecBuf(const char *ifilename, long int ts_num, int filetype_int, int apply_znorm,
                            isax_index *index, int kn);
root_mask_type isax_pRecBuf_index_insert_inmemory(isax_index *index, 
                                    sax_type * sax,
                                    file_position_type * pos, pthread_mutex_t *lock_firstnode,
                                    int workernumber, int total_workernumber, int kn);
isax_node * add_record_to_node_inmemory(isax_index *index,
                                 isax_node *tree_node,
                                 isax_node_record *record,
                                 const char leaf_size_check, int kn);
enum response flush_subtree_leaf_buffers_inmemory (isax_index *index, isax_node *node);
isax_index * isax_index_init_inmemory(isax_index_settings *settings);

typedef struct buffer_data_inmemory
{
	isax_index *index;
	int start_number,stop_number;	
	ts_type * ts;
	pthread_mutex_t *lock_firstnode;
	pthread_mutex_t *lock_fft_plan;
	int workernumber;
	int total_workernumber;
	pthread_barrier_t *lock_barrier1, *lock_barrier2;
	int *node_counter;
	int kn;
}buffer_data_inmemory;

extern ts_type *rawfile;

typedef struct node_list
{
	isax_node **nlist;
	int node_amount;
}node_list;
#endif
