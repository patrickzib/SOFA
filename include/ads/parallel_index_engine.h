#ifndef parallel_parallel_index_engine_h
#define parallel_parallel_index_engine_h
#include "config.h"
#include "../../globals.h"
#include "sax/ts.h"
#include "sax/sax.h"	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include "isax_index.h"
#include "isax_query_engine.h"

void isax_index_binary_file_m(const char *ifilename, int ts_num, isax_index *index,int calculate_thread);
void isax_index_binary_file_m_new(const char *ifilename, int ts_num, isax_index *index,int calculate_thread);
isax_node * insert_to_pRecBuf(parallel_first_buffer_layer *fbl,
									sax_type *sax,
									file_position_type *pos,
									root_mask_type mask, 
									isax_index *index, 
									pthread_mutex_t *lock_firstnode,
									int workernumber,
									int total_workernumber);
extern int read_block_length;
#endif
