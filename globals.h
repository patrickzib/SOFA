//
//  defines.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "config.h"
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include "include/ads/pthread_barrier.h"

#ifndef isax_globals_h
#define isax_globals_h

#define STORE_ANSWER


#pragma GCC diagnostic ignored "-Wunused-result" 
#pragma GCC diagnostic ignored "-Wunused-variable" 
 

#define PAGE_SIZE 4096
#define PROGRESS_CALCULATE_THREAD_NUMBER 12
#define PROGRESS_FLUSH_THREAD_NUMBER 12
#define QUERIES_THREAD_NUMBER 4
#define DISK_BUFFER_SIZE 8192
#define LOCK_SIZE 65536 
///// TYPES /////
typedef unsigned char sax_type;
typedef unsigned char file_type;
typedef float ts_type;
typedef unsigned long long file_position_type;
typedef unsigned long long root_mask_type;

/* Dense input encoding.  Keep byte input stored as unsigned bytes and apply
 * signedness only while converting to ts_type; this avoids duplicating every
 * I/O buffer and preserves the historical --filetype-int behavior. */
enum file_input_type {
    FILE_INPUT_FLOAT32 = 0,
    FILE_INPUT_UINT8 = 1,
    FILE_INPUT_INT8 = 2
};

static inline ts_type file_value_to_ts(file_type value, int input_type) {
    return input_type == FILE_INPUT_INT8
        ? (ts_type) (int8_t) value
        : (ts_type) value;
}

enum response {OUT_OF_MEMORY_FAILURE, FAILURE, SUCCESS};
enum insertion_mode {PARTIAL = 1, 
                     TMP = 2, 
                     FULL = 4,
                     NO_TMP = 8};

enum buffer_cleaning_mode {FULL_CLEAN, TMP_ONLY_CLEAN, TMP_AND_TS_CLEAN};
enum node_cleaning_mode {DO_NOT_INCLUDE_CHILDREN = 0,
                         INCLUDE_CHILDREN = 1};

///// DEFINITIONS /////
#define MINVAL -2000000
#define MAXVAL 2000000
#define DELIMITER ' '
#define TRUE 1
#define FALSE 0
#define BUFFER_REALLOCATION_RATE  2 

///// GLOBAL VARIABLES /////
int FLUSHES;
unsigned long BYTES_ACCESSED;
float APPROXIMATE;

void* LOGFILE;
/* Shared priority-queue count used by iSAX and trie query scheduling. */
extern int N_PQUEUE;
/* Console-only query-row sampling.  CSV logging and aggregate statistics are
 * unaffected.  MESSI.c initializes this to 10 unless overridden. */
extern int query_report_interval;

#define SHOULD_REPORT_QUERY(query_index, query_count) \
    (query_report_interval > 0 && \
     ((query_index) == 0 || (query_index) == (query_count) - 1 || \
      (((query_index) + 1) % query_report_interval == 0)))

#define INCREASE_BYTES_ACCESSED(new_bytes) \
	    BYTES_ACCESSED += (unsigned long) new_bytes;
#define RESET_BYTES_ACCESSED \
		BYTES_ACCESSED = 0;
#define SET_APPROXIMATE(approximate)\
		APPROXIMATE = approximate;

#define SET_LOGFILE(logfile) \
    LOGFILE = logfile;

///// MACROS /////
#define CREATE_MASK(mask, index, sax_array)\
	int mask__i; \
	for (mask__i=0; mask__i < index->settings->n_segments; mask__i++) \
		if(index->settings->bit_masks[index->settings->sax_bit_cardinality - 1] & sax_array[mask__i]) \
			mask |= index->settings->bit_masks[index->settings->n_segments - mask__i - 1];  

#define CREATE_MASK_SFAD(mask, index, sax_array, kn)\
	int mask__i2, mask__j; \
	for (mask__i2=0; mask__i2 < index->settings->n_segments/kn; mask__i2++) \
		for (mask__j=0; mask__j < kn; mask__j++) \
			if(index->settings->bit_masks[index->settings->sax_bit_cardinality - 1-mask__j] & sax_array[mask__i2]) \
				mask |= index->settings->bit_masks[index->settings->n_segments - mask__i2*kn-mask__j - 1];

///// BENCHMARKING /////
#ifdef BENCHMARK
		#include <time.h>
		#include <sys/time.h>
	   
        double tS;
        double tE;
        int added_tree_node;
        struct timeval total_time_start;
        struct timeval queue_time_start;
        struct timeval parse_time_start;
        struct timeval input_time_start;
        struct timeval input2_time_start;
        struct timeval indexing_time_start;
        struct timeval querying_time_start;
        struct timeval cal_time_start;
        struct timeval output_time_start;
        struct timeval output2_time_start;
        struct timeval load_node_start;
        struct timeval binning_time_start;
        struct timeval current_time;
        struct timeval fetch_start;
        struct timeval fetch_check_start;
        double total_input_time;
        double total_queue_time;
        double total_input2_time;
        double total_cal_time;
        double load_node_time;
        double total_output_time;
        double total_output2_time;
        double total_parse_time;
        double total_time;
        double total_binning_time;
        double total_indexing_time;
        double total_querying_time;

        double total_querying_time_all;
        unsigned long bytes_accessed_all;
        float approximate_all;
        double result_distance_all;
        double total_time_all;

        struct timeval init_time_start;
        struct timeval tree_pass_time_start;

        double total_init_time;
        double total_tree_pass_time;
        unsigned long int TOTAL_PQ_INSERT_TIME;
        unsigned long int TOTAL_PQ_REMOVE_TIME;
        unsigned long int TOTAL_LB_DIST_CALC_TIME;
        unsigned long int TOTAL_MBR_DIST_CALC_TIME;
        unsigned long int TOTAL_RECORD_LB_DIST_CALC_TIME;
        unsigned long int TOTAL_REAL_DIST_CALC_TIME;
        /* Opt-in, direct query-phase profiling.  These are accumulated worker
         * times (rather than wall-clock time), so parallel work can sum to
         * more than the query wall time. */
        unsigned long int TOTAL_TREE_TRAVERSAL_TIME;
        unsigned long int TOTAL_TRIE_FRONTIER_TIME;
        unsigned long int TOTAL_TRIE_QUEUE_TIME;
        unsigned long int TOTAL_TRIE_HEAP_TIME;
        unsigned long int TOTAL_TRIE_SYNC_TIME;
        int profile_query_phases;

        unsigned long int TOTAL_INDEXING_PART_TIME;
        unsigned long int TOTAL_TRANSFORMATION_PART_TIME;

        double total_init_time_all;
        double total_tree_pass_time_all;
        double total_pq_insert_time_all;
        double total_pq_remove_time_all;
        double total_lb_dist_calc_time_all;
        double total_mbr_dist_calc_time_all;
        double total_record_lb_dist_calc_time_all;
        double total_real_dist_calc_time_all;
        double total_trie_frontier_time_all;
        double total_trie_queue_time_all;
        double total_trie_heap_time_all;
        double total_trie_sync_time_all;

        int total_tree_nodes;
        int loaded_nodes;
        int checked_nodes;
        int loaded_nodes_all;
        int checked_nodes_all;
        int stats_header_printed;
        file_position_type loaded_records;
        unsigned long int LBDcalculationnumber;
        unsigned long int RDcalculationnumber;
        unsigned long int LBDcalculationnumber_all;
        unsigned long int RDcalculationnumber_all;
        /* Trie leaf IVF diagnostics; separate from record lower-bound
         * counts because a pruned cluster skips many record checks. */
        unsigned long int trie_cluster_bounds;
        unsigned long int trie_cluster_pruned;
        unsigned long int trie_cluster_records_pruned;
        unsigned long int trie_cluster_symbolic_pruned;
        unsigned long int trie_cluster_symbolic_records_pruned;
        unsigned long int trie_cluster_raw_ball_bounds;
        unsigned long int trie_cluster_raw_ball_pruned;
        unsigned long int trie_cluster_raw_ball_records_pruned;
        unsigned long int trie_radial_candidates;
        unsigned long int trie_radial_pruned;
        unsigned long int trie_cluster_bounds_all;
        unsigned long int trie_cluster_pruned_all;
        unsigned long int trie_cluster_records_pruned_all;
        unsigned long int trie_cluster_symbolic_pruned_all;
        unsigned long int trie_cluster_symbolic_records_pruned_all;
        unsigned long int trie_cluster_raw_ball_bounds_all;
        unsigned long int trie_cluster_raw_ball_pruned_all;
        unsigned long int trie_cluster_raw_ball_records_pruned_all;
        unsigned long int trie_radial_candidates_all;
        unsigned long int trie_radial_pruned_all;

        #define INIT_STATS() total_input_time = 0;\
                            total_output_time = 0;\
                            total_time = 0;\
                            total_parse_time = 0;\
                            total_tree_nodes = 0;\
                            total_binning_time = 0;\
                            total_indexing_time = 0;\
                            total_querying_time = 0;\
                            loaded_nodes = 0;\
                            checked_nodes = 0;\
                            loaded_nodes_all = 0;\
                            checked_nodes_all = 0;\
                            load_node_time = 0;\
                            loaded_records = 0; \
                            APPROXIMATE=0;\
                            BYTES_ACCESSED=0;\
                            LBDcalculationnumber=0;\
                            RDcalculationnumber=0;\
                            LBDcalculationnumber_all=0;\
                            RDcalculationnumber_all=0;\
                            trie_cluster_bounds=0;\
                            trie_cluster_pruned=0;\
                            trie_cluster_records_pruned=0;\
                            trie_cluster_symbolic_pruned=0;\
                            trie_cluster_symbolic_records_pruned=0;\
                            trie_cluster_raw_ball_bounds=0;\
                            trie_cluster_raw_ball_pruned=0;\
                            trie_cluster_raw_ball_records_pruned=0;\
                            trie_radial_candidates=0;\
                            trie_radial_pruned=0;\
                            trie_cluster_bounds_all=0;\
                            trie_cluster_pruned_all=0;\
                            trie_cluster_records_pruned_all=0;\
                            trie_cluster_symbolic_pruned_all=0;\
                            trie_cluster_symbolic_records_pruned_all=0;\
                            trie_cluster_raw_ball_bounds_all=0;\
                            trie_cluster_raw_ball_pruned_all=0;\
                            trie_cluster_raw_ball_records_pruned_all=0;\
                            trie_radial_candidates_all=0;\
                            trie_radial_pruned_all=0;\
                            total_querying_time_all=0.0;\
                            bytes_accessed_all=0;\
                            approximate_all = 0.0;\
                            result_distance_all=0.0;\
                            total_init_time=0.0;\
                            total_tree_pass_time=0.0;\
                            TOTAL_PQ_INSERT_TIME=0;\
                            TOTAL_PQ_REMOVE_TIME=0;\
                            TOTAL_LB_DIST_CALC_TIME=0;\
                            TOTAL_MBR_DIST_CALC_TIME=0;\
                            TOTAL_RECORD_LB_DIST_CALC_TIME=0;\
                            TOTAL_REAL_DIST_CALC_TIME=0;\
                            TOTAL_TREE_TRAVERSAL_TIME=0;\
                            TOTAL_TRIE_FRONTIER_TIME=0;\
                            TOTAL_TRIE_QUEUE_TIME=0;\
                            TOTAL_TRIE_HEAP_TIME=0;\
                            TOTAL_TRIE_SYNC_TIME=0;\
                            profile_query_phases=0;\
                            total_init_time_all=0.0;\
                            total_tree_pass_time_all=0.0;\
                            total_pq_insert_time_all=0.0;\
                            total_pq_remove_time_all=0.0;\
                            total_lb_dist_calc_time_all=0.0;\
                            total_mbr_dist_calc_time_all=0.0;\
                            total_record_lb_dist_calc_time_all=0.0;\
                            total_real_dist_calc_time_all=0.0;\
                            total_trie_frontier_time_all=0.0;\
                            total_trie_queue_time_all=0.0;\
                            total_trie_heap_time_all=0.0;\
                            total_trie_sync_time_all=0.0;\
                            TOTAL_INDEXING_PART_TIME = 0.0;\
                            TOTAL_TRANSFORMATION_PART_TIME = 0.0;\
                            stats_header_printed = 0;
        #define PRINT_STATS(result_distance) do { \
            if (!stats_header_printed) { \
                fprintf(stderr, "%4s %8s %13s %20s %12s %12s\n", \
                       "idx:", "nodes", "checked_nodes", "approximate_distance", "distance", "cumulative_ms"); \
                stats_header_printed = 1; \
            } \
            fprintf(stderr, "%8d %13d %20.3f %12.3f %12.3f\n", \
                   total_tree_nodes, checked_nodes, APPROXIMATE, \
                   result_distance, total_time / 1000.0); \
        } while (0);
        #define PRINT_STATS_HEADER() do { \
            if (!stats_header_printed) { \
                fprintf(stderr, "%4s %8s %13s %20s %12s %12s\n", \
                       "idx:", "nodes", "checked_nodes", "approximate_distance", "distance", "cumulative_ms"); \
                stats_header_printed = 1; \
            } \
        } while (0);
        //#define PRINT_STATS(result_distance) printf("%d\t",loaded_nodes);
        /* One self-contained index-construction record.  All durations are
         * microseconds.  In-memory benchmarks do not persist an index file,
         * so deliberately do not advertise a fictitious size. */
        #define INIT_INDEX_STATS_FILE(ifile) \
            fprintf((ifile), \
                    "layout,function type,indexed records,series length,symbolic dimensions,record lower-bound dimensions,binning us,index construction wall us,total construction wall us\n" \
                    "%s,%d,%ld,%d,%d,%d,%.0f,%.0f,%.0f\n", \
                    (idx)->settings->index_type == MESSI_INDEX_TRIE ? "trie" : "isax", \
                    (idx)->settings->function_type, (long) (idx)->total_records, \
                    (idx)->settings->timeseries_size, (idx)->settings->n_segments, \
                    (idx)->settings->trie_bound_dimensions, total_binning_time, \
                    total_indexing_time, total_binning_time + total_indexing_time)
        #define INIT_SAVE_FILE(ifile) fprintf(ifile, "querying time, init time, tree pass time, pq insert time, pq remove time, lb dist time, real dist time, nodes, checked_nodes, bytes_accessed, loaded_nodes, loaded_records, real dist calc, lb dist calc, approximate_distance, distance, total\n");
        #define SAVE_STATS(result_distance) fprintf(LOGFILE,"%.0f, %.0f, %.0f, %lu, %lu, %lu, %lu, %d, %d, %ld, %d, %lld, %lu, %lu, %lf, %lf, %.0f\n", \
        total_querying_time, total_init_time, \
        total_tree_pass_time, TOTAL_PQ_INSERT_TIME, TOTAL_PQ_REMOVE_TIME,\
        TOTAL_LB_DIST_CALC_TIME, TOTAL_REAL_DIST_CALC_TIME,\
        total_tree_nodes, checked_nodes, \
        BYTES_ACCESSED, loaded_nodes, \
        loaded_records, RDcalculationnumber, \
        LBDcalculationnumber, APPROXIMATE, \
        result_distance, total_time);\
        total_querying_time_all += total_querying_time;\
        total_time_all = total_time;\
        bytes_accessed_all += BYTES_ACCESSED;\
        approximate_all += APPROXIMATE;\
        result_distance_all += result_distance;\
        RDcalculationnumber_all += RDcalculationnumber;\
        RDcalculationnumber = 0;\
        LBDcalculationnumber_all += LBDcalculationnumber;\
        LBDcalculationnumber = 0;\
        trie_cluster_bounds_all += trie_cluster_bounds;\
        trie_cluster_bounds = 0;\
        trie_cluster_pruned_all += trie_cluster_pruned;\
        trie_cluster_pruned = 0;\
        trie_cluster_records_pruned_all += trie_cluster_records_pruned;\
        trie_cluster_records_pruned = 0;\
        trie_cluster_symbolic_pruned_all += trie_cluster_symbolic_pruned;\
        trie_cluster_symbolic_pruned = 0;\
        trie_cluster_symbolic_records_pruned_all += trie_cluster_symbolic_records_pruned;\
        trie_cluster_symbolic_records_pruned = 0;\
        trie_cluster_raw_ball_bounds_all += trie_cluster_raw_ball_bounds;\
        trie_cluster_raw_ball_bounds = 0;\
        trie_cluster_raw_ball_pruned_all += trie_cluster_raw_ball_pruned;\
        trie_cluster_raw_ball_pruned = 0;\
        trie_cluster_raw_ball_records_pruned_all += trie_cluster_raw_ball_records_pruned;\
        trie_cluster_raw_ball_records_pruned = 0;\
        trie_radial_candidates_all += trie_radial_candidates;\
        trie_radial_candidates = 0;\
        trie_radial_pruned_all += trie_radial_pruned;\
        trie_radial_pruned = 0;\
        checked_nodes_all += checked_nodes;\
        checked_nodes = 0;\
        loaded_nodes_all += loaded_nodes;\
        loaded_nodes = 0;\
        total_init_time_all += total_init_time;\
        total_init_time=0.0;\
        total_tree_pass_time_all += total_tree_pass_time;\
        total_tree_pass_time=0.0;\
        total_pq_insert_time_all += TOTAL_PQ_INSERT_TIME;\
        TOTAL_PQ_INSERT_TIME=0;\
        total_pq_remove_time_all += TOTAL_PQ_REMOVE_TIME;\
        TOTAL_PQ_REMOVE_TIME=0;\
        total_lb_dist_calc_time_all += TOTAL_LB_DIST_CALC_TIME;\
        TOTAL_LB_DIST_CALC_TIME=0;\
        total_mbr_dist_calc_time_all += TOTAL_MBR_DIST_CALC_TIME;\
        TOTAL_MBR_DIST_CALC_TIME=0;\
        total_record_lb_dist_calc_time_all += TOTAL_RECORD_LB_DIST_CALC_TIME;\
        TOTAL_RECORD_LB_DIST_CALC_TIME=0;\
        total_real_dist_calc_time_all += TOTAL_REAL_DIST_CALC_TIME;\
        TOTAL_REAL_DIST_CALC_TIME=0;\
        total_trie_frontier_time_all += TOTAL_TRIE_FRONTIER_TIME;\
        TOTAL_TRIE_FRONTIER_TIME=0;\
        total_trie_queue_time_all += TOTAL_TRIE_QUEUE_TIME;\
        TOTAL_TRIE_QUEUE_TIME=0;\
        total_trie_heap_time_all += TOTAL_TRIE_HEAP_TIME;\
        TOTAL_TRIE_HEAP_TIME=0;\
        total_trie_sync_time_all += TOTAL_TRIE_SYNC_TIME;\
        TOTAL_TRIE_SYNC_TIME=0;
        #define SAVE_STATS_TOTAL(ifile, queries_size) fprintf(ifile,"%lf, %lf, %lf, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", \
        (total_querying_time_all / queries_size), (total_init_time_all / queries_size),\
        (total_tree_pass_time_all / queries_size), (total_pq_insert_time_all / queries_size),(total_pq_remove_time_all / queries_size),\
        (total_lb_dist_calc_time_all / queries_size), (total_real_dist_calc_time_all / queries_size),\
        total_tree_nodes, ((double)checked_nodes_all / queries_size), \
        ((double)bytes_accessed_all / queries_size), ((double)loaded_nodes_all / queries_size), \
        ((double)loaded_records / queries_size), ((double)RDcalculationnumber_all / queries_size), \
        ((double)LBDcalculationnumber_all / queries_size), (approximate_all / queries_size),\
        (result_distance_all / queries_size), total_time_all);
        #define min(x,y)  ( x<y?x:y )
        #define COUNT_NEW_NODE() __sync_fetch_and_add(&total_tree_nodes,1);
        #define COUNT_LOADED_NODE() loaded_nodes++;
        #define COUNT_CHECKED_NODE() checked_nodes++;
        #define COUNT_LOADED_RECORD() loaded_records++;

        #define COUNT_INIT_TIME_START gettimeofday(&init_time_start, NULL);
        #define COUNT_TREE_PASS_TIME_START gettimeofday(&tree_pass_time_start, NULL);
        #define COUNT_PQ_INSERT_TIME_START gettimeofday(&pq_insert_time_start, NULL);

        #define COUNT_INIT_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = init_time_start.tv_sec*1000000 + (init_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_init_time += (tE - tS);
        #define COUNT_TREE_PASS_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = tree_pass_time_start.tv_sec*1000000 + (tree_pass_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_tree_pass_time += (tE - tS);                        

        #define COUNT_INPUT_TIME_START gettimeofday(&input_time_start, NULL);
        #define COUNT_QUEUE_TIME_START gettimeofday(&queue_time_start, NULL);
        #define COUNT_INDEXING_TIME_START gettimeofday(&indexing_time_start, NULL);
        #define COUNT_QUERYING_TIME_START gettimeofday(&querying_time_start, NULL); \
          total_querying_time = 0.0;
        #define COUNT_CAL_TIME_START gettimeofday(&cal_time_start, NULL); 
        #define COUNT_INPUT2_TIME_START gettimeofday(&input2_time_start, NULL); 
        #define COUNT_OUTPUT_TIME_START gettimeofday(&output_time_start, NULL); 
        #define COUNT_OUTPUT2_TIME_START gettimeofday(&output2_time_start, NULL);
        #define COUNT_TOTAL_TIME_START gettimeofday(&total_time_start, NULL);   
        #define COUNT_PARSE_TIME_START gettimeofday(&parse_time_start, NULL);   
        #define COUNT_LOAD_NODE_START gettimeofday(&load_node_start, NULL);
        #define COUNT_BINNING_TIME_START gettimeofday(&binning_time_start, NULL);
        #define COUNT_BINNING_TIME_END gettimeofday(&current_time, NULL); \
                                      tS = binning_time_start.tv_sec*1000000 + (binning_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_binning_time += (tE - tS);
        #define COUNT_INDEXING_TIME_END gettimeofday(&current_time, NULL); \
                                      tS = indexing_time_start.tv_sec*1000000 + (indexing_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_indexing_time += (tE - tS);
        #define COUNT_QUERYING_TIME_END gettimeofday(&current_time, NULL); \
                                      tS = querying_time_start.tv_sec*1000000 + (querying_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_querying_time += (tE - tS);
        #define COUNT_INPUT_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = input_time_start.tv_sec*1000000 + (input_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_input_time += (tE - tS); 
        #define COUNT_INPUT2_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = input2_time_start.tv_sec*1000000 + (input2_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_input2_time += (tE - tS); 
        #define COUNT_CAL_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = cal_time_start.tv_sec*1000000 + (cal_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_cal_time += (tE - tS); 
        #define COUNT_OUTPUT_TIME_END gettimeofday(&current_time, NULL); \
                                      tS = output_time_start.tv_sec*1000000 + (output_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_output_time += (tE - tS); 
        #define COUNT_OUTPUT2_TIME_END gettimeofday(&current_time, NULL); \
                                      tS = output2_time_start.tv_sec*1000000 + (output2_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_output2_time += (tE - tS); 
        #define COUNT_TOTAL_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = total_time_start.tv_sec*1000000 + (total_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_time += (tE - tS); 
        #define COUNT_PARSE_TIME_END  gettimeofday(&current_time, NULL);  \
                                      tS = parse_time_start.tv_sec*1000000 + (parse_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_parse_time += (tE - tS); 
        #define COUNT_LOAD_NODE_END   gettimeofday(&current_time, NULL);  \
                                      tS = load_node_start.tv_sec*1000000 + (load_node_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      load_node_time += (tE - tS); 
        #define COUNT_QUEUE_TIME_END  gettimeofday(&current_time, NULL);  \
                                      tS = queue_time_start.tv_sec*1000000 + (queue_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_queue_time += (tE - tS); 
    #else
        #define INIT_STATS() ;
        #define PRINT_STATS() ;
        #define COUNT_NEW_NODE() ;
        #define COUNT_CHECKED_NODE();
        #define COUNT_LOADED_NODE() ;
        #define COUNT_LOADED_RECORD() ;
        #define COUNT_INPUT_TIME_START ;
        #define COUNT_INPUT_TIME_END ;
        #define COUNT_OUTPUT_TIME_START ;
        #define COUNT_OUTPUT_TIME_END ;
        #define COUNT_TOTAL_TIME_START ;
        #define COUNT_TOTAL_TIME_END ;
        #define COUNT_PARSE_TIME_START ;
        #define COUNT_PARSE_TIME_END ;
    #endif
#endif
