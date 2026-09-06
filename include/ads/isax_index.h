//
//  isax_index.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//

#ifndef isaxlib_isax_index_h
#define isaxlib_isax_index_h
#include <pthread.h>
#include <stdbool.h>
#include <time.h>
#include <sys/stat.h>
#include <stdio.h>

#include "config.h"
#include "../../globals.h"
#include "isax_node.h"
#include "isax_node_record.h"
#include "sax/ts.h"
#include "pqueue.h"
typedef struct {
    unsigned long mem_tree_structure;
    unsigned long mem_data; 
    unsigned long mem_summaries;
    unsigned long disk_data_full;
    unsigned long disk_data_partial;
} meminfo;

typedef enum {
    MESSI_INDEX_ISAX = 0,
    MESSI_INDEX_TRIE = 1
} messi_index_type;


typedef struct {
    char new_index;

    char * raw_filename;
    const char* root_directory;
    int initial_fbl_buffer_size;
	sax_type *max_sax_cardinalities;
    
    // ALWAYS: TIMESERIES_SIZE = TS_VALUES_PER_PAA_SEGMENT * PAA_SEGMENTS
    int timeseries_size;
    int ts_values_per_paa_segment;
    int n_segments;
	
	int tight_bound;
	int aggressive_check;
    
    int sax_byte_size;
    int position_byte_size;
    int ts_byte_size;
    
    int full_record_size;
    int partial_record_size;
    
    // ALWAYS: SAX_ALPHABET_CARDINALITY = 2^SAX_BIT_CARDINALITY
    int sax_bit_cardinality;
    root_mask_type * bit_masks;
    int sax_alphabet_cardinality;
    
    int max_leaf_size;
    int min_leaf_size;
    int initial_leaf_buffer_size;
    int max_total_buffer_size;
    int max_total_full_buffer_size;
    
    int max_filename_size;
    
    float mindist_sqrt;
    int root_nodes_size;
    
    int total_loaded_leaves;

    char SIMD_flag;
    char is_norm;
    
    int function_type;
    unsigned int sample_size;
    int histogram_type;
    int sample_type;
    unsigned int sampling_seed;
    int n_coefficients;
    int node_split_criterion;
    messi_index_type index_type;
    /* Optional iSAX SOFA v2 bounds.  The primary iSAX word/tree always
     * remains n_segments wide; MBR suffix dimensions live in a sidecar. */
    char isax_node_mbr;
    char isax_record_mbr_suffix_bound;
    char isax_record_lb_table;
    int isax_mbr_dimensions;
    /* Trie-only: public --n-segments remains the lower-bound dimensionality,
     * while the trie may materialize a wider symbolic word for splitting. */
    int trie_bound_dimensions;
    /* Trie-only: prefix of the symbolic word considered while selecting a
     * split.  MBRs always retain the full symbolic word. */
    int trie_split_dimensions;
    char trie_record_mbr_suffix_bound;
    /* Refine passing records immediately instead of ordering them in a
     * per-leaf lower-bound heap.  This is the trie default; clearing it
     * restores the heap for A/B benchmarking. */
    char trie_streaming_leaf_scan;
    /* Optional flat IVF/MRB directory inside large terminal trie leaves. */
    int trie_leaf_ivf;
    /* Certified raw-space centroid/radius bound for trie IVF clusters. */
    char trie_leaf_ivf_raw_ball_bound;
    /* Opt-in per-record centroid-radius triangle bound. */
    char trie_leaf_ivf_radial_bound;
    /* Calibrate the radial bound per query and keep it only when at least
     * one quarter of sampled candidates are rejected. */
    char trie_leaf_ivf_radial_bound_auto;
    /* Trie-only symbolic partition fanout.  Fixed mode accepts 2, 4, or 8;
     * dynamic mode derives per-dimension fanouts up to the 8-bit alphabet. */
    int trie_fanout;
    char trie_dynamic_alphabet;
    int trie_min_bits;
    int trie_max_bits;
    int trie_alphabet_budget_bits;
    double *symbolic_variances;

    /* Root-only partitioning for in-memory symbolic indexes.  When enabled,
     * this contains the number of leading SAX bits used per dimension. */
    char dynamic_root_split_variance;
    sax_type *root_bit_cardinalities;

    // int filetype_int;

} isax_index_settings;

typedef struct {
    meminfo memory_info;

    FILE *sax_file;
    sax_type *sax_cache;
    unsigned long sax_cache_size;
    /* One byte per enabled suffix dimension and record.  This is allocated
     * only for SOFA v2 suffix MBRs, never for the default iSAX path. */
    sax_type *mbr_suffix_sax_cache;
    unsigned long mbr_suffix_sax_cache_size;

    unsigned long long allocated_memory;
    unsigned long root_nodes;
    unsigned long long total_records;
    unsigned long long loaded_records;
    
    int * locations;
    struct isax_node *first_node;
    isax_index_settings *settings;
    
    char has_wedges;

    struct first_buffer_layer *fbl;

    ts_type *answer;

    ts_type **bins;
    ts_type *binsv;
    ts_type norm_factor;

    int * coefficients;
    ts_type *pca_mean;
    ts_type *pca_components;
    double *pca_bias;
    double *pca_explained_variance;
    int pca_components_count;
    int pca_dim;

    struct symbolic_trie_index *trie;

} isax_index;



//TODO: Put sanity check for variables (cardinalities etc.)

isax_index * isax_index_init(isax_index_settings *settings);
isax_index_settings * isax_index_settings_init (const char * root_directory,
                                                int timeseries_size, 
                                                int n_segments, 
                                                int sax_bit_cardinality,
                                                int max_leaf_size,
                                                int min_leaf_size,
                                                int initial_leaf_buffer_size,
                                                int max_total_buffer_size, 
                                                int initial_fbl_buffer_size,
                                                int total_loaded_leaves,
												int tight_bound, int aggressive_check, int new_index,
                                                int function_type, char inmemory_flag, char SIMD_flag,
                                                int sample_size, char is_norm, int histogram_type,
                                                int sample_type, int n_coefficients,
                                                messi_index_type index_type);
void print_settings(isax_index_settings *settings, int query_workers, int trie_query_batch);

isax_node * add_record_to_node(isax_index *index, isax_node *node,
                               isax_node_record *record,
                               const char leaf_size_check, int kn);
root_mask_type isax_fbl_index_insert(isax_index *index, 
                                    sax_type * sax,
                                    file_position_type * pos);
enum response isax_index_insert(isax_index *index, sax_type *sax, 
                                file_position_type *pos);
enum response flush_subtree_leaf_buffers (isax_index *index, isax_node *node);
enum response flush_subtree_leaf_buffers_m (isax_index *index, isax_node *node, pthread_mutex_t *lock_index, pthread_mutex_t *lock_write);
enum response flush_subtree_leaf_buffers_m1(isax_index *index, isax_node *node, pthread_mutex_t *lock_index, pthread_mutex_t *lock_write);
enum response flush_all_leaf_buffers(isax_index *index, enum buffer_cleaning_mode buffer_clean_mode);
enum response create_node_filename(isax_index *index,
                                   isax_node *node,
                                   isax_node_record *record,
                                   int kn);
void isax_index_clear_node_buffers(isax_index *index, isax_node *node, 
                                   enum node_cleaning_mode node_cleaning_mode,
                                   enum buffer_cleaning_mode buffer_clean_mode);
void isax_index_destroy(isax_index *index, isax_node *node);
void MESSI2_index_destroy(isax_index *index, isax_node *node);

void isax_index_pRecBuf_destroy(isax_index *index, isax_node *node,int prewokernumber);
int comp(const void * a, const void * b);
void load_random_leaf(isax_index *index);

float calculate_node_distance (isax_index *index, isax_node *node, ts_type *query, float bsf);
void calculate_node_topk (isax_index *index, isax_node *node, ts_type *query, pqueue_bsf *pq_bsf);
float calculate_node_distance_SIMD (isax_index *index, isax_node *node, ts_type *query, float bsf);
float isax_index_load_node(isax_index *index, isax_node *c_node, ts_type * query, float bsf);
void isax_index_load_node_topk(isax_index *index, isax_node *c_node, ts_type * query, pqueue_bsf *pq_bsf);

void complete_index(isax_index *index, int ts_num);
void complete_index_leafs(isax_index *index);
void complete_subtree_leafs(isax_index *index, isax_node *node);

int find_siblings(isax_node *c_node, isax_node **nodes_to_load, 
                  int *number_of_leaves, int *offset);
float calculate_minimum_distance (isax_index *index, isax_node *node, ts_type *raw_query, ts_type *query);
float calculate_minimum_distance_SIMD (isax_index *index, isax_node *node, ts_type *raw_query, ts_type *query);
void cache_sax_file(isax_index *index);



void index_write(isax_index *index);
void index_mRecBuf_write(isax_index *index);
void node_write(isax_index *index, isax_node *node, FILE *file);
isax_index * index_read(const char * root_directory);
void get_index_size(isax_index *index, struct stat *stat_index, struct stat *stat_adaptive);
isax_node * node_read(isax_index *index, FILE *file);

//// WEDGES FUNCTIONALITY ////
void create_wedges(isax_index *index, isax_node *node);
void clear_wedges(isax_index *index, isax_node *node);
int compare_file_positions (const void * a, const void * b);
void print_mem_info(isax_index *index);
meminfo get_memory_utilization_info(isax_index *index);
#endif
