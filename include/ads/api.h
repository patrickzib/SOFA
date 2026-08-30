#ifndef MESSI_API_H
#define MESSI_API_H

#include "isax_index.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct messi_index messi_index;

typedef struct {
    const char *root_directory;
    int timeseries_size;
    int n_segments;
    int sax_bit_cardinality;
    int max_leaf_size;
    int min_leaf_size;
    int initial_leaf_buffer_size;
    int max_total_buffer_size;
    int initial_fbl_buffer_size;
    int total_loaded_leaves;
    int tight_bound;
    int aggressive_check;
    int function_type;
    char simd;
    int sample_size;
    char is_norm;
    int histogram_type;
    int sample_type;
    int n_coefficients;
    int filetype_int;
    int max_query_threads;
    int queue_count;
    int index_type;
    unsigned int sampling_seed;
    int node_split_criterion;
    char isax_node_mbr;
    char isax_record_mbr_suffix_bound;
    char isax_record_lb_table;
    int isax_mbr_dimensions;
    int trie_bound_dimensions;
    int trie_split_dimensions;
    char trie_record_mbr_suffix_bound;
    int trie_query_engine;
    int trie_leaf_ivf;
    int trie_fanout;
    char trie_dynamic_alphabet;
    int trie_min_fanout;
    int trie_max_fanout;
    int trie_alphabet_budget_bits;
    char dynamic_root_split_variance;
} messi_index_params;

messi_index *messi_index_create(const messi_index_params *params);
void messi_index_destroy(messi_index *index);
int messi_index_add_file(messi_index *index, const char *path, long ts_num, int dynamic_index);
int messi_index_search(messi_index *index, const float *queries, size_t nq, size_t dim, size_t k,
                       float *distances, long *labels, int dynamic_index);
int messi_index_pca_transform(messi_index *index,
                              const float *queries,
                              size_t nq,
                              size_t dim,
                              float *out,
                              size_t out_dim);

#ifdef __cplusplus
}
#endif

#endif /* MESSI_API_H */
