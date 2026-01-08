cdef extern from "ads/api.h":
    ctypedef struct messi_index
    ctypedef struct messi_index_params:
        const char *root_directory
        int timeseries_size
        int n_segments
        int sax_bit_cardinality
        int max_leaf_size
        int min_leaf_size
        int initial_leaf_buffer_size
        int max_total_buffer_size
        int initial_fbl_buffer_size
        int total_loaded_leaves
        int tight_bound
        int aggressive_check
        int function_type
        char simd
        int sample_size
        char is_norm
        int histogram_type
        int sample_type
        int n_coefficients
        int max_query_threads
        int queue_count

    messi_index *messi_index_create(const messi_index_params *params)
    void messi_index_destroy(messi_index *index)
    int messi_index_add_file(messi_index *index, const char *path, long ts_num)
    int messi_index_search(messi_index *index,
                           const float *queries,
                           size_t nq,
                           size_t dim,
                           size_t k,
                           float *distances,
                           long *labels)
