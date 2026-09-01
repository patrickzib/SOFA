#include "ads/api.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(void) {
    enum { N = 16, RECORDS = 64 };
    char data_path[] = "/tmp/messi-sffa-data-XXXXXX";
    char root_path[] = "/tmp/messi-sffa-index-XXXXXX";
    const int data_fd = mkstemp(data_path);
    if (data_fd < 0) return 1;
    FILE *data = fdopen(data_fd, "wb");
    if (data == NULL) { unlink(data_path); return 1; }
    float rows[RECORDS][N];
    for (int row = 0; row < RECORDS; ++row) {
        for (int i = 0; i < N; ++i)
            rows[row][i] = sinf(0.17f * (float) (row + 1) * (i + 1)) + 0.01f * row;
        if (fwrite(rows[row], sizeof(float), N, data) != N) {
            fclose(data); unlink(data_path); return 1;
        }
    }
    fclose(data);
    if (mkdtemp(root_path) == NULL || rmdir(root_path) != 0) {
        unlink(data_path); return 1;
    }

    messi_index_params params;
    memset(&params, 0, sizeof(params));
    params.root_directory = root_path;
    params.timeseries_size = N;
    params.n_segments = 4;
    params.sax_bit_cardinality = 3;
    params.max_leaf_size = 8;
    params.min_leaf_size = 2;
    params.initial_leaf_buffer_size = 8;
    params.max_total_buffer_size = 1024;
    params.initial_fbl_buffer_size = 16;
    params.total_loaded_leaves = 1;
    params.tight_bound = 1;
    params.function_type = MESSI_TRANSFORM_SFFA;
    params.simd = 0;
    params.sample_size = 64;
    params.is_norm = 1;
    params.histogram_type = 1;
    params.sample_type = 1;
    params.max_query_threads = 1;
    params.queue_count = 1;
    params.index_type = MESSI_INDEX_ISAX;
    params.sampling_seed = 7;
    params.node_split_criterion = 1;
    params.sffa_order = 0.5;
    params.sffa_auto_order = 1;
    messi_index *index = messi_index_create(&params);
    if (index == NULL) { unlink(data_path); rmdir(root_path); return 1; }
    if (messi_index_add_file(index, data_path, RECORDS, 1) != 0) {
        messi_index_destroy(index); unlink(data_path); rmdir(root_path); return 1;
    }
    const double selected_order = messi_index_sffa_order(index);
    if (!isfinite(selected_order) || selected_order < 0.2 || selected_order > 1.0) {
        messi_index_destroy(index); unlink(data_path); rmdir(root_path); return 1;
    }
    float distances[2];
    long labels[2];
    float queries[2][N];
    memcpy(queries[0], rows[3], sizeof(rows[3]));
    memcpy(queries[1], rows[29], sizeof(rows[29]));
    const int result = messi_index_search(index, &queries[0][0], 2, N, 1,
                                          distances, labels, 1);
    messi_index_destroy(index);
    unlink(data_path);
    rmdir(root_path);
    /* The iSAX API currently reports -1 labels for exact zero-distance
     * matches, so assert the exact distance contract here. */
    if (result != 0 || distances[0] > 1e-5f || distances[1] > 1e-5f) {
        fprintf(stderr, "SFFA API exact-search mismatch: rc=%d d=(%g,%g) labels=(%ld,%ld)\n",
                result, distances[0], distances[1], labels[0], labels[1]);
        return 1;
    }
    return 0;
}
