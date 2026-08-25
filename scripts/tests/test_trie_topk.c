#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ads/api.h"

static int compare_float(const void *left, const void *right) {
    const float a = *(const float *) left, b = *(const float *) right;
    return a < b ? -1 : a > b;
}

int main(void) {
    enum { RECORDS = 12, LENGTH = 16, K = 4 };
    float data[RECORDS][LENGTH];
    float query[LENGTH] = {0};
    float expected[RECORDS];
    char path[] = "/tmp/messi-trie-topk-XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) return 1;
    FILE *file = fdopen(fd, "wb");
    if (file == NULL) { close(fd); unlink(path); return 1; }
    for (int record = 0; record < RECORDS; ++record) {
        for (int d = 0; d < LENGTH; ++d)
            /* SAX bounds assume normalized-scale values; keep the fixture in
             * that range while retaining distinct raw Euclidean neighbours. */
            data[record][d] = ((float) record - 6.0f) * 0.05f +
                              (float) (d % 3) * 0.02f;
        expected[record] = 0.0f;
        for (int d = 0; d < LENGTH; ++d) {
            float delta = data[record][d] - query[d];
            expected[record] += delta * delta;
        }
    }
    /* Deliberate equal-distance records verify deterministic ID ordering. */
    memcpy(data[5], data[6], sizeof(data[5]));
    expected[5] = expected[6];
    if (fwrite(data, sizeof(data[0]), RECORDS, file) != RECORDS || fclose(file) != 0) {
        unlink(path); return 1;
    }

    messi_index_params params;
    memset(&params, 0, sizeof(params));
    params.root_directory = "/tmp";
    params.timeseries_size = LENGTH;
    params.n_segments = LENGTH;
    params.sax_bit_cardinality = 8;
    params.max_leaf_size = 3;
    params.min_leaf_size = 1;
    params.initial_leaf_buffer_size = 16;
    params.initial_fbl_buffer_size = 16;
    params.max_total_buffer_size = 1024;
    params.function_type = 3;
    params.sample_size = RECORDS;
    params.is_norm = 1;
    params.histogram_type = 1;
    params.sample_type = 1;
    params.max_query_threads = 1;
    params.queue_count = 1;
    params.index_type = MESSI_INDEX_TRIE;
    params.sampling_seed = 1;
    params.node_split_criterion = 1;
    params.trie_bound_dimensions = LENGTH;
    params.trie_split_dimensions = LENGTH;
    params.trie_fanout = 2;

    messi_index *index = messi_index_create(&params);
    float distances[K], parallel_distances[K];
    long labels[K], parallel_labels[K];
    qsort(expected, RECORDS, sizeof(expected[0]), compare_float);
    int ok = index != NULL && messi_index_add_file(index, path, RECORDS, 1) == 0 &&
             messi_index_search(index, query, 1, LENGTH, K, distances, labels, 1) == 0;
    if (ok) for (int rank = 0; rank < K; ++rank) {
        if (fabsf(distances[rank] - expected[rank]) > 1e-4f ||
            (rank > 0 && (distances[rank] < distances[rank - 1] ||
                          (distances[rank] == distances[rank - 1] && labels[rank] < labels[rank - 1])))) {
            fprintf(stderr, "top-k mismatch at rank %d: got (%g, %ld), expected distance %g\n",
                    rank, distances[rank], labels[rank], expected[rank]);
            ok = 0;
            break;
        }
    }
    if (index != NULL) messi_index_destroy(index);

    /* The leaf-worker implementation must return exactly the serial order. */
    params.max_query_threads = 4;
    params.queue_count = 4;
    index = messi_index_create(&params);
    if (index == NULL || messi_index_add_file(index, path, RECORDS, 1) != 0 ||
        messi_index_search(index, query, 1, LENGTH, K, parallel_distances, parallel_labels, 1) != 0) {
        ok = 0;
    } else for (int rank = 0; rank < K; ++rank) {
        if (parallel_distances[rank] != distances[rank] || parallel_labels[rank] != labels[rank]) {
            fprintf(stderr, "parallel top-k mismatch at rank %d: serial=(%g,%ld) parallel=(%g,%ld)\n",
                    rank, distances[rank], labels[rank], parallel_distances[rank], parallel_labels[rank]);
            ok = 0;
            break;
        }
    }
    if (index != NULL) messi_index_destroy(index);
    unlink(path);
    return ok ? 0 : 1;
}
