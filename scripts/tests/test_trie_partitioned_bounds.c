#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ads/api.h"
#include "ads/calc_utils.h"
#include "ads/lower_bound_simd.h"
#include "ads/spartan/pca.h"
#include "ads/spartan/spartan.h"

void sfa_from_fft(isax_index *index, const ts_type *cur_transform, sax_type *cur_sfa_word);

static unsigned int state = 1;

extern int maxquerythread;

static unsigned int next_random(void) {
    state = state * 1664525U + 1013904223U;
    return state;
}

static sax_type linear_spartan_symbol(const ts_type *boundaries, int count, ts_type value) {
    int symbol = 0;
    while (symbol < count && value >= boundaries[symbol]) ++symbol;
    return (sax_type) symbol;
}

static int test_parallel_leaf_list_ivf_ranges(void) {
    enum { RECORDS = 5000, LENGTH = 16, QUERY_RECORD = 4321 };
    char path[] = "/tmp/messi-leaf-list-XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) return 0;
    FILE *file = fdopen(fd, "wb");
    float *dataset = malloc((size_t) RECORDS * LENGTH * sizeof(*dataset));
    if (file == NULL || dataset == NULL) {
        if (file != NULL) fclose(file); else close(fd);
        free(dataset); unlink(path); return 0;
    }
    for (int record = 0; record < RECORDS; ++record)
        for (int dimension = 0; dimension < LENGTH; ++dimension)
            dataset[(size_t) record * LENGTH + dimension] =
                (float) record * 0.001f + (float) dimension * 0.01f +
                sinf((float) (record * 17 + dimension * 13)) * 0.0001f;
    if (fwrite(dataset, sizeof(*dataset), (size_t) RECORDS * LENGTH, file) !=
            (size_t) RECORDS * LENGTH || fclose(file) != 0) {
        free(dataset); unlink(path); return 0;
    }

    messi_index_params params;
    memset(&params, 0, sizeof(params));
    params.root_directory = "/tmp";
    params.timeseries_size = LENGTH;
    params.n_segments = 8;
    params.sax_bit_cardinality = 8;
    params.max_leaf_size = 10000;
    params.min_leaf_size = 1;
    params.initial_leaf_buffer_size = 10000;
    params.max_total_buffer_size = 1000000;
    params.initial_fbl_buffer_size = 100;
    params.total_loaded_leaves = 1;
    params.tight_bound = 1;
    params.function_type = 4;
    params.sample_size = 1000;
    params.is_norm = 1;
    params.histogram_type = 2;
    params.sample_type = 1;
    params.n_coefficients = 8;
    params.max_query_threads = 8;
    params.queue_count = 8;
    params.index_type = MESSI_INDEX_TRIE;
    params.sampling_seed = 1;
    params.node_split_criterion = 1;
    params.trie_bound_dimensions = 8;
    params.trie_split_dimensions = 8;
    params.trie_record_mbr_suffix_bound = 1;
    params.trie_query_engine = MESSI_TRIE_QUERY_LEAF_LIST;
    params.trie_leaf_ivf = 2;
    params.trie_fanout = 8;

    messi_index *trie = messi_index_create(&params);
    if (trie == NULL || messi_index_add_file(trie, path, RECORDS, 0) != 0) {
        messi_index_destroy(trie); free(dataset); unlink(path); return 0;
    }
    float serial_distance, parallel_distance;
    long serial_label, parallel_label;
    maxquerythread = 1;
    int serial_status = messi_index_search(trie, dataset + (size_t) QUERY_RECORD * LENGTH,
                                            1, LENGTH, 1, &serial_distance, &serial_label, 0);
    maxquerythread = 8;
    int parallel_status = messi_index_search(trie, dataset + (size_t) QUERY_RECORD * LENGTH,
                                              1, LENGTH, 1, &parallel_distance, &parallel_label, 0);
    messi_index_destroy(trie);
    free(dataset);
    unlink(path);
    if (serial_status != 0 || parallel_status != 0 || serial_label != QUERY_RECORD ||
        parallel_label != QUERY_RECORD || serial_distance != 0.0f || parallel_distance != 0.0f) {
        fprintf(stderr,
                "leaf-list IVF range mismatch: serial=(%g,%ld,%d) parallel=(%g,%ld,%d)\n",
                serial_distance, serial_label, serial_status,
                parallel_distance, parallel_label, parallel_status);
        return 0;
    }
    return 1;
}

int main(void) {
    enum { DIMENSIONS = 128, RECORD_DIMENSIONS = 16, BINS = 255 };
    isax_index index;
    isax_index_settings settings;
    float values[DIMENSIONS];
    sax_type sax_min[DIMENSIONS], sax_max[DIMENSIONS], cardinalities[DIMENSIONS];

    memset(&index, 0, sizeof(index));
    memset(&settings, 0, sizeof(settings));
    settings.n_segments = DIMENSIONS;
    settings.sax_bit_cardinality = 8;
    settings.sax_alphabet_cardinality = 256;
    settings.is_norm = 1;
    index.settings = &settings;
    index.bins = calloc(DIMENSIONS, sizeof(*index.bins));
    if (index.bins == NULL) return 1;
    for (int dimension = 0; dimension < DIMENSIONS; ++dimension) {
        index.bins[dimension] = malloc(BINS * sizeof(*index.bins[dimension]));
        if (index.bins[dimension] == NULL) {
            for (int i = 0; i < dimension; ++i) free(index.bins[i]);
            free(index.bins);
            return 1;
        }
        for (int symbol = 0; symbol < BINS; ++symbol)
            index.bins[dimension][symbol] = ((float) symbol - 127.0f) / 16.0f;
    }

    /* Binary SPARTAN quantization must preserve the former strict-boundary
     * linear scan for values below, at, between, and above bin boundaries. */
    for (int trial = 0; trial < 1000; ++trial) {
        sax_type actual[DIMENSIONS];
        for (int dimension = 0; dimension < DIMENSIONS; ++dimension) {
            const int symbol = (int) (next_random() % 256U);
            const float boundary = symbol == 255 ? 8.0f : index.bins[dimension][symbol];
            const int mode = trial % 4;
            values[dimension] = mode == 0 ? boundary - 0.001f
                              : mode == 1 ? boundary
                              : mode == 2 ? boundary + 0.001f
                                          : ((float) (next_random() & 0xffff) / 4096.0f) - 8.0f;
        }
        spartan_from_pca(&index, values, actual);
        for (int dimension = 0; dimension < DIMENSIONS; ++dimension) {
            const sax_type expected = linear_spartan_symbol(index.bins[dimension], BINS, values[dimension]);
            if (actual[dimension] != expected) {
                fprintf(stderr, "SPARTAN symbol mismatch: dimension=%d value=%g expected=%u actual=%u\n",
                        dimension, values[dimension], (unsigned int) expected, (unsigned int) actual[dimension]);
                for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
                free(index.bins);
                return 1;
            }
        }

        /* PISA delegates its final symbolization to this same SFA function. */
        sfa_from_fft(&index, values, actual);
        for (int dimension = 0; dimension < DIMENSIONS; ++dimension) {
            const sax_type expected = linear_spartan_symbol(index.bins[dimension], BINS, values[dimension]);
            if (actual[dimension] != expected) {
                fprintf(stderr, "SFA/PISA symbol mismatch: dimension=%d value=%g expected=%u actual=%u\n",
                        dimension, values[dimension], (unsigned int) expected, (unsigned int) actual[dimension]);
                for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
                free(index.bins);
                return 1;
            }
        }
    }

    /* PISA projects a distance-preserving packed real FFT through an
     * orthonormal PCA basis.  Its coordinates are therefore ordinary
     * Euclidean dimensions: unlike SFA's unpaired FFT coordinates, they must
     * not receive a further factor of two or any DC-coordinate exception. */
    settings.function_type = 6;
    settings.n_coefficients = DIMENSIONS;
    for (int i = 0; i < DIMENSIONS; ++i) {
        values[i] = 10.0f;
        sax_min[i] = sax_max[i] = 127;
        cardinalities[i] = 8;
    }
    const float one_dimension = values[0] - index.bins[0][127];
    const float expected_pisa_full = DIMENSIONS * one_dimension * one_dimension;
    const float pisa_range = messi_minidist_range_raw(&index, values, sax_min, sax_max,
                                                       cardinalities, FLT_MAX);
    const float pisa_record = messi_minidist_raw(&index, values, sax_min, cardinalities, FLT_MAX);
    float pisa_table[MESSI_RECORD_LB_MAX_DIMENSIONS][256];
    if (!messi_build_record_lb_table(&index, values, MESSI_RECORD_LB_MAX_DIMENSIONS, pisa_table)) {
        fprintf(stderr, "failed to build PISA record lower-bound table\n");
        for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
        free(index.bins);
        return 1;
    }
    const float pisa_table_sum = messi_record_lb_table_sum(
        pisa_table, sax_min, MESSI_RECORD_LB_MAX_DIMENSIONS, FLT_MAX, 0);
    const float expected_pisa_record = MESSI_RECORD_LB_MAX_DIMENSIONS * one_dimension * one_dimension;
    if (fabsf(pisa_range - expected_pisa_full) > 1e-4f * expected_pisa_full ||
        fabsf(pisa_record - expected_pisa_full) > 1e-4f * expected_pisa_full ||
        fabsf(pisa_table_sum - expected_pisa_record) > 1e-4f * expected_pisa_record) {
        fprintf(stderr, "PISA lower-bound scale mismatch: range=%g record=%g table=%g\n",
                pisa_range, pisa_record, pisa_table_sum);
        for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
        free(index.bins);
        return 1;
    }

    for (int trial = 0; trial < 1000; ++trial) {
        for (int i = 0; i < DIMENSIONS; ++i) {
            const sax_type first = (sax_type) (next_random() >> 24);
            const sax_type second = (sax_type) (next_random() >> 24);
            values[i] = ((float) (next_random() & 0xffff) / 4096.0f) - 8.0f;
            sax_min[i] = first < second ? first : second;
            sax_max[i] = first < second ? second : first;
            cardinalities[i] = 8;
        }
        for (int function_type = 4; function_type <= 6; ++function_type) {
            float suffix = 0.0f;
            settings.function_type = function_type;
            const float full = messi_minidist_range_raw(&index, values, sax_min, sax_max,
                                                         cardinalities, FLT_MAX);
            const float partitioned = messi_minidist_range_raw_partitioned(
                &index, values, sax_min, sax_max, cardinalities, FLT_MAX,
                RECORD_DIMENSIONS, &suffix);
            if (fabsf(full - partitioned) > 1e-4f * fmaxf(1.0f, full) ||
                suffix < 0.0f || suffix > full + 1e-4f * fmaxf(1.0f, full)) {
                fprintf(stderr, "partitioned MBR mismatch: type=%d full=%g suffix=%g\n",
                        function_type, full, suffix);
                for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
                free(index.bins);
                return 1;
            }
        }
    }

    /* Ensure the suffix really includes dimensions above the old 64-bit mask
     * limit: make only the suffix lie outside its MBR. */
    for (int i = 0; i < DIMENSIONS; ++i) {
        values[i] = i < RECORD_DIMENSIONS ? 0.0f : 100.0f;
        sax_min[i] = 120;
        sax_max[i] = 136;
        cardinalities[i] = 8;
    }
    for (int function_type = 4; function_type <= 6; ++function_type) {
        float suffix = 0.0f;
        settings.function_type = function_type;
        const float full = messi_minidist_range_raw_partitioned(
            &index, values, sax_min, sax_max, cardinalities, FLT_MAX,
            RECORD_DIMENSIONS, &suffix);
        if (fabsf(full - suffix) > 1e-4f * fmaxf(1.0f, full)) {
            fprintf(stderr, "128-D suffix mismatch: type=%d full=%g suffix=%g\n",
                    function_type, full, suffix);
            for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
            free(index.bins);
            return 1;
        }
    }

    /* The optional CBLAS batch projection must agree with the scalar PCA
     * projection used for queries and for builds without CBLAS. */
    enum { BATCH_ROWS = 3 };
    ts_type batch_input[BATCH_ROWS * DIMENSIONS];
    ts_type batch_output[BATCH_ROWS * DIMENSIONS];
    ts_type scalar_output[DIMENSIONS];
    index.pca_dim = DIMENSIONS;
    index.pca_components_count = DIMENSIONS;
    index.pca_components = calloc((size_t) DIMENSIONS * DIMENSIONS,
                                  sizeof(*index.pca_components));
    index.pca_bias = calloc(DIMENSIONS, sizeof(*index.pca_bias));
    if (index.pca_components == NULL || index.pca_bias == NULL) {
        pca_free(&index);
        for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
        free(index.bins);
        return 1;
    }
    for (int i = 0; i < DIMENSIONS; ++i) {
        index.pca_components[i * DIMENSIONS + i] = 1.0f;
        index.pca_bias[i] = (double) (i % 7) * 0.125;
    }
    for (int row = 0; row < BATCH_ROWS; ++row)
        for (int i = 0; i < DIMENSIONS; ++i)
            batch_input[row * DIMENSIONS + i] = (ts_type) (row * 100 + i) / 16.0f;
    if (pca_project_batch(&index, batch_input, BATCH_ROWS, batch_output, BATCH_ROWS) != SUCCESS) {
        fprintf(stderr, "PCA batch projection failed\n");
        pca_free(&index);
        for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
        free(index.bins);
        return 1;
    }
    for (int row = 0; row < BATCH_ROWS; ++row) {
        if (pca_from_ts(&index, batch_input + row * DIMENSIONS, scalar_output) != SUCCESS) {
            fprintf(stderr, "PCA scalar projection failed\n");
            pca_free(&index);
            for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
            free(index.bins);
            return 1;
        }
        for (int i = 0; i < DIMENSIONS; ++i) {
            if (fabsf(batch_output[row * DIMENSIONS + i] - scalar_output[i]) > 1e-5f) {
                fprintf(stderr, "PCA batch mismatch: row=%d dimension=%d batch=%g scalar=%g\n",
                        row, i, batch_output[row * DIMENSIONS + i], scalar_output[i]);
                pca_free(&index);
                for (int j = 0; j < DIMENSIONS; ++j) free(index.bins[j]);
                free(index.bins);
                return 1;
            }
        }
    }
    pca_free(&index);
    for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
    free(index.bins);
    if (!test_parallel_leaf_list_ivf_ranges()) return 1;
    return 0;
}
