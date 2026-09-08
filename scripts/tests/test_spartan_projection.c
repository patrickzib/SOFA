#include <assert.h>
#include <stdlib.h>
#include <string.h>

static int fail_allocation;
static void *projection_calloc(size_t count, size_t size) {
    return fail_allocation ? NULL : calloc(count, size);
}

/* Exercise the private helper and the actual histogram builder. */
#define calloc projection_calloc
#include "../../src/ads/spartan/spartan.c"
#undef calloc

int main(void) {
    enum { ROWS = 103, INPUT = 17, OUTPUT = 8, SYMBOLS = 4 };
    isax_index_settings settings = {0};
    isax_index index = {0};
    ts_type samples[ROWS * INPUT], components[OUTPUT * INPUT];
    double bias[OUTPUT];
    ts_type expected[OUTPUT][ROWS], actual[OUTPUT][ROWS], scratch[OUTPUT];
    ts_type bins[OUTPUT][SYMBOLS - 1], reference_bins[OUTPUT][SYMBOLS - 1];
    ts_type *expected_columns[OUTPUT], *actual_columns[OUTPUT], *bin_columns[OUTPUT];
    settings.n_segments = OUTPUT;
    settings.timeseries_size = INPUT;
    settings.sax_alphabet_cardinality = SYMBOLS;
    index.settings = &settings;
    index.pca_dim = INPUT;
    index.pca_components_count = OUTPUT - 1; /* Also check padded coefficients. */
    index.pca_components = components;
    index.pca_bias = bias;
    index.bins = bin_columns;
    for (int i = 0; i < ROWS * INPUT; ++i) samples[i] = (float) ((i * 37 % 101) - 50) / 19.0f;
    for (int i = 0; i < OUTPUT * INPUT; ++i) components[i] = (float) ((i * 13 % 41) - 20) / 23.0f;
    for (int k = 0; k < OUTPUT; ++k) {
        bias[k] = (k - 3) / 17.0;
        expected_columns[k] = expected[k];
        actual_columns[k] = actual[k];
        bin_columns[k] = bins[k];
    }
    const int counts[] = {1, 7, ROWS};
    const int workers[] = {1, 2, 4, 8};
    for (unsigned int n = 0; n < sizeof(counts) / sizeof(counts[0]); ++n) {
        settings.sample_size = counts[n];
        for (unsigned int w = 0; w < sizeof(workers) / sizeof(workers[0]); ++w) {
            for (int histogram = 1; histogram <= 2; ++histogram) {
                settings.histogram_type = histogram;
                for (int i = 0; i < counts[n]; ++i) {
                    assert(pca_from_ts(&index, samples + i * INPUT, scratch) == SUCCESS);
                    for (int k = 0; k < OUTPUT; ++k) expected[k][i] = scratch[k];
                }
                assert(spartan_project_samples(&index, samples, counts[n], actual_columns, workers[w]) == SUCCESS);
                for (int k = 0; k < OUTPUT; ++k)
                    assert(memcmp(expected[k], actual[k], counts[n] * sizeof(ts_type)) == 0);
                spartan_bins_data data = {&index, expected_columns, 0, OUTPUT};
                spartan_order_divide_worker(&data);
                memcpy(reference_bins, bins, sizeof(bins));
                data.coeff_mem_array = actual_columns;
                spartan_order_divide_worker(&data);
                assert(memcmp(reference_bins, bins, sizeof(bins)) == 0);
            }
        }
    }
    fail_allocation = 1;
    assert(spartan_project_samples(&index, samples, ROWS, actual_columns, 4) == FAILURE);
    fail_allocation = 0;
    index.pca_components = NULL;
    assert(spartan_project_samples(&index, samples, ROWS, actual_columns, 4) == FAILURE);
    puts("SPARTAN projection: bitwise coefficients/bins and failure paths passed");
    return 0;
}
