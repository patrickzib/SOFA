#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ads/calc_utils.h"
#include "ads/spartan/spartan.h"

void sfa_from_fft(isax_index *index, const ts_type *cur_transform, sax_type *cur_sfa_word);

static unsigned int state = 1;

static unsigned int next_random(void) {
    state = state * 1664525U + 1013904223U;
    return state;
}

static sax_type linear_spartan_symbol(const ts_type *boundaries, int count, ts_type value) {
    int symbol = 0;
    while (symbol < count && value >= boundaries[symbol]) ++symbol;
    return (sax_type) symbol;
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
    for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
    free(index.bins);
    return 0;
}
