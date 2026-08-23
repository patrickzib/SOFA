#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ads/calc_utils.h"

static unsigned int state = 1;

static unsigned int next_random(void) {
    state = state * 1664525U + 1013904223U;
    return state;
}

int main(void) {
    enum { DIMENSIONS = 32, RECORD_DIMENSIONS = 16, BINS = 255 };
    isax_index index;
    isax_index_settings settings;
    float values[DIMENSIONS];
    sax_type sax_min[DIMENSIONS], sax_max[DIMENSIONS], cardinalities[DIMENSIONS];
    const uint64_t record_dimensions = (UINT64_C(1) << RECORD_DIMENSIONS) - 1;

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
            float prefix = 0.0f;
            settings.function_type = function_type;
            const float full = messi_minidist_range_raw(&index, values, sax_min, sax_max,
                                                         cardinalities, FLT_MAX);
            const float partitioned = messi_minidist_range_raw_partitioned(
                &index, values, sax_min, sax_max, cardinalities, FLT_MAX,
                record_dimensions, &suffix);
            (void) messi_minidist_range_raw_partitioned(
                &index, values, sax_min, sax_max, cardinalities, FLT_MAX,
                ~record_dimensions, &prefix);
            if (fabsf(full - partitioned) > 1e-4f * fmaxf(1.0f, full) ||
                fabsf(full - (prefix + suffix)) > 1e-4f * fmaxf(1.0f, full)) {
                fprintf(stderr, "partitioned MBR mismatch: type=%d full=%g prefix=%g suffix=%g\n",
                        function_type, full, prefix, suffix);
                for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
                free(index.bins);
                return 1;
            }
        }
    }
    for (int i = 0; i < DIMENSIONS; ++i) free(index.bins[i]);
    free(index.bins);
    return 0;
}
