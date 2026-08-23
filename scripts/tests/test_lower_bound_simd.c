#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ads/lower_bound_simd.h"

static unsigned int state = 1;

static unsigned int next_random(void) {
    state = state * 1664525U + 1013904223U;
    return state;
}

int main(void) {
    isax_index index;
    isax_index_settings settings;
    memset(&index, 0, sizeof(index));
    memset(&settings, 0, sizeof(settings));
    settings.sax_bit_cardinality = 8;
    settings.sax_alphabet_cardinality = 256;
    index.settings = &settings;
    index.binsv = malloc(sizeof(*index.binsv) * 16 * 255);
    if (index.binsv == NULL) return 1;
    for (int dimension = 0; dimension < 16; ++dimension)
        for (int symbol = 0; symbol < 255; ++symbol)
            index.binsv[dimension * 255 + symbol] = ((float) symbol - 127.0f) / 16.0f;

    for (int trial = 0; trial < 1000; ++trial) {
        float values[16];
        sax_type sax[16], cardinalities[16];
        for (int i = 0; i < 16; ++i) {
            values[i] = ((float) (next_random() & 0xffff) / 4096.0f) - 8.0f;
            sax[i] = (sax_type) (next_random() >> 24);
            cardinalities[i] = (sax_type) (4 + (next_random() % 5));
        }
        for (int factor_index = 0; factor_index < 2; ++factor_index) {
            const float factor = factor_index == 0 ? 1.0f : 2.0f;
            const float reference = messi_lower_bound_16_scalar(
                &index, values, sax, cardinalities, FLT_MAX, factor);
            const float actual = messi_lower_bound_16(
                &index, values, sax, cardinalities, FLT_MAX, factor);
            if (fabsf(reference - actual) > 1e-4f * fmaxf(1.0f, reference)) {
                fprintf(stderr, "lower-bound mismatch: reference=%g actual=%g\n", reference, actual);
                free(index.binsv);
                return 1;
            }
        }
    }
    free(index.binsv);
    return 0;
}
