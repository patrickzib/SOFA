#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ads/lower_bound_simd.h"
#include "ads/sax/sax.h"

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
    index.binsv = malloc(sizeof(*index.binsv) * 64 * 255);
    if (index.binsv == NULL) return 1;
    for (int dimension = 0; dimension < 64; ++dimension)
        for (int symbol = 0; symbol < 255; ++symbol)
            index.binsv[dimension * 255 + symbol] = ((float) symbol - 127.0f) / 16.0f;

    const int learned_widths[] = {16, 32, 64};
    for (int trial = 0; trial < 1000; ++trial) {
        float values[64];
        sax_type sax[64], cardinalities[64];
        for (int i = 0; i < 64; ++i) {
            values[i] = ((float) (next_random() & 0xffff) / 4096.0f) - 8.0f;
            sax[i] = (sax_type) (next_random() >> 24);
            cardinalities[i] = (sax_type) (4 + (next_random() % 5));
        }
        for (size_t width = 0; width < sizeof(learned_widths) / sizeof(learned_widths[0]); ++width) {
            for (int factor_index = 0; factor_index < 2; ++factor_index) {
                const float factor = factor_index == 0 ? 1.0f : 2.0f;
                const float limit = trial % 2 == 0 ? FLT_MAX : 5.0f;
                const float reference = messi_lower_bound_scalar(
                    &index, values, sax, cardinalities, learned_widths[width], limit, factor);
                const float actual = messi_lower_bound_simd(
                    &index, values, sax, cardinalities, learned_widths[width], limit, factor);
                if ((reference > limit) != (actual > limit) ||
                    (reference <= limit &&
                     fabsf(reference - actual) > 1e-4f * fmaxf(1.0f, reference))) {
                    fprintf(stderr,
                            "lower-bound mismatch: dims=%d reference=%g actual=%g limit=%g\n",
                            learned_widths[width], reference, actual, limit);
                    free(index.binsv);
                    return 1;
                }
            }
        }
    }

    /* The direct SAX kernel must consume the complete symbolic word, not
     * silently stop after the historical first sixteen dimensions. */
    for (size_t width = 0; width < sizeof(learned_widths) / sizeof(learned_widths[0]); ++width) {
        const int dimensions = learned_widths[width];
        float paa[64];
        sax_type sax[64], cardinalities[64];
        settings.n_segments = dimensions;
        settings.mindist_sqrt = 1.0f;
        for (int trial = 0; trial < 1000; ++trial) {
            for (int i = 0; i < dimensions; ++i) {
                paa[i] = ((float) (next_random() & 0xffff) / 8192.0f) - 4.0f;
                sax[i] = (sax_type) (next_random() >> 24);
                cardinalities[i] = (sax_type) (4 + (next_random() % 5));
            }
            const float reference = minidist_paa_to_isax(
                paa, sax, cardinalities, &settings, 1);
            const float actual = minidist_paa_to_isax_raw_SIMD(
                paa, sax, cardinalities, &settings);
            if (fabsf(reference - actual) > 1e-4f * fmaxf(1.0f, reference)) {
                fprintf(stderr, "SAX lower-bound mismatch: dims=%d reference=%g actual=%g\n",
                        dimensions, reference, actual);
                free(index.binsv);
                return 1;
            }
        }
    }

    /* Query-local record tables must agree for every supported trie prefix,
     * both with the scalar path and with target-specific gather SIMD. */
    float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256];
    sax_type word[MESSI_RECORD_LB_MAX_DIMENSIONS];
    const int widths[] = {16, 17, 31, 32, 48, 64};
    for (int trial = 0; trial < 1000; ++trial) {
        for (int dimension = 0; dimension < MESSI_RECORD_LB_MAX_DIMENSIONS; ++dimension) {
            word[dimension] = (sax_type) (next_random() >> 24);
            for (int symbol = 0; symbol < 256; ++symbol)
                table[dimension][symbol] = (float) (next_random() & 0xffU) / 4096.0f;
        }
        for (size_t width = 0; width < sizeof(widths) / sizeof(widths[0]); ++width) {
            const float limit = trial % 2 == 0 ? FLT_MAX : 0.75f;
            const float reference = messi_record_lb_table_scalar(table, word, widths[width], limit);
            const float actual = messi_record_lb_table_sum(table, word, widths[width], limit, 1);
            /* SIMD checks a whole vector before testing the BSF, so it may
             * overshoot a bounded scalar sum. Both must agree on pruning; if
             * neither crosses the limit, their complete sums must agree. */
            if ((reference > limit) != (actual > limit) ||
                (reference <= limit && fabsf(reference - actual) > 1e-5f * fmaxf(1.0f, reference))) {
                fprintf(stderr, "record-table mismatch: dims=%d reference=%g actual=%g\n",
                        widths[width], reference, actual);
                free(index.binsv);
                return 1;
            }
        }
    }
    free(index.binsv);
    return 0;
}
