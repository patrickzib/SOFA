#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "ads/sax/ts.h"

static unsigned int state = 1;

static unsigned int next_random(void) {
    state = state * 1664525U + 1013904223U;
    return state;
}

static float reference_ed(const float *left, const float *right, int length) {
    float sum = 0.0f;
    for (int i = 0; i < length; ++i) {
        const float difference = left[i] - right[i];
        sum += difference * difference;
    }
    return sum;
}

static int nearly_equal(float expected, float actual) {
    return fabsf(expected - actual) <= 2e-5f * fmaxf(1.0f, expected);
}

int main(void) {
    enum { MAX_LENGTH = 129, TRIALS = 100 };
    float *left_raw = malloc((MAX_LENGTH + 17) * sizeof(*left_raw));
    float *right_raw = malloc((MAX_LENGTH + 17) * sizeof(*right_raw));
    if (left_raw == NULL || right_raw == NULL) return 1;
    float *left = (float *) (((uintptr_t) left_raw + 63U) & ~(uintptr_t) 63U);
    float *right = (float *) (((uintptr_t) right_raw + 63U) & ~(uintptr_t) 63U);

    /* Run both aligned and unaligned input pairs, including arbitrary tails. */
    for (int trial = 0; trial < TRIALS; ++trial) {
        for (int i = 0; i <= MAX_LENGTH; ++i) {
            left[i] = ((float) (next_random() & 0xffff) / 4096.0f) - 8.0f;
            right[i] = ((float) (next_random() & 0xffff) / 4096.0f) - 8.0f;
        }
        for (int offset = 0; offset <= 1; ++offset) {
            for (int length = 0; length <= MAX_LENGTH; ++length) {
                const float expected = reference_ed(left + offset, right + offset, length);
                const float actual = ts_euclidean_distance_SIMD(left + offset, right + offset, length, FLT_MAX);
                if (!nearly_equal(expected, actual)) {
                    fprintf(stderr, "ED mismatch (offset=%d length=%d): reference=%g SIMD=%g\n",
                            offset, length, expected, actual);
                    free(right_raw);
                    free(left_raw);
                    return 1;
                }

                if (expected > 0.0f) {
                    const float bound = expected * 0.5f;
                    const float abandoned = ts_euclidean_distance_SIMD(left + offset, right + offset, length, bound);
                    if (abandoned < bound) {
                        fprintf(stderr, "ED did not respect bound (offset=%d length=%d): bound=%g result=%g\n",
                                offset, length, bound, abandoned);
                        free(right_raw);
                        free(left_raw);
                        return 1;
                    }
                }
            }
        }
    }
    free(right_raw);
    free(left_raw);
    return 0;
}
