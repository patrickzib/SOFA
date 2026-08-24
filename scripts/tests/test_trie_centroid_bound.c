#include <float.h>
#include <math.h>
#include <stdio.h>

#include "ads/sax/ts.h"

/* Validate the exact rule used by leaf-centroid pruning.  A reported prune
 * must never discard a member whose squared ED is below the current BSF. */
static unsigned int state = 1;

static float next_value(void) {
    state = state * 1664525U + 1013904223U;
    return ((float) (state & 0xffffU) / 32768.0f) - 1.0f;
}

int main(void) {
    enum { MEMBERS = 37, LENGTH = 128, TRIALS = 200 };
    ts_type records[MEMBERS][LENGTH], centroid[LENGTH], query[LENGTH];
    for (int trial = 0; trial < TRIALS; ++trial) {
        for (int i = 0; i < MEMBERS; ++i)
            for (int d = 0; d < LENGTH; ++d)
                records[i][d] = next_value() + (float) (i % 5) * 0.2f;
        for (int d = 0; d < LENGTH; ++d) {
            double sum = 0.0;
            for (int i = 0; i < MEMBERS; ++i) sum += records[i][d];
            centroid[d] = (ts_type) (sum / MEMBERS);
            query[d] = next_value() * 2.0f;
        }
        double radius_double = 0.0;
        for (int i = 0; i < MEMBERS; ++i) {
            double distance = 0.0;
            for (int d = 0; d < LENGTH; ++d) {
                const double delta = (double) records[i][d] - centroid[d];
                distance += delta * delta;
            }
            if (distance > radius_double) radius_double = distance;
        }
        const float radius_sq = nextafterf((float) radius_double, INFINITY);
        for (int b = 1; b <= 20; ++b) {
            const float bsf = (float) b * 4.0f;
            float threshold = sqrtf(bsf) + sqrtf(radius_sq);
            threshold = nextafterf(threshold * threshold, INFINITY);
            const float centroid_distance = ts_ed(query, centroid, LENGTH, threshold);
            if (centroid_distance < threshold) continue;
            for (int i = 0; i < MEMBERS; ++i) {
                const float actual = ts_ed(query, records[i], LENGTH, FLT_MAX);
                if (actual < bsf) {
                    fprintf(stderr, "unsafe centroid prune: trial=%d bsf=%g actual=%g\n",
                            trial, bsf, actual);
                    return 1;
                }
            }
        }
    }
    return 0;
}
