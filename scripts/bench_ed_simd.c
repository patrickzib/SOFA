/* Build with `cmake --build build --target bench_ed_simd` and run the result
   on the CPU of interest.  It reports scalar and SIMD squared-ED throughput. */
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ads/sax/ts.h"

static double elapsed_seconds(clock_t start, clock_t end) {
    return (double) (end - start) / (double) CLOCKS_PER_SEC;
}

static double benchmark(float (*function)(float *, float *, int, float),
                        float *left, float *right, int length, float bound, int iterations) {
    volatile float sink = 0.0f;
    const clock_t start = clock();
    for (int i = 0; i < iterations; ++i) sink += function(left, right, length, bound);
    const double seconds = elapsed_seconds(start, clock());
    if (sink == -1.0f) fprintf(stderr, "ignore: %g\n", sink);
    return seconds;
}

int main(int argc, char **argv) {
    const int length = argc > 1 ? atoi(argv[1]) : 256;
    const int iterations = argc > 2 ? atoi(argv[2]) : 1000000;
    if (length <= 0 || iterations <= 0) return 2;
    float *left = malloc((size_t) length * sizeof(*left));
    float *right = malloc((size_t) length * sizeof(*right));
    if (left == NULL || right == NULL) return 1;
    for (int i = 0; i < length; ++i) {
        left[i] = (float) (i % 17);
        right[i] = (float) ((i * 7) % 19);
    }

    const double scalar = benchmark(ts_euclidean_distance, left, right, length, FLT_MAX, iterations);
    const double simd = benchmark(ts_euclidean_distance_SIMD, left, right, length, FLT_MAX, iterations);
    printf("length=%d iterations=%d scalar=%.3f ns/call SIMD=%.3f ns/call speedup=%.2fx\n",
           length, iterations, scalar * 1e9 / iterations, simd * 1e9 / iterations, scalar / simd);
    free(right);
    free(left);
    return 0;
}
