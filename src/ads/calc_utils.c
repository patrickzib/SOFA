#include "../../config.h"
#include "../../globals.h"

#include <stdio.h>
#include "math.h"
#if defined(__x86_64__)
#include "immintrin.h"
#endif

#include "ads/calc_utils.h"

////// Utility functions ////
// Function to calculate mean of an array of floats
float calculateMean(ts_type *data, int n) {
    float sum = 0.0;

    for (int i = 0; i < n; i++) {
        sum += data[i];
    }

    return sum / n;
}

#if defined(__x86_64__)
float calculateMean_SIMD(ts_type *data, int n) {
    float sum = 0.0;
    float sum_a[8];
    __m256 sum_v = _mm256_setzero_ps();

    int i = 0;

    for (; i < n - (n % 8); i += 8) {
        sum_v = _mm256_add_ps(_mm256_loadu_ps(&data[i]), sum_v);
    }

    _mm256_storeu_ps(sum_a, sum_v);
    sum = sum_a[0] + sum_a[1] + sum_a[2] + sum_a[3] + sum_a[4] + sum_a[5] + sum_a[6] + sum_a[7];

    for (; i < n; i++) {
        sum += data[i];
    }

    return sum / n;
}
#endif

// Function to calculate standard deviation of an array of floats
float calculateStdDev(ts_type *data, int n, ts_type mean) {
    ts_type sumSquaredDiffs = 0.0;

    for (int i = 0; i < n; i++) {
        sumSquaredDiffs += (data[i] - mean) * (data[i] - mean);
    }

    return sqrt(sumSquaredDiffs / n);
}

#if defined(__x86_64__)
float calculateStdDev_SIMD(ts_type *data, int n, ts_type mean) {
    ts_type sumSquaredDiffs = 0.0;
    int i = 0;

    __m256 datavec;
    __m256 meanvec = _mm256_set1_ps(mean);
    __m256 squaredDiffVec = _mm256_setzero_ps();
    float sumSquaredDiffVec[8];

    for (; i < n - (n % 8); i += 8) {
        datavec = _mm256_loadu_ps(&data[i]);
        datavec = _mm256_sub_ps(datavec, meanvec);
        squaredDiffVec = _mm256_fmadd_ps(datavec, datavec, squaredDiffVec);
    }

    squaredDiffVec = _mm256_hadd_ps(squaredDiffVec, squaredDiffVec);
    _mm256_storeu_ps(sumSquaredDiffVec,
        _mm256_hadd_ps(squaredDiffVec, squaredDiffVec));
    sumSquaredDiffs = sumSquaredDiffVec[0] + sumSquaredDiffVec[4];

    for (; i < n; i++) {
        sumSquaredDiffs += (data[i] - mean) * (data[i] - mean);
    }

    return sqrt(sumSquaredDiffs / n);
}
#endif

#if defined(__x86_64__)
__attribute__((target("avx512f")))
float calculateStdDev_SIMD512F(ts_type *data, int n, ts_type mean) {
    ts_type sumSquaredDiffs = 0.0;
    int i = 0;

    __m512 datavec_v;
    __m512 meanvec_v = _mm512_set1_ps(mean);
    __m512 squaredDiffVec_v = _mm512_setzero_ps();

    for (; i < n - (n % 16); i += 16) {
        datavec_v = _mm512_loadu_ps(&data[i]);
        datavec_v = _mm512_sub_ps(datavec_v, meanvec_v);
        squaredDiffVec_v = _mm512_fmadd_ps(datavec_v, datavec_v, squaredDiffVec_v);
    }

    sumSquaredDiffs = _mm512_reduce_add_ps(squaredDiffVec_v);

    __m256 datavec;
    __m256 meanvec = _mm256_set1_ps(mean);
    __m256 squaredDiffVec = _mm256_setzero_ps();
    float sumSquaredDiffVec[8];

    for (; i < n - (n % 8); i += 8) {
        datavec = _mm256_loadu_ps(&data[i]);
        datavec = _mm256_sub_ps(datavec, meanvec);
        squaredDiffVec = _mm256_fmadd_ps(datavec, datavec, squaredDiffVec);
    }

    squaredDiffVec = _mm256_hadd_ps(squaredDiffVec, squaredDiffVec);
    _mm256_storeu_ps(sumSquaredDiffVec,
            _mm256_hadd_ps(squaredDiffVec, squaredDiffVec));
    sumSquaredDiffs += sumSquaredDiffVec[0] + sumSquaredDiffVec[4];

    for (; i < n; i++) {
        sumSquaredDiffs += (data[i] - mean) * (data[i] - mean);
    }

    return sqrt(sumSquaredDiffs / n);
}
#endif

// Function to perform zero mean normalization
void znorm(ts_type *data, int n) {
    // Calculate mean
    ts_type mean = calculateMean(data, n);

    // Calculate standard deviation
    ts_type stdDev = calculateStdDev(data, n, mean);

    // Normalize each data point
    if (stdDev < 1e-8) {
        stdDev = 1.0;
    }

    for (int i = 0; i < n; i++) {
        data[i] = (data[i] - mean) / stdDev;
    }
}

#if defined(__x86_64__)
void znorm_SIMD(ts_type *data, int n) {
    // Calculate mean
    ts_type mean = calculateMean_SIMD(data, n);

    // Calculate standard deviation
    ts_type stdDev = calculateStdDev_SIMD(data, n, mean);

    // Normalize each data point
    if (stdDev < 1e-8) {
        stdDev = 1.0;
    }

    int i = 0;

    __m256 stdDev_v = _mm256_set1_ps(stdDev);
    __m256 mean_v = _mm256_set1_ps(mean);

    for (; i < n - (n % 8); i += 8) {
        _mm256_storeu_ps(&data[i], _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(&data[i]), mean_v), stdDev_v));
    }

    for (; i < n; i++) {
        data[i] = (data[i] - mean) / stdDev;
    }
}
#endif

#if defined(__x86_64__)
__attribute__((target("avx512f")))
void znorm_SIMD512F(ts_type *data, int n) {
    // Calculate mean
    ts_type mean = calculateMean_SIMD(data, n);

    // Calculate standard deviation
    ts_type stdDev = calculateStdDev_SIMD512F(data, n, mean);

    // Normalize each data point
    if (stdDev < 1e-8) {
        stdDev = 1.0;
    }

    int i = 0;

    __m512 stdDev_vv = _mm512_set1_ps(stdDev);
    __m512 mean_vv = _mm512_set1_ps(mean);

    for (; i < n - (n % 16); i += 16) {
        _mm512_storeu_ps(&data[i], _mm512_div_ps(_mm512_sub_ps(_mm512_loadu_ps(&data[i]), mean_vv), stdDev_vv));
    }

    __m256 stdDev_v = _mm256_set1_ps(stdDev);
    __m256 mean_v = _mm256_set1_ps(mean);

    for (; i < n - (n % 8); i += 8) {
        _mm256_storeu_ps(&data[i], _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(&data[i]), mean_v), stdDev_v));
    }

    for (; i < n; i++) {
        data[i] = (data[i] - mean) / stdDev;
    }
}
#endif
