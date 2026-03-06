//
//  ts.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "../../../config.h"
#include "../../../globals.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#if defined(__x86_64__)
#include "immintrin.h"
#endif
#include "ads/sax/ts.h"

/**
 This function converts a string of floats seperated by a delimeter into a ts 
 record of a size ts_size.
 @param char ts_str[]
 @param int ts_size
 @param const char * delims
 @return *ts
 */
void ts_parse_str(char ts_str[], ts_type *ts_out, int ts_size, const char *delims) {
    int index = 0;
    char *result = strtok(ts_str, delims);

    while (result != NULL) {
        ts_out[index] = atof(result);
        result = strtok(NULL, delims);

#ifdef SANITY_CHECK
        if (index >= ts_size) {
            fprintf(stderr, "sanity error: Time series bigger than limit of %d", ts_size);

            exit(-1);
        }
#endif
        index++;
    }

    free(result);
}

float ts_euclidean_distance_dot_product(ts_type *X, ts_type *Y, int size) {
    float dot = 0;

    for (int count = 0; count < size; count++) {
        dot += X[count] * Y[count];
    }

    return 2.0f * (float)size * (1.0f - (dot / (float)size));
}

#if defined(__x86_64__)
float ts_euclidean_distance_dot_product_SIMD(float *X, float *Y, int m) {
    __m256 sum_vec = _mm256_setzero_ps();  // Initialize sum vector to zero
    int count = 0;

    // Process 8 elements at a time
    for (; count < m - 7; count += 8) {
        __m256 x_vec = _mm256_loadu_ps(&X[count]);
        __m256 y_vec = _mm256_loadu_ps(&Y[count]);
        sum_vec = _mm256_fmadd_ps(x_vec, y_vec, sum_vec);
    }

    // Sum the elements of the sum_vec
    float temp[8];
    _mm256_storeu_ps(temp, sum_vec);
    float dot = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];

    // Remaining values, if length is not divisible by 8!
    for (; count < m; count++) {
        dot += X[count] * Y[count];
    }

    float size = (float)m;

    return 2.0f * size * (1.0f - (dot / size));
}
#endif

float ts_euclidean_distance(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0;

    while (size > 0 && distance < bound) {
        size--;
        distance += (t[size] - s[size]) * (t[size] - s[size]);
    }

    return distance;
}

#if defined(__x86_64__)
float ts_euclidean_distance_SIMD(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0;
    float distancef[8];
    int count = 0;
    int size2 = size;

    __m256 v_t;
    __m256 v_s;
    __m256 v_d;
    __m256 distancev = _mm256_setzero_ps();

    while (size >= 8 && distance < bound) {
        v_t = _mm256_loadu_ps(&t[count]);
        v_s = _mm256_loadu_ps(&s[count]);
        v_d = _mm256_sub_ps(v_t, v_s);
        distancev = _mm256_fmadd_ps(v_d, v_d, distancev);
        size -= 8;
        count += 8;
        __m256 intermed_dist = _mm256_hadd_ps(distancev, distancev);
        intermed_dist = _mm256_hadd_ps(intermed_dist, intermed_dist);
        _mm256_storeu_ps(distancef, intermed_dist);
        distance = distancef[0] + distancef[4];
    }

    // Remaining values, if length is not divisible by 8!
    while (count < size2 && distance < bound) {
        distance += (t[count] - s[count]) * (t[count] - s[count]);
        count++;
    }

    return distance;
}
#endif

#if defined(__x86_64__)
__attribute__((target("avx512dq")))
float ts_euclidean_distance_SIMD512DQ(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0;
    int count = 0;
    float distancef[8];
    int size2 = size;

    __m512 vv_t;
    __m512 vv_s;
    __m512 vv_d;
    __m512 sumv = _mm512_setzero_ps();

    while (size >= 16 && distance < bound) {
        vv_t = _mm512_loadu_ps(&t[count]);
        vv_s = _mm512_loadu_ps(&s[count]);
        vv_d = _mm512_sub_ps(vv_t, vv_s);
        sumv = _mm512_fmadd_ps(vv_d, vv_d, sumv);
        size -= 16;
        count += 16;
        distance = _mm512_reduce_add_ps(sumv);
    }

    __m256 v_t;
    __m256 v_s;
    __m256 v_d;
    __m256 distancev = _mm256_setzero_ps();

    if (size >= 8 && distance < bound) {
        v_t = _mm256_loadu_ps(&t[count]);
        v_s = _mm256_loadu_ps(&s[count]);
        v_d = _mm256_sub_ps(v_t, v_s);
        distancev = _mm256_fmadd_ps(v_d, v_d, distancev);
        count += 8;
        distancev = _mm256_hadd_ps(distancev, distancev);
        distancev = _mm256_hadd_ps(distancev, distancev);
        _mm256_storeu_ps(distancef, distancev);
        distance += distancef[0] + distancef[4];
    }

    // Remaining values, if length is not divisible by 16!
    while (count < size2 && distance < bound) {
        distance += (t[count] - s[count]) * (t[count] - s[count]);
        count++;
    }

    return distance;
}
#endif

float ts_ed(ts_type * t, ts_type * s, int size, float bound, char is_simd, char is_norm) {
#if defined(__x86_64__)
    if (is_simd) {
        return __builtin_cpu_supports("avx512dq") ?
                ts_euclidean_distance_SIMD512DQ(t, s, size, bound) :
                ts_euclidean_distance_SIMD(t, s, size, bound);
    } else {
#endif
        return ts_euclidean_distance(t, s, size, bound);
#if defined(__x86_64__)
    }
#endif
}

/** 
 This function prints a ts record of a size.
 @param ts *ts
 @param int size
*/
void ts_print(ts_type *ts, int size) {
    for (int count = 0; count < size; count++) {
        printf("%lf", ts[count]);
    }

    printf("\n");
}
