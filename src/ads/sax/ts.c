//
//  ts.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//
#include "config.h"
#include "../../../globals.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#if ADS_HAVE_AVX2 || defined(__AVX512F__)
#include "immintrin.h"
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
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
        if (index >= ts_size)
        {
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
    for (int i = 0; i < size; i++) {
        dot += X[i] * Y[i];
    }

    float result = 2.0f * (float) size * (1.0f - (dot / (float) size));
    // result = sqrtf(result);
    return result;
}

#if ADS_HAVE_AVX2
float ts_euclidean_distance_dot_product_SIMD(float *X, float *Y, int m) {
    __m256 sum_vec = _mm256_setzero_ps(); // Initialize sum vector to zero

    int i = 0;
    // Process 8 elements at a time
    for (; i < m-7; i += 8) {
        __m256 x_vec = _mm256_loadu_ps(&X[i]);
        __m256 y_vec = _mm256_loadu_ps(&Y[i]);
        __m256 prod_vec = _mm256_mul_ps(x_vec, y_vec);
        sum_vec = _mm256_add_ps(sum_vec, prod_vec);
    }

    // Sum the elements of the sum_vec
    float temp[8];
    _mm256_storeu_ps(temp, sum_vec);
    float dot = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];

    // Remaining values, if length is not divisible by 8!
    for (; i < m; i++) {
        dot += X[i] * Y[i];
    }

    float size = (float) m;
    float result = 2.0f * size * (1.0f - (dot / size));
    // result= sqrtf(result);
    return result;
}
#else
float ts_euclidean_distance_dot_product_SIMD(float *X, float *Y, int m) {
    return ts_euclidean_distance_dot_product(X, Y, m);
}
#endif

float ts_euclidean_distance(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0;
    while (size > 0 && distance < bound) {
        size--;
        distance += (t[size] - s[size]) * (t[size] - s[size]);
    }
    //    distance = sqrtf(distance);
    return distance;
}

#if defined(__AVX512F__)
float ts_euclidean_distance_SIMD(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0.0f;
    const int original_size = size;
    int i = 0;
    while (size >= 16 && distance < bound) {
        const int aligned = (((uintptr_t) &t[i] | (uintptr_t) &s[i]) & 63U) == 0;
        __m512 vt = aligned ? _mm512_load_ps(&t[i]) : _mm512_loadu_ps(&t[i]);
        __m512 vs = aligned ? _mm512_load_ps(&s[i]) : _mm512_loadu_ps(&s[i]);
        __m512 diff = _mm512_sub_ps(vt, vs);
        __m512 squared = _mm512_mul_ps(diff, diff);
        distance += _mm512_reduce_add_ps(squared);
        i += 16;
        size -= 16;
    }
    while (i < original_size && distance < bound) {
        distance += (t[i] - s[i]) * (t[i] - s[i]);
        ++i;
    }
    return distance;
}
#elif ADS_HAVE_AVX2
float ts_euclidean_distance_SIMD(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0;
    int i = 0;
    float distancef[8];
    int size2 = size;

    __m256 v_t, v_s, v_d, distancev;
    while (size >= 8 && distance < bound) {
        const int aligned = (((uintptr_t) &t[i] | (uintptr_t) &s[i]) & 31U) == 0;
        v_t = aligned ? _mm256_load_ps(&t[i]) : _mm256_loadu_ps(&t[i]);
        v_s = aligned ? _mm256_load_ps(&s[i]) : _mm256_loadu_ps(&s[i]);

        v_d = _mm256_sub_ps(v_t, v_s);
#if defined(__FMA__)
        distancev = _mm256_fmadd_ps(v_d, v_d, _mm256_setzero_ps());
#else
        v_d = _mm256_mul_ps(v_d, v_d);
        distancev = v_d;
#endif
        size -= 8;

        i = i + 8;
        distancev = _mm256_hadd_ps(distancev, distancev);
        distancev = _mm256_hadd_ps(distancev, distancev);
        _mm256_storeu_ps(distancef, distancev);
        distance += distancef[0] + distancef[4];
    }

    // Remaining values, if length is not divisible by 8!
    while (i < size2 && distance < bound) {
        distance += (t[i] - s[i]) * (t[i] - s[i]);
        i++;
    }

    //    distance = sqrtf(distance);
    return distance;
}
#else
float ts_euclidean_distance_SIMD(ts_type *t, ts_type *s, int size, float bound) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    float distance = 0.0f;
    int i = 0;
    while (size >= 4 && distance < bound) {
        float32x4_t vt = vld1q_f32(&t[i]);
        float32x4_t vs = vld1q_f32(&s[i]);
        float32x4_t diff = vsubq_f32(vt, vs);
        distance += vaddvq_f32(vmulq_f32(diff, diff));
        i += 4;
        size -= 4;
    }
    while (size > 0 && distance < bound) {
        distance += (t[i] - s[i]) * (t[i] - s[i]);
        ++i;
        --size;
    }
    return distance;
#else
    return ts_euclidean_distance(t, s, size, bound);
#endif
}
#endif


float ts_ed(ts_type * t, ts_type * s, int size, float bound) {
    return ts_euclidean_distance_SIMD(t, s, size, bound);
}


/** 
 This function prints a ts record of a size.
 @param ts *ts
 @param int size
*/
void ts_print(ts_type *ts, int size) {
    int i;
    for (i = 0; i < size; i++) {
        printf("%lf", ts[i]);
    }
    printf("\n");
}
