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
#include "immintrin.h"
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

float dot_product(ts_type *X, ts_type *Y, int m) {
    float dot = 0;
    for (int i = 0; i < m; i++) {
        dot += X[i] * Y[i];
    }
    return dot;
}

float ts_euclidean_distance_dot_product(ts_type *X, ts_type *Y, int size) {
    float dot = dot_product(X, Y, size);
    float result = 2.0f * (float) size * (1.0f - (dot / (float) size));
    // result = sqrtf(result);
    return result;
}

float dot_product_simd(ts_type *X, ts_type *Y, int m) {
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
    while (i < m) {
        dot += X[i] * Y[i];
        i++;
    }

    return dot;
}

float ts_euclidean_distance_dot_product_SIMD(float *X, float *Y, int size) {
    float dot = dot_product_simd(X, Y, size);
    float result = 2.0f * (float) size * (1.0f - (dot / (float) size));
    // result= sqrtf(result);
    return result;
}


float ts_euclidean_distance(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0;
    while (size > 0 && distance < bound) {
        size--;
        distance += (t[size] - s[size]) * (t[size] - s[size]);
    }
//    distance = sqrtf(distance);
    return distance;
}

float ts_euclidean_distance_SIMD(ts_type *t, ts_type *s, int size, float bound) {
    float distance = 0;
    int i = 0;
    float distancef[8];
    int size2 = size;

    __m256 v_t, v_s, v_d, distancev;
    while (size >= 8 && distance < bound) {
        v_t = _mm256_loadu_ps(&t[i]);
        v_s = _mm256_loadu_ps(&s[i]);

        v_d = _mm256_sub_ps(v_t, v_s);
        v_d = _mm256_mul_ps(v_d, v_d);
        size -= 8;

        i = i + 8;
        distancev = _mm256_hadd_ps(v_d, v_d);
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


float ts_ed(ts_type * t, ts_type * s, int size, float bound, char is_simd, char is_norm) {
    if (is_simd) {
        if (is_norm) {
            return ts_euclidean_distance_dot_product_SIMD(t, s, size);
        }
        else {
            return ts_euclidean_distance_SIMD(t, s, size, bound);
        }
    }
    else {
        if (is_norm) {
            return ts_euclidean_distance_dot_product(t, s, size);
        }
        else {
            return ts_euclidean_distance(t, s, size, bound);
        }
    }

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
