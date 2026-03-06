//
//  dft.c
//  sfa C version for MESSI
//
//  Based on dft code by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "config.h"
#include "globals.h"
#include "ads/isax_index.h"
#include "math.h"
#if defined(__x86_64__)
#include "x86intrin.h"
#include "immintrin.h"
#endif
#include <fftw3.h>
#include "ads/sfa/dft.h"

#if defined(__x86_64__)
__attribute__((target("avx512f")))
int transform_vec_SIMD512F(ts_type norm_factor, int i, int coeff_number, ts_type *transform) {
    __m512 norm_factor_vec_h = _mm512_set1_ps(norm_factor);
    __m512 transform_vec_h = _mm512_loadu_ps(&transform[i]);

    for (; i < coeff_number - (coeff_number % 16); i += 16) {
        _mm512_storeu_ps(&transform[i],
                _mm512_mul_ps(transform_vec_h, norm_factor_vec_h));
        transform_vec_h = _mm512_loadu_ps(&transform[i + 16]);
    }

    return i;
}
#endif

/*
    This function calculates the FFT coefficients for a given time series
*/
void fft_from_ts(
        isax_index *index, ts_type *ts,
        int coeff_number, int best_only,
        fftwf_complex *ts_out, ts_type *transform, fftwf_plan plan_forward) {
    unsigned long ts_length = index->settings->timeseries_size;
    fftwf_execute(plan_forward);

    // Image part of first (DC) coefficient
    ts_out[0][1] = 0;
    int j = 0;

    // if normalized, ignore first coeff and start with offset 1
    int start_offset = index->settings->is_norm ? 1 : 0;

    if (best_only) {
        for (int k = 0; k < coeff_number / 2; ++k, j += 2) {
            int coeff = index->coefficients[k] + start_offset;
            transform[j] = ts_out[coeff][0];
            transform[j + 1] = ts_out[coeff][1] * -1;
        }
    } else {
        for (int k = start_offset; k < coeff_number / 2 + start_offset; ++k, j += 2) {
            transform[j] = ts_out[k][0];
            transform[j + 1] = ts_out[k][1] * -1;
        }
    }

    // normalizing fft result in frequency domain to allow for lower bounding
    ts_type norm_factor = index->norm_factor;
    int i = 0;

#if defined(__x86_64__)
    if (index->settings->SIMD_flag) {
//        if (__builtin_cpu_supports("avx512f")) {  // this actually is slower than scalar...
//            i = transform_vec_SIMD512F(norm_factor, i, coeff_number, transform);
//        }

        __m256 norm_factor_vec = _mm256_set1_ps(norm_factor);
        __m256 transform_vec = _mm256_loadu_ps(&transform[i]);

        for (; i < coeff_number - (coeff_number % 8); i += 8) {
            _mm256_storeu_ps(&transform[i], _mm256_mul_ps(transform_vec, norm_factor_vec));
            transform_vec = _mm256_loadu_ps(&transform[i + 8]);
        }
    }
#endif

    for (; i < coeff_number; ++i) {
        transform[i] *= norm_factor;
    }

    return;
}

/*
    This function discretized FFT coefficients with the intervals from MCB
    The current transform is pointed to by dft_mem_array
*/
void sfa_from_fft(isax_index *index, ts_type *cur_transform, unsigned char *cur_sfa_word) {
    int paa_segments = index->settings->paa_segments;
    int cardinality = index->settings->sax_alphabet_cardinality - 1;

    for (int k = 0; k < paa_segments; ++k) {
        unsigned int c;

        for (c = 0; c < cardinality; c++) {
            if (cur_transform[k] < index->bins[k][c]) {
                break;
            }
        }

        cur_sfa_word[k] = (unsigned char) (c);
    }
}

/*
    This function discretized FFT coefficients with the intervals from MCB
    The current transform is pointed to by dft_mem_array
    A sax-cardinality as a multiple of 8 is assumed
*/
#if defined(__x86_64__)
void sfa_from_fft_SIMD(isax_index *index, ts_type *cur_transform, unsigned char *cur_sfa_word) {
    int paa_segments = index->settings->paa_segments;
    int cardinality = index->settings->sax_alphabet_cardinality;

    for (int k = 0; k < paa_segments; ++k) {
        int idx = cardinality;
        __m256 comparison;

        for (int count = 0; count < cardinality; count += 8) {
            comparison = _mm256_cmp_ps(
                _mm256_set1_ps(cur_transform[k]),
                _mm256_loadu_ps(&index->bins[k][count]),
                17);
            int compMask = _mm256_movemask_ps(comparison);
            idx = count + _tzcnt_u32(compMask);

            if (_popcnt32(compMask)) {
                break;
            }
        }

        cur_sfa_word[k] = (unsigned char)(idx > cardinality - 1 ? cardinality - 1 : idx);
    }
}
#endif

#if defined(__x86_64__)
__attribute__((target("avx512vl")))
void sfa_from_fft_SIMD512_256b(isax_index *index, ts_type *cur_transform, unsigned char *cur_sfa_word) {
    int paa_segments = index->settings->paa_segments;
    int cardinality = index->settings->sax_alphabet_cardinality;

    for (int k = 0; k < paa_segments; ++k) {
        int idx = cardinality;
        __m256 comparison;

        for (int count = 0; count < cardinality; count += 8) {
            int compMask = _mm256_cmp_ps_mask(
                _mm256_set1_ps(cur_transform[k]),
                _mm256_loadu_ps(&index->bins[k][count]),
                17);
            idx = count + (_tzcnt_u32(compMask));

            if (_popcnt32(compMask)) {
                break;
            }
        }

        cur_sfa_word[k] = (unsigned char)(idx > cardinality - 1 ? cardinality - 1 : idx);
    }
}
#endif

#if defined(__x86_64__)
__attribute__((target("avx512vl")))
void sfa_from_fft_SIMD512_512b(isax_index *index, ts_type *cur_transform, unsigned char *cur_sfa_word) {
    int paa_segments = index->settings->paa_segments;
    int cardinality = index->settings->sax_alphabet_cardinality;

    for (int k = 0; k < paa_segments; ++k) {
        int idx = cardinality;
        __m256 comparison;

        for (int count = 0; count < cardinality; count += 16) {
            int compMask = _mm512_cmp_ps_mask(
                _mm512_set1_ps(cur_transform[k]),
                _mm512_loadu_ps(&index->bins[k][count]),
                17);
            idx = count + (_tzcnt_u32(compMask));

            if (_popcnt32(compMask)) {
                break;
            }
        }

        cur_sfa_word[k] = (unsigned char)(idx > cardinality - 1 ? cardinality - 1 : idx);
    }
}
#endif

/*
    This function creates an SFA representation of a time series 
*/
enum response sfa_from_ts(isax_index *index, ts_type *ts_in, sax_type *sax_out, fftwf_complex *ts_out, ts_type *transform,
            fftwf_plan plan_forward) {
    int use_best = index->settings->coeff_number != 0;
    fft_from_ts(index, ts_in, index->settings->paa_segments, use_best, ts_out, transform, plan_forward);
    ts_type *cur_coeff_line = (ts_type*)calloc(index->settings->paa_segments, sizeof(ts_type));

    for (int i = 0; i < index->settings->paa_segments; ++i) {
        cur_coeff_line[i] = transform[i];
    }

#if defined(__x86_64__)
    if (index->settings->SIMD_flag && __builtin_cpu_supports("avx512vl")) {
        if (index->settings->sax_alphabet_cardinality >= 16) {
            sfa_from_fft_SIMD512_512b(index, cur_coeff_line, sax_out);
        } else {
            sfa_from_fft_SIMD512_256b(index, cur_coeff_line, sax_out);
        }
    } else if (index->settings->SIMD_flag) {
        sfa_from_fft_SIMD(index, cur_coeff_line, sax_out);
    } else {
#endif
        sfa_from_fft(index, cur_coeff_line, sax_out);
#if defined(__x86_64__)
    }
#endif

    free(cur_coeff_line);

    if (sax_out != NULL) {
        return SUCCESS;
    } else {
        fprintf(stderr, "SFA error");
    }
}
