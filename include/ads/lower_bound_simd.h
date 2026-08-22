#ifndef MESSI_LOWER_BOUND_SIMD_H
#define MESSI_LOWER_BOUND_SIMD_H

#include "config.h"
#include "../../globals.h"
#include "ads/isax_index.h"
#include <stdint.h>

#if defined(__AVX512F__)
#include <immintrin.h>
#elif ADS_HAVE_AVX2
#include <immintrin.h>
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

static inline void messi_lower_bound_16_bounds(const isax_index *index,
                                                const sax_type *sax,
                                                const sax_type *cardinalities,
                                                float lower[16], float upper[16]) {
    const int max_bits = index->settings->sax_bit_cardinality;
    const int alphabet = index->settings->sax_alphabet_cardinality;
    const int stride = alphabet - 1;
    for (int i = 0; i < 16; ++i) {
        const int shift = max_bits - cardinalities[i];
        const int low = (sax[i] >> shift) << shift;
        const int high = low | (shift == 0 ? 0 : ((1 << shift) - 1));
        lower[i] = low == 0 ? MINVAL : index->binsv[i * stride + low - 1];
        upper[i] = high == alphabet - 1 ? MAXVAL : index->binsv[i * stride + high];
    }
}

static inline float messi_lower_bound_16_distance(const float values[16],
                                                   const float lower[16],
                                                   const float upper[16],
                                                   float bsf, float factor) {
#if defined(__AVX512F__)
    const __m512 zero = _mm512_setzero_ps();
    const __m512 v0 = ((uintptr_t) values & 63U) == 0 ? _mm512_load_ps(values) : _mm512_loadu_ps(values);
    const __m512 lo0 = _mm512_loadu_ps(lower);
    const __m512 hi0 = _mm512_loadu_ps(upper);
    __mmask16 below = _mm512_cmp_ps_mask(lo0, v0, _CMP_GT_OS);
    __mmask16 inside = _mm512_kand(_mm512_cmp_ps_mask(lo0, v0, _CMP_LE_OS),
                                   _mm512_cmp_ps_mask(hi0, v0, _CMP_LT_OS));
    __m512 dlo = _mm512_sub_ps(lo0, v0);
    __m512 dhi = _mm512_sub_ps(hi0, v0);
 #if defined(__FMA__)
    const __m512 low_squared = _mm512_fmadd_ps(dlo, dlo, zero);
    const __m512 high_squared = _mm512_fmadd_ps(dhi, dhi, zero);
 #else
    const __m512 low_squared = _mm512_mul_ps(dlo, dlo);
    const __m512 high_squared = _mm512_mul_ps(dhi, dhi);
 #endif
    __m512 squared = _mm512_maskz_mov_ps(below, low_squared);
    squared = _mm512_mask_mov_ps(squared, inside, high_squared);
    const float first = _mm512_reduce_add_ps(_mm512_maskz_mov_ps(0x00ff, squared)) * factor;
    if (first > bsf) return first;
    return (first + _mm512_reduce_add_ps(_mm512_maskz_mov_ps(0xff00, squared)) * factor);
#elif ADS_HAVE_AVX2
    const __m256 zero = _mm256_setzero_ps();
    const __m256 v0 = ((uintptr_t) values & 31U) == 0 ? _mm256_load_ps(values) : _mm256_loadu_ps(values);
    const __m256 v1 = ((uintptr_t) &values[8] & 31U) == 0 ? _mm256_load_ps(&values[8]) : _mm256_loadu_ps(&values[8]);
    const __m256 lo0 = _mm256_loadu_ps(lower);
    const __m256 lo1 = _mm256_loadu_ps(&lower[8]);
    const __m256 hi0 = _mm256_loadu_ps(upper);
    const __m256 hi1 = _mm256_loadu_ps(&upper[8]);
    const __m256 d00 = _mm256_sub_ps(lo0, v0);
    const __m256 d01 = _mm256_sub_ps(hi0, v0);
    const __m256 d10 = _mm256_sub_ps(lo1, v1);
    const __m256 d11 = _mm256_sub_ps(hi1, v1);
    const __m256 m00 = _mm256_cmp_ps(lo0, v0, _CMP_GT_OS);
    const __m256 m01 = _mm256_and_ps(_mm256_cmp_ps(lo0, v0, _CMP_LE_OS), _mm256_cmp_ps(hi0, v0, _CMP_LT_OS));
    const __m256 m10 = _mm256_cmp_ps(lo1, v1, _CMP_GT_OS);
    const __m256 m11 = _mm256_and_ps(_mm256_cmp_ps(lo1, v1, _CMP_LE_OS), _mm256_cmp_ps(hi1, v1, _CMP_LT_OS));
#if defined(__FMA__)
    __m256 sum0 = _mm256_and_ps(m00, _mm256_fmadd_ps(d00, d00, zero));
    sum0 = _mm256_or_ps(sum0, _mm256_and_ps(m01, _mm256_fmadd_ps(d01, d01, zero)));
    __m256 sum1 = _mm256_and_ps(m10, _mm256_fmadd_ps(d10, d10, zero));
    sum1 = _mm256_or_ps(sum1, _mm256_and_ps(m11, _mm256_fmadd_ps(d11, d11, zero)));
#else
    __m256 sum0 = _mm256_or_ps(_mm256_and_ps(m00, _mm256_mul_ps(d00, d00)), _mm256_and_ps(m01, _mm256_mul_ps(d01, d01)));
    __m256 sum1 = _mm256_or_ps(_mm256_and_ps(m10, _mm256_mul_ps(d10, d10)), _mm256_and_ps(m11, _mm256_mul_ps(d11, d11)));
#endif
    sum0 = _mm256_hadd_ps(sum0, sum0);
    sum0 = _mm256_hadd_ps(sum0, sum0);
    float first_result[8];
    _mm256_storeu_ps(first_result, sum0);
    const float first = (first_result[0] + first_result[4]) * factor;
    if (first > bsf) return first;
    sum1 = _mm256_hadd_ps(sum1, sum1);
    sum1 = _mm256_hadd_ps(sum1, sum1);
    float second_result[8];
    _mm256_storeu_ps(second_result, sum1);
    return first + (second_result[0] + second_result[4]) * factor;
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    float first = 0.0f;
    for (int i = 0; i < 8; i += 4) {
        float32x4_t value = vld1q_f32(&values[i]);
        float32x4_t low = vld1q_f32(&lower[i]);
        float32x4_t high = vld1q_f32(&upper[i]);
        uint32x4_t below = vcgtq_f32(low, value);
        uint32x4_t inside = vandq_u32(vcleq_f32(low, value), vcltq_f32(high, value));
        float32x4_t dlow = vsubq_f32(low, value);
        float32x4_t dhigh = vsubq_f32(high, value);
        first += vaddvq_f32(vbslq_f32(below, vmulq_f32(dlow, dlow), vdupq_n_f32(0.0f)));
        first += vaddvq_f32(vbslq_f32(inside, vmulq_f32(dhigh, dhigh), vdupq_n_f32(0.0f)));
    }
    first *= factor;
    if (first > bsf) return first;
    float second = 0.0f;
    for (int i = 8; i < 16; i += 4) {
        float32x4_t value = vld1q_f32(&values[i]);
        float32x4_t low = vld1q_f32(&lower[i]);
        float32x4_t high = vld1q_f32(&upper[i]);
        uint32x4_t below = vcgtq_f32(low, value);
        uint32x4_t inside = vandq_u32(vcleq_f32(low, value), vcltq_f32(high, value));
        float32x4_t dlow = vsubq_f32(low, value);
        float32x4_t dhigh = vsubq_f32(high, value);
        second += vaddvq_f32(vbslq_f32(below, vmulq_f32(dlow, dlow), vdupq_n_f32(0.0f)));
        second += vaddvq_f32(vbslq_f32(inside, vmulq_f32(dhigh, dhigh), vdupq_n_f32(0.0f)));
    }
    return first + second * factor;
#else
    float result = 0.0f;
    for (int i = 0; i < 16 && result * factor <= bsf; ++i) {
        const float difference = values[i] < lower[i] ? lower[i] - values[i] :
                                 values[i] > upper[i] ? upper[i] - values[i] : 0.0f;
        result += difference * difference;
    }
    return result * factor;
#endif
}

static inline float messi_lower_bound_16(const isax_index *index,
                                         const float *values,
                                         const sax_type *sax,
                                         const sax_type *cardinalities,
                                         float bsf, float factor) {
    float lower[16], upper[16];
    messi_lower_bound_16_bounds(index, sax, cardinalities, lower, upper);
    return messi_lower_bound_16_distance(values, lower, upper, bsf, factor);
}

#endif
