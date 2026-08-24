#ifndef MESSI_LOWER_BOUND_SIMD_H
#define MESSI_LOWER_BOUND_SIMD_H

#include "config.h"
#include "../../globals.h"
#include "ads/isax_index.h"
#include <stdint.h>

#ifndef MESSI_RECORD_LB_MAX_DIMENSIONS
#define MESSI_RECORD_LB_MAX_DIMENSIONS 64
#endif

#if defined(__AVX512F__) || ADS_HAVE_AVX2
#include <immintrin.h>
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

/* The values in this table are already squared lower-bound contributions.
 * Gather them by the record's 8-bit symbols; this avoids reconstructing bin
 * intervals for every record in a leaf. */
static inline float messi_record_lb_table_scalar(
        const float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256], const sax_type *word,
        int dimensions, float limit) {
    float distance = 0.0f;
    for (int dimension = 0; dimension < dimensions && distance <= limit; ++dimension)
        distance += table[dimension][word[dimension]];
    return distance;
}

#if ADS_HAVE_AVX2
static inline float messi_record_lb_table_avx2(
        const float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256], const sax_type *word,
        int dimensions, float limit) {
    const float *base = &table[0][0];
    const __m256i offsets = _mm256_setr_epi32(0, 256, 512, 768, 1024, 1280, 1536, 1792);
    float distance = 0.0f;
    int dimension = 0;
    for (; dimension + 8 <= dimensions && distance <= limit; dimension += 8) {
        const __m128i symbols8 = _mm_loadl_epi64((const __m128i *) (word + dimension));
        const __m256i symbols = _mm256_cvtepu8_epi32(symbols8);
        const __m256i indices = _mm256_add_epi32(offsets, symbols);
        const __m256 values = _mm256_i32gather_ps(base + (size_t) dimension * 256, indices, 4);
        __m128 sum = _mm_add_ps(_mm256_castps256_ps128(values), _mm256_extractf128_ps(values, 1));
        sum = _mm_hadd_ps(sum, sum);
        sum = _mm_hadd_ps(sum, sum);
        distance += _mm_cvtss_f32(sum);
    }
    for (; dimension < dimensions && distance <= limit; ++dimension)
        distance += table[dimension][word[dimension]];
    return distance;
}
#endif

#if defined(__AVX512F__)
static inline float messi_record_lb_table_avx512(
        const float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256], const sax_type *word,
        int dimensions, float limit) {
    const float *base = &table[0][0];
    const __m512i offsets = _mm512_setr_epi32(
        0, 256, 512, 768, 1024, 1280, 1536, 1792,
        2048, 2304, 2560, 2816, 3072, 3328, 3584, 3840);
    float distance = 0.0f;
    int dimension = 0;
    for (; dimension + 16 <= dimensions && distance <= limit; dimension += 16) {
        const __m128i symbols16 = _mm_loadu_si128((const __m128i *) (word + dimension));
        const __m512i symbols = _mm512_cvtepu8_epi32(symbols16);
        const __m512i indices = _mm512_add_epi32(offsets, symbols);
        const __m512 values = _mm512_i32gather_ps(indices, base + (size_t) dimension * 256, 4);
        distance += _mm512_reduce_add_ps(values);
    }
    for (; dimension < dimensions && distance <= limit; ++dimension)
        distance += table[dimension][word[dimension]];
    return distance;
}
#endif

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
static inline float messi_record_lb_table_neon(
        const float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256], const sax_type *word,
        int dimensions, float limit) {
    float distance = 0.0f;
    int dimension = 0;
    for (; dimension + 4 <= dimensions && distance <= limit; dimension += 4) {
        const float values[4] = {
            table[dimension][word[dimension]], table[dimension + 1][word[dimension + 1]],
            table[dimension + 2][word[dimension + 2]], table[dimension + 3][word[dimension + 3]]
        };
        distance += vaddvq_f32(vld1q_f32(values));
    }
    for (; dimension < dimensions && distance <= limit; ++dimension)
        distance += table[dimension][word[dimension]];
    return distance;
}
#endif

static inline float messi_record_lb_table_sum(
        const float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256], const sax_type *word,
        int dimensions, float limit, int simd_enabled) {
    if (!simd_enabled)
        return messi_record_lb_table_scalar(table, word, dimensions, limit);
#if defined(__AVX512F__)
    return messi_record_lb_table_avx512(table, word, dimensions, limit);
#elif ADS_HAVE_AVX2
    return messi_record_lb_table_avx2(table, word, dimensions, limit);
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    return messi_record_lb_table_neon(table, word, dimensions, limit);
#else
    return messi_record_lb_table_scalar(table, word, dimensions, limit);
#endif
}

static inline float messi_lower_bound_16_scalar(const isax_index *index,
                                                 const float values[16],
                                                 const sax_type sax[16],
                                                 const sax_type cardinalities[16],
                                                 float bsf, float factor) {
    const int max_bits = index->settings->sax_bit_cardinality;
    const int alphabet = index->settings->sax_alphabet_cardinality;
    const int stride = alphabet - 1;
    float distance = 0.0f;
    for (int i = 0; i < 16 && distance * factor <= bsf; ++i) {
        const int shift = max_bits - cardinalities[i];
        const int low = (sax[i] >> shift) << shift;
        const int high = low | (shift == 0 ? 0 : ((1 << shift) - 1));
        const float lower = low == 0 ? MINVAL : index->binsv[i * stride + low - 1];
        const float upper = high == alphabet - 1 ? MAXVAL : index->binsv[i * stride + high];
        const float difference = values[i] < lower ? lower - values[i] :
                                 values[i] > upper ? upper - values[i] : 0.0f;
        distance += difference * difference;
    }
    return distance * factor;
}

#if ADS_HAVE_AVX2
static inline float messi_lower_bound_avx2_block(const isax_index *index,
                                                  const float *values,
                                                  const sax_type *sax,
                                                  const sax_type *cardinalities,
                                                  int start) {
    const int max_bits = index->settings->sax_bit_cardinality;
    const int alphabet = index->settings->sax_alphabet_cardinality;
    const int stride = alphabet - 1;
    const __m128i symbols8 = _mm_loadl_epi64((const __m128i *) (sax + start));
    const __m128i cards8 = _mm_loadl_epi64((const __m128i *) (cardinalities + start));
    const __m256i symbols = _mm256_cvtepu8_epi32(symbols8);
    const __m256i cards = _mm256_cvtepu8_epi32(cards8);
    const __m256i shifts = _mm256_sub_epi32(_mm256_set1_epi32(max_bits), cards);
    const __m256i low = _mm256_sllv_epi32(_mm256_srlv_epi32(symbols, shifts), shifts);
    const __m256i one = _mm256_set1_epi32(1);
    const __m256i high = _mm256_or_si256(low, _mm256_sub_epi32(_mm256_sllv_epi32(one, shifts), one));
    const __m256i offsets = _mm256_setr_epi32(start * stride, (start + 1) * stride,
                                                (start + 2) * stride, (start + 3) * stride,
                                                (start + 4) * stride, (start + 5) * stride,
                                                (start + 6) * stride, (start + 7) * stride);
    const __m256i zeroi = _mm256_setzero_si256();
    const __m256i lower_zero = _mm256_cmpeq_epi32(low, zeroi);
    const __m256i upper_max = _mm256_cmpeq_epi32(high, _mm256_set1_epi32(alphabet - 1));
    const __m256i lower_index = _mm256_sub_epi32(_mm256_add_epi32(offsets, low), one);
    const __m256i upper_index = _mm256_add_epi32(offsets, high);
    const __m256i safe_lower = _mm256_blendv_epi8(lower_index, zeroi, lower_zero);
    const __m256i safe_upper = _mm256_blendv_epi8(upper_index, zeroi, upper_max);
    __m256 lower = _mm256_i32gather_ps(index->binsv, safe_lower, 4);
    __m256 upper = _mm256_i32gather_ps(index->binsv, safe_upper, 4);
    lower = _mm256_blendv_ps(lower, _mm256_set1_ps(MINVAL), _mm256_castsi256_ps(lower_zero));
    upper = _mm256_blendv_ps(upper, _mm256_set1_ps(MAXVAL), _mm256_castsi256_ps(upper_max));
    const __m256 value = ((uintptr_t) values & 31U) == 0 ? _mm256_load_ps(values) : _mm256_loadu_ps(values);
    const __m256 below = _mm256_cmp_ps(lower, value, _CMP_GT_OS);
    const __m256 above = _mm256_cmp_ps(upper, value, _CMP_LT_OS);
    const __m256 dlow = _mm256_sub_ps(lower, value);
    const __m256 dhigh = _mm256_sub_ps(upper, value);
#if defined(__FMA__)
    __m256 squared = _mm256_and_ps(below, _mm256_fmadd_ps(dlow, dlow, _mm256_setzero_ps()));
    squared = _mm256_or_ps(squared, _mm256_and_ps(above, _mm256_fmadd_ps(dhigh, dhigh, _mm256_setzero_ps())));
#else
    __m256 squared = _mm256_or_ps(_mm256_and_ps(below, _mm256_mul_ps(dlow, dlow)),
                                   _mm256_and_ps(above, _mm256_mul_ps(dhigh, dhigh)));
#endif
    squared = _mm256_hadd_ps(squared, squared);
    squared = _mm256_hadd_ps(squared, squared);
    float sums[8];
    _mm256_storeu_ps(sums, squared);
    return sums[0] + sums[4];
}
#endif

#if defined(__AVX512F__)
static inline float messi_lower_bound_avx512_block(const isax_index *index,
                                                    const float *values,
                                                    const sax_type *sax,
                                                    const sax_type *cardinalities,
                                                    int start) {
    const int max_bits = index->settings->sax_bit_cardinality;
    const int alphabet = index->settings->sax_alphabet_cardinality;
    const int stride = alphabet - 1;
    int lower_index[16] = {0}, upper_index[16] = {0};
    __mmask16 lower_zero = 0, upper_max = 0;
    for (int lane = 0; lane < 8; ++lane) {
        const int shift = max_bits - cardinalities[start + lane];
        const int low = (sax[start + lane] >> shift) << shift;
        const int high = low | (shift == 0 ? 0 : ((1 << shift) - 1));
        if (low == 0) lower_zero |= ((__mmask16) 1 << lane);
        else lower_index[lane] = (start + lane) * stride + low - 1;
        if (high == alphabet - 1) upper_max |= ((__mmask16) 1 << lane);
        else upper_index[lane] = (start + lane) * stride + high;
    }
    const __mmask16 active = 0x00ff;
    const __m512 lower_gathered = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active,
        _mm512_loadu_si512((const void *) lower_index), index->binsv, 4);
    const __m512 upper_gathered = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active,
        _mm512_loadu_si512((const void *) upper_index), index->binsv, 4);
    const __m512 lower = _mm512_mask_mov_ps(lower_gathered, lower_zero, _mm512_set1_ps(MINVAL));
    const __m512 upper = _mm512_mask_mov_ps(upper_gathered, upper_max, _mm512_set1_ps(MAXVAL));
    const __m512 value = _mm512_maskz_loadu_ps(active, values);
    const __mmask16 below = _mm512_mask_cmp_ps_mask(active, lower, value, _CMP_GT_OS);
    const __mmask16 above = _mm512_mask_cmp_ps_mask(active, upper, value, _CMP_LT_OS);
    const __m512 dlow = _mm512_sub_ps(lower, value);
    const __m512 dhigh = _mm512_sub_ps(upper, value);
#if defined(__FMA__)
    __m512 squared = _mm512_maskz_mov_ps(below, _mm512_fmadd_ps(dlow, dlow, _mm512_setzero_ps()));
    squared = _mm512_mask_mov_ps(squared, above, _mm512_fmadd_ps(dhigh, dhigh, _mm512_setzero_ps()));
#else
    __m512 squared = _mm512_maskz_mov_ps(below, _mm512_mul_ps(dlow, dlow));
    squared = _mm512_mask_mov_ps(squared, above, _mm512_mul_ps(dhigh, dhigh));
#endif
    return _mm512_reduce_add_ps(squared);
}
#endif

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
static inline float messi_lower_bound_neon_block(const isax_index *index,
                                                 const float *values,
                                                 const sax_type *sax,
                                                 const sax_type *cardinalities,
                                                 int start) {
    const int max_bits = index->settings->sax_bit_cardinality;
    const int alphabet = index->settings->sax_alphabet_cardinality;
    const int stride = alphabet - 1;
    float lower[8], upper[8];
    for (int lane = 0; lane < 8; ++lane) {
        const int shift = max_bits - cardinalities[start + lane];
        const int low = (sax[start + lane] >> shift) << shift;
        const int high = low | (shift == 0 ? 0 : ((1 << shift) - 1));
        lower[lane] = low == 0 ? MINVAL : index->binsv[(start + lane) * stride + low - 1];
        upper[lane] = high == alphabet - 1 ? MAXVAL : index->binsv[(start + lane) * stride + high];
    }
    float32x4_t total = vdupq_n_f32(0.0f);
    for (int lane = 0; lane < 8; lane += 4) {
        const float32x4_t value = vld1q_f32(values + lane);
        const float32x4_t low = vld1q_f32(lower + lane);
        const float32x4_t high = vld1q_f32(upper + lane);
        const uint32x4_t below = vcgtq_f32(low, value);
        const uint32x4_t above = vcltq_f32(high, value);
        const float32x4_t dlow = vsubq_f32(low, value);
        const float32x4_t dhigh = vsubq_f32(high, value);
        total = vaddq_f32(total, vbslq_f32(below, vmulq_f32(dlow, dlow), vdupq_n_f32(0.0f)));
        total = vaddq_f32(total, vbslq_f32(above, vmulq_f32(dhigh, dhigh), vdupq_n_f32(0.0f)));
    }
    return vaddvq_f32(total);
}
#endif

static inline float messi_lower_bound_16(const isax_index *index,
                                         const float *values,
                                         const sax_type *sax,
                                         const sax_type *cardinalities,
                                         float bsf, float factor) {
    float first;
#if defined(__AVX512F__)
    first = messi_lower_bound_avx512_block(index, values, sax, cardinalities, 0);
#elif ADS_HAVE_AVX2
    first = messi_lower_bound_avx2_block(index, values, sax, cardinalities, 0);
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    first = messi_lower_bound_neon_block(index, values, sax, cardinalities, 0);
#else
    return messi_lower_bound_16_scalar(index, values, sax, cardinalities, bsf, factor);
#endif
    first *= factor;
    if (first > bsf) return first;
#if defined(__AVX512F__)
    return first + messi_lower_bound_avx512_block(index, values + 8, sax, cardinalities, 8) * factor;
#elif ADS_HAVE_AVX2
    return first + messi_lower_bound_avx2_block(index, values + 8, sax, cardinalities, 8) * factor;
#else
    return first + messi_lower_bound_neon_block(index, values + 8, sax, cardinalities, 8) * factor;
#endif
}

#endif
