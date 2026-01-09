//
//  sax.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//

#ifndef isaxlib_sax_h
#define isaxlib_sax_h
#include "config.h"
#include "../../../globals.h"
#include "ts.h"
#include "../isax_index.h"

enum response sax_from_ts(ts_type *ts_in, sax_type *sax_out, isax_index_settings *settings);
enum response sax_from_ts_new(ts_type *ts_in, sax_type *sax_out, isax_index_settings *settings);
void sax_print(sax_type *sax, int segments, int cardinality);
void printbin(unsigned long long n, int size);
void serial_printbin (unsigned long long n, int size);
int compare(const void *a, const void *b);
float minidist_paa_to_isax(float *paa, sax_type *sax, sax_type *sax_cardinalities,
                           const isax_index_settings *settings);

float minidist_paa_to_isax_raw(float *paa, sax_type *sax, sax_type *sax_cardinalities,
                               const isax_index_settings *settings);
#if ADS_HAVE_AVX2
float   minidist_paa_to_isax_raw_SIMD(float *paa, sax_type *sax,
                           sax_type *sax_cardinalities,
                           const isax_index_settings *settings);
#else
static inline float minidist_paa_to_isax_raw_SIMD(float *paa, sax_type *sax,
                           sax_type *sax_cardinalities,
                           const isax_index_settings *settings) {
    return minidist_paa_to_isax_raw(paa, sax, sax_cardinalities, settings);
}
#endif
#if ADS_HAVE_AVX2
float   minidist_paa_to_isax_rawa_SIMD(float *paa, sax_type *sax,
                           sax_type *sax_cardinalities,
                           const isax_index_settings *settings);
#else
static inline float minidist_paa_to_isax_rawa_SIMD(float *paa, sax_type *sax,
                           sax_type *sax_cardinalities,
                           const isax_index_settings *settings) {
    return minidist_paa_to_isax(paa, sax, sax_cardinalities, settings);
}
#endif

enum response paa_from_ts (ts_type *ts_in, ts_type *paa_out, isax_index_settings *settings);
enum response sax_from_paa (ts_type *paa, sax_type *sax, isax_index_settings *settings);
//float* bsearchsimdpaa(float paaseg, float* sax_breakpoints, int candinalityoffset);
#endif
