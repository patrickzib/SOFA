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
                           sax_type max_bit_cardinality,
                           int max_cardinality,
                           int number_of_segments, 
                           int min_val, int max_val, float ratio_sqrt);

float minidist_paa_to_isax_raw(float *paa, sax_type *sax, sax_type *sax_cardinalities,
						                              sax_type max_bit_cardinality,
						                              int max_cardinality,
						                              int number_of_segments, 
						                              int min_val, int max_val, float ratio_sqrt);
float   minidist_paa_to_isax_raw_SIMD(float *paa, sax_type *sax, 
                           sax_type *sax_cardinalities,
                           sax_type max_bit_cardinality,
                           int max_cardinality,
                           int number_of_segments,
                           int min_val,
                           int max_val,
                           float ratio_sqrt);
float   minidist_paa_to_isax_rawa_SIMD(float *paa, sax_type *sax, 
                           sax_type *sax_cardinalities,
                           sax_type max_bit_cardinality,
                           int max_cardinality,
                           int number_of_segments,
                           int min_val,
                           int max_val,
                           float ratio_sqrt);
#if ADS_HAVE_AVX2
#else
#define minidist_paa_to_isax_raw_SIMD minidist_paa_to_isax
#define minidist_paa_to_isax_rawa_SIMD minidist_paa_to_isax
#endif

enum response paa_from_ts (ts_type *ts_in, ts_type *paa_out, isax_index_settings *settings);
enum response sax_from_paa (ts_type *paa, sax_type *sax, isax_index_settings *settings);
//float* bsearchsimdpaa(float paaseg, float* sax_breakpoints, int candinalityoffset);
#endif
