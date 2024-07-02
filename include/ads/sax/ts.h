//
//  ts.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//

#ifndef isaxlib_ts_h
#define isaxlib_ts_h
#include "../../../config.h"
#include "../../../globals.h"
void ts_parse_str(char ts_str[], ts_type *ts_out, int ts_size, const char * delims);
void ts_print(ts_type *ts, int size);

float ts_euclidean_distance_dot_product(ts_type *X, ts_type *Y, int size);
float ts_euclidean_distance_dot_product_SIMD(float *X, float *Y, int size);

float ts_euclidean_distance(ts_type * t, ts_type * s, int size, float bound);
float ts_euclidean_distance_SIMD(ts_type * t, ts_type * s, int size, float bound);

float ts_ed(ts_type * t, ts_type * s, int size, float bound, char is_simd, char is_norm);
#endif
