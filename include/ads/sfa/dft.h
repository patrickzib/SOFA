#ifndef sfalib_dft_h
#define sfalib_dft_h

#include "../../../config.h"
#include "../../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <fftw3.h>
#include <sys/types.h>

void
fft_from_ts(isax_index *index, ts_type *ts, int coeff_number, int best_only, fftwf_complex *ts_out, ts_type *transform,
            fftwf_plan plan_forward);

// void fft_from_ts_coeff(isax_index *index, ts_type *ts,  fftwf_complex *ts_out, ts_type *transform, fftwf_plan plan_forward);
// void fft_from_ts_all_coeff(isax_index *index, ts_type *ts,  fftwf_complex *ts_out, ts_type *transform, fftwf_plan plan_forward);
void sfa_from_fft(isax_index *index, ts_type *cur_transform, unsigned char *cur_sfa_word);

enum response
sfa_from_ts(isax_index *index, ts_type *ts_in, sax_type *sax_out, fftwf_complex *ts_out, ts_type *transform,
            fftwf_plan plan_forward);

#endif
