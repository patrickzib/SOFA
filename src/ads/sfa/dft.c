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
#include "config.h"
#include "globals.h"
#include "ads/isax_index.h"
#include "math.h"

#include <fftw3.h>

#include "ads/sfa/dft.h"

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
        for (int k = 0; k < coeff_number / 2 + start_offset; ++k) {
            int coeff = index->coefficients[k] + start_offset;
            transform[j] = ts_out[coeff][0];
            transform[j + 1] = ts_out[coeff][1];
            j += 2;
        }
    } else {
        for (int k = start_offset; k < coeff_number / 2 + start_offset; ++k) {
            transform[j] = ts_out[k][0];
            transform[j + 1] = ts_out[k][1];
            j += 2;
        }
    }

    // normalizing fft result in frequency domain
    int sign = 1;
    ts_type norm_factor = index->norm_factor;

    for (int i = 0; i < coeff_number; ++i) {
        transform[i] *= norm_factor * sign;
        sign *= -1;
    }
    return;
}


/*
    This function discretized FFT coefficients with the intervals from MCB
    The current transform is pointed to by dft_mem_array
*/
void sfa_from_fft(isax_index *index, ts_type *cur_transform, unsigned char *cur_sfa_word) {
    unsigned long ts_length = index->settings->timeseries_size;
    int paa_segments = index->settings->paa_segments;

    for (int k = 0; k < paa_segments; ++k) {
        unsigned int c;
        for (c = 0; c < index->settings->sax_alphabet_cardinality - 1; c++) {
            if (cur_transform[k] < index->bins[k][c]) {
                break;
            }
        }
        cur_sfa_word[k] = (unsigned char) (c);
    }
}

/*
    This function creates an SFA representation of a time series 
*/
enum response sfa_from_ts(isax_index *index, ts_type *ts_in, sax_type *sax_out, fftwf_complex *ts_out, ts_type *transform,
            fftwf_plan plan_forward) {

    int use_best = index->settings->coeff_number != 0;
    fft_from_ts(index, ts_in, index->settings->paa_segments, use_best, ts_out, transform, plan_forward);

    ts_type *cur_coeff_line = calloc(index->settings->paa_segments, sizeof(ts_type));

    for (int i = 0; i < index->settings->paa_segments; ++i) {
        // cur_coeff_line[i] = (ts_type) roundf(transform[i] * 100.0) / 100.0;
        cur_coeff_line[i] = transform[i];
    }

    sfa_from_fft(index, cur_coeff_line, sax_out);

    free(cur_coeff_line);

    if (sax_out != NULL) return SUCCESS;
    else {
        fprintf(stderr, "SFA error");
    }
}
