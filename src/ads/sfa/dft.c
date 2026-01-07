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
        isax_index *index,
        int coeff_number, int best_only,
        fftw_workspace *fftw) {
    fftwf_execute(fftw->plan_forward);

    // Image part of first (DC) coefficient
    fftw->ts_out[0][1] = 0;

    int j = 0;

    // if normalized, ignore first coeff and start with offset 1
    int start_offset = index->settings->is_norm ? 1 : 0;

    if (best_only) {
        for (int k = 0; k < coeff_number / 2; ++k, j+= 2) {
            int coeff = index->coefficients[k] + start_offset;
            fftw->transform[j] = fftw->ts_out[coeff][0];
            fftw->transform[j + 1] = fftw->ts_out[coeff][1] * -1;
        }
    } else {
        for (int k = start_offset; k < coeff_number / 2 + start_offset; ++k, j+= 2) {
            fftw->transform[j] = fftw->ts_out[k][0];
            fftw->transform[j + 1] = fftw->ts_out[k][1] * -1;
        }
    }

    // normalizing fft result in frequency domain to allow for lower bounding
    ts_type norm_factor = index->norm_factor;
    for (int i = 0; i < coeff_number; ++i) {
        fftw->transform[i] *= norm_factor;
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
    int cardinality = index->settings->sax_alphabet_cardinality;
    int offset = ((cardinality - 1) * (cardinality - 2)) / 2;

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
enum response sfa_from_ts(isax_index *index, sax_type *sax_out, fftw_workspace *fftw) {

    int use_best = index->settings->coeff_number != 0;
    fft_from_ts(index, index->settings->paa_segments, use_best, fftw);

    ts_type *cur_coeff_line = calloc(index->settings->paa_segments, sizeof(ts_type));
    memcpy(cur_coeff_line, fftw->transform, sizeof(ts_type) * index->settings->paa_segments);

    sfa_from_fft(index, cur_coeff_line, sax_out);

    free(cur_coeff_line);

    if (sax_out != NULL) return SUCCESS;
    else {
        fprintf(stderr, "SFA error");
    }
}
