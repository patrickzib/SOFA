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
enum response fft_from_ts(
        isax_index *index,
        int n_coefficients, int best_only,
        fftw_workspace *fftw) {
    if (index == NULL || index->settings == NULL || fftw == NULL || fftw->plan_forward == NULL ||
        fftw->ts_out == NULL || fftw->transform == NULL || n_coefficients <= 0 ||
        n_coefficients % 2 != 0 || n_coefficients > index->settings->timeseries_size ||
        (best_only && index->coefficients == NULL)) {
        return FAILURE;
    }

    fftwf_execute(fftw->plan_forward);

    // Image part of first (DC) coefficient
    fftw->ts_out[0][1] = 0;

    int j = 0;

    // if normalized, ignore first coeff and start with offset 1
    int start_offset = index->settings->is_norm ? 1 : 0;

    if (best_only) {
        for (int k = 0; k < n_coefficients / 2; ++k, j+= 2) {
            int coeff = index->coefficients[k] + start_offset;
            fftw->transform[j] = fftw->ts_out[coeff][0];
            fftw->transform[j + 1] = fftw->ts_out[coeff][1] * -1;
        }
    } else {
        for (int k = start_offset; k < n_coefficients / 2 + start_offset; ++k, j+= 2) {
            fftw->transform[j] = fftw->ts_out[k][0];
            fftw->transform[j + 1] = fftw->ts_out[k][1] * -1;
        }
    }

    // normalizing fft result in frequency domain to allow for lower bounding
    ts_type norm_factor = index->norm_factor;
    for (int i = 0; i < n_coefficients; ++i) {
        fftw->transform[i] *= norm_factor;
    }
    return SUCCESS;
}

/* Return the first bin boundary strictly above value, preserving the legacy
 * linear scan's behavior for values equal to a boundary. */
static unsigned int sfa_symbol_from_value(const ts_type *boundaries, int count,
                                          ts_type value) {
    int low = 0, high = count;
    while (low < high) {
        const int middle = low + (high - low) / 2;
        if (value < boundaries[middle]) high = middle;
        else low = middle + 1;
    }
    return (unsigned int) low;
}

/*
    This function discretized FFT coefficients with the intervals from MCB
    The current transform is pointed to by dft_mem_array
*/
void sfa_from_fft(isax_index *index, const ts_type *cur_transform, sax_type *cur_sfa_word) {
    const int n_segments = index->settings->n_segments;
    const int boundary_count = index->settings->sax_alphabet_cardinality - 1;
    for (int k = 0; k < n_segments; ++k)
        cur_sfa_word[k] = (sax_type) sfa_symbol_from_value(index->bins[k], boundary_count,
                                                            cur_transform[k]);
}

/*
    This function creates an SFA representation of a time series 
*/
enum response sfa_from_ts(isax_index *index, const ts_type *ts, sax_type *sax_out,
                          fftw_workspace *fftw) {

    int use_best = index->settings->n_coefficients != 0;
    memcpy(fftw->ts, ts, sizeof(ts_type) * index->settings->timeseries_size);
    if (fft_from_ts(index, index->settings->n_segments, use_best, fftw) != SUCCESS) {
        return FAILURE;
    }

    sfa_from_fft(index, fftw->transform, sax_out);

    if (sax_out != NULL) return SUCCESS;
    else {
        fprintf(stderr, "SFA error");
    }
    return FAILURE;
}
