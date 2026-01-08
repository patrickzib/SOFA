#ifndef MESSI_SFA_CALC_UTILS_H
#define MESSI_SFA_CALC_UTILS_H

#include "config.h"
#include "../../globals.h"
#include "isax_index.h"
#include "sax/sax.h"
#include "sfa/sfa.h"

////// Utility functions ////
// Function to calculate mean of an array of floats
float calculateMean(ts_type *data, int n);

// Function to calculate standard deviation of an array of floats
float calculateStdDev(ts_type *data, int n, ts_type mean);

// Function to perform zero mean normalization
void znorm(ts_type *data, int n);

// Shared minidist dispatch for SAX/SFA.
static inline ts_type messi_minidist_raw(isax_index *index,
                                         float *paa_or_fft,
                                         sax_type *sax,
                                         sax_type *sax_cardinalities,
                                         float bsf) {
    if (index->settings->SIMD_flag) {
        if (index->settings->function_type == 4) {
            return minidist_fft_to_sfa_rawe_SIMD(index, paa_or_fft, sax, sax_cardinalities, bsf);
        }
        return minidist_paa_to_isax_raw_SIMD(paa_or_fft, sax, sax_cardinalities,
                                             index->settings->sax_bit_cardinality,
                                             index->settings->sax_alphabet_cardinality,
                                             index->settings->paa_segments, MINVAL, MAXVAL,
                                             index->settings->mindist_sqrt);
    }

    if (index->settings->function_type == 4) {
        return minidist_fft_to_sfa_raw(index, paa_or_fft, sax, sax_cardinalities, bsf);
    }
    return minidist_paa_to_isax_raw(paa_or_fft, sax, sax_cardinalities,
                                    index->settings->sax_bit_cardinality,
                                    index->settings->sax_alphabet_cardinality,
                                    index->settings->paa_segments, MINVAL, MAXVAL,
                                    index->settings->mindist_sqrt);
}

static inline ts_type messi_minidist(isax_index *index,
                                     float *paa_or_fft,
                                     sax_type *sax,
                                     sax_type *sax_cardinalities,
                                     float bsf) {
    if (index->settings->function_type == 4) {
        return minidist_fft_to_sfa(index, paa_or_fft, sax, sax_cardinalities, bsf);
    }
    return minidist_paa_to_isax(paa_or_fft, sax, sax_cardinalities,
                                index->settings->sax_bit_cardinality,
                                index->settings->sax_alphabet_cardinality,
                                index->settings->paa_segments, MINVAL, MAXVAL,
                                index->settings->mindist_sqrt);
}

#endif //MESSI_SFA_CALC_UTILS_H
