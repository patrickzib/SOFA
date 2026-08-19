#ifndef MESSI_SFA_CALC_UTILS_H
#define MESSI_SFA_CALC_UTILS_H

#include "config.h"
#include "../../globals.h"
#include "isax_index.h"
#include "sax/sax.h"
#include "sfa/sfa.h"
#include "spartan/spartan.h"

////// Utility functions ////
// Function to calculate mean of an array of floats
float calculateMean(ts_type *data, int n);

// Function to calculate standard deviation of an array of floats
float calculateStdDev(ts_type *data, int n, ts_type mean);

// Function to perform zero mean normalization
void znorm(ts_type *data, int n);

void isax_node_mbb_sax_update(isax_node *node, const sax_type *sax, int size);
ts_type messi_minidist_range_raw(isax_index *index,
                                 float *paa_or_fft,
                                 sax_type *sax_min,
                                 sax_type *sax_max,
                                 sax_type *sax_cardinalities,
                                 float bsf);

// Shared minidist dispatch for SAX/SFA.
static inline ts_type messi_minidist_raw(isax_index *index,
                                         float *paa_or_fft,
                                         sax_type *sax,
                                         sax_type *sax_cardinalities,
                                         float bsf) {
    if (index->settings->n_segments == 16 && sizeof(sax_type) == 1) {
        if (index->settings->function_type == 4 || index->settings->function_type == 6) {
            return minidist_fft_to_sfa_rawe_SIMD(index, paa_or_fft, sax, sax_cardinalities, bsf);
        }
        if (index->settings->function_type == 5) {
            return minidist_pca_to_spartan_rawe_SIMD(index, paa_or_fft, sax, sax_cardinalities, bsf);
        }
        return minidist_paa_to_isax_raw_SIMD(paa_or_fft, sax, sax_cardinalities, index->settings);
    }

    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        return minidist_fft_to_sfa(index, paa_or_fft, sax, sax_cardinalities, bsf);
    }
    if (index->settings->function_type == 5) {
        return minidist_pca_to_spartan(index, paa_or_fft, sax, sax_cardinalities, bsf);
    }
    return minidist_paa_to_isax(paa_or_fft, sax, sax_cardinalities, index->settings, 1);
}

static inline ts_type messi_minidist(isax_index *index,
                                     float *paa_or_fft,
                                     sax_type *sax,
                                     sax_type *sax_cardinalities,
                                     float bsf) {
    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        return minidist_fft_to_sfa(index, paa_or_fft, sax, sax_cardinalities, bsf);
    }
    if (index->settings->function_type == 5) {
        return minidist_pca_to_spartan(index, paa_or_fft, sax, sax_cardinalities, bsf);
    }
    return minidist_paa_to_isax(paa_or_fft, sax, sax_cardinalities, index->settings, 0);
}

#endif //MESSI_SFA_CALC_UTILS_H
