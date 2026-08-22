#ifndef MESSI_SFA_CALC_UTILS_H
#define MESSI_SFA_CALC_UTILS_H

#include "config.h"
#include "../../globals.h"
#include "isax_index.h"
#include "sax/sax.h"
#include "sfa/sfa.h"
#include "spartan/spartan.h"

/* Monotonic wall-clock time for concise construction/query phase reporting. */
double messi_monotonic_seconds(void);

////// Utility functions ////
// Function to calculate mean of an array of floats
float calculateMean(ts_type *data, int n);

// Function to calculate standard deviation of an array of floats
float calculateStdDev(ts_type *data, int n, ts_type mean);

// Function to perform zero mean normalization
void znorm(ts_type *data, int n);

void isax_node_mbb_sax_update(isax_node *node, const sax_type *sax, int size);
enum response configure_dynamic_bit_allocation(isax_index *index,
                                                 const double *variance,
                                                 int dimensions,
                                                 int budget,
                                                 int min_bit_val,
                                                 int max_bit_val);
root_mask_type isax_root_mask_from_sax(const isax_index *index,
                                       const sax_type *sax,
                                       int uniform_kn);
ts_type messi_minidist_range_raw(isax_index *index,
                                 float *paa_or_fft,
                                 sax_type *sax_min,
                                 sax_type *sax_max,
                                 sax_type *sax_cardinalities,
                                 float bsf);
int messi_build_record_lb_table(const isax_index *index,
                                const float *paa_or_fft,
                                int dimensions,
                                float table[16][256]);

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
