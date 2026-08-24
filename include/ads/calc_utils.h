#ifndef MESSI_SFA_CALC_UTILS_H
#define MESSI_SFA_CALC_UTILS_H

#include "config.h"
#include "../../globals.h"
#include "isax_index.h"

/* Trie record bounds may use a prefix of up to 64 transform dimensions. */
#define MESSI_RECORD_LB_MAX_DIMENSIONS 64

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
/* Return the full node-MBR lower bound and, when requested, the contribution
 * from dimensions starting at unselected_start_dimension. */
ts_type messi_minidist_range_raw_partitioned(isax_index *index,
                                             float *paa_or_fft,
                                             sax_type *sax_min,
                                             sax_type *sax_max,
                                             sax_type *sax_cardinalities,
                                             float bsf,
                                             int unselected_start_dimension,
                                             float *unselected_distance);
int messi_build_record_lb_table(const isax_index *index,
                                const float *paa_or_fft,
                                int dimensions,
                                float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256]);

// Shared minidist dispatch. Representation-specific dependencies stay private
// to calc_utils.c so users of this header do not import SFA or SPARTAN.
ts_type messi_minidist_raw(isax_index *index, float *paa_or_fft, sax_type *sax,
                           sax_type *sax_cardinalities, float bsf);
ts_type messi_minidist(isax_index *index, float *paa_or_fft, sax_type *sax,
                       sax_type *sax_cardinalities, float bsf);

#endif //MESSI_SFA_CALC_UTILS_H
