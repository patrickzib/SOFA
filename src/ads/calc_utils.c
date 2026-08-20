#include "config.h"
#include "../../globals.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "math.h"

#include "ads/calc_utils.h"
#include "ads/sax/sax_breakpoints.h"


////// Utility functions ////
// Function to calculate mean of an array of floats
float calculateMean(ts_type *data, int n) {
    float sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += data[i];
    }
    return sum / n;
}

// Function to calculate standard deviation of an array of floats
float calculateStdDev(ts_type *data, int n, ts_type mean) {
    ts_type sumSquaredDiffs = 0.0;
    for (int i = 0; i < n; i++) {
        sumSquaredDiffs += (data[i] - mean) * (data[i] - mean);
    }
    return sqrt(sumSquaredDiffs / n);
}

// Function to perform zero mean normalization
void znorm(ts_type *data, int n) {
    // printf("No nor applied");
    // Calculate mean
    ts_type mean = calculateMean(data, n);

    // Calculate standard deviation
    ts_type stdDev = calculateStdDev(data, n, mean);

    // Normalize each data point
    if (stdDev < 1e-8) {
        stdDev = 1.0;
    }

    for (int i = 0; i < n; i++) {
        data[i] = (data[i] - mean) / stdDev;
    }
}

root_mask_type isax_root_mask_from_sax(const isax_index *index,
                                       const sax_type *sax,
                                       int uniform_kn) {
    const isax_index_settings *settings = index->settings;
    root_mask_type mask = 0;

    if (settings->root_bit_cardinalities != NULL) {
        int bit = settings->n_segments - 1;
        for (int i = 0; i < settings->n_segments && bit >= 0; ++i) {
            int count = settings->root_bit_cardinalities[i];
            if (count > settings->sax_bit_cardinality) {
                count = settings->sax_bit_cardinality;
            }
            for (int j = 0; j < count && bit >= 0; ++j, --bit) {
                if (sax[i] & settings->bit_masks[settings->sax_bit_cardinality - 1 - j]) {
                    mask |= ((root_mask_type) 1 << bit);
                }
            }
        }
        return mask;
    }

    if (uniform_kn < 1) {
        uniform_kn = 1;
    }
    for (int i = 0; i < settings->n_segments / uniform_kn; ++i) {
        for (int j = 0; j < uniform_kn; ++j) {
            if (sax[i] & settings->bit_masks[settings->sax_bit_cardinality - 1 - j]) {
                mask |= ((root_mask_type) 1 << (settings->n_segments - i * uniform_kn - j - 1));
            }
        }
    }
    return mask;
}

enum response isax_configure_variance_root_split(isax_index *index,
                                                 ts_type *const *coefficients,
                                                 unsigned int sample_size) {
    isax_index_settings *settings = index->settings;
    if (!settings->dynamic_root_split_variance) {
        return SUCCESS;
    }
    if (coefficients == NULL || sample_size == 0) {
        return FAILURE;
    }

    const int dimensions = settings->n_segments;
    const int budget = dimensions < (int) (sizeof(root_mask_type) * 8)
                           ? dimensions
                           : (int) (sizeof(root_mask_type) * 8);
    sax_type *bits = calloc((size_t) dimensions, sizeof(*bits));
    double *variance = calloc((size_t) dimensions, sizeof(*variance));
    if (bits == NULL || variance == NULL) {
        free(bits);
        free(variance);
        return FAILURE;
    }

    for (int i = 0; i < dimensions; ++i) {
        double mean = 0.0;
        for (unsigned int j = 0; j < sample_size; ++j) {
            mean += coefficients[i][j];
        }
        mean /= sample_size;
        for (unsigned int j = 0; j < sample_size; ++j) {
            double delta = coefficients[i][j] - mean;
            variance[i] += delta * delta;
        }
        variance[i] /= sample_size;
        if (!isfinite(variance[i]) || variance[i] < 0.0) {
            variance[i] = 0.0;
        }
    }

    for (int assigned = 0; assigned < budget; ++assigned) {
        int best = -1;
        double best_gain = -1.0;
        for (int i = 0; i < dimensions; ++i) {
            if (bits[i] >= settings->sax_bit_cardinality) {
                continue;
            }
            double gain = ldexp(variance[i], -(int) bits[i] - 1);
            if (gain > best_gain) {
                best_gain = gain;
                best = i;
            }
        }
        if (best < 0) {
            free(bits);
            free(variance);
            return FAILURE;
        }
        ++bits[best];
    }

    free(settings->root_bit_cardinalities);
    settings->root_bit_cardinalities = bits;
    free(variance);

    fprintf(stderr, ">>> variance root split (%d bits):", budget);
    for (int i = 0; i < dimensions; ++i) {
        fprintf(stderr, " %u", (unsigned int) bits[i]);
    }
    fprintf(stderr, "\n");
    return SUCCESS;
}

void isax_node_mbb_sax_update(isax_node *node, const sax_type *sax, int size) {
    if (node == NULL || sax == NULL || size <= 0) {
        return;
    }
    if (node->mbb_sax_min == NULL) {
        node->mbb_sax_min = malloc(sizeof(sax_type) * (size_t) size);
    }
    if (node->mbb_sax_max == NULL) {
        node->mbb_sax_max = malloc(sizeof(sax_type) * (size_t) size);
    }
    if (node->mbb_sax_min == NULL || node->mbb_sax_max == NULL) {
        fprintf(stderr, "error: could not allocate SAX MBB arrays.\n");
        return;
    }
    if (!node->mbb_sax_valid) {
        for (int i = 0; i < size; ++i) {
            node->mbb_sax_min[i] = sax[i];
            node->mbb_sax_max[i] = sax[i];
        }
        node->mbb_sax_valid = 1;
        return;
    }
    for (int i = 0; i < size; ++i) {
        sax_type value = sax[i];
        if (value < node->mbb_sax_min[i]) {
            node->mbb_sax_min[i] = value;
        }
        if (value > node->mbb_sax_max[i]) {
            node->mbb_sax_max[i] = value;
        }
    }
}

static ts_type messi_minidist_range_bins(const ts_type *bins,
                                         float value,
                                         sax_type sax_min,
                                         sax_type sax_max,
                                         int max_cardinality,
                                         float factor) {
    sax_type region_lower = sax_min;
    sax_type region_upper = sax_max;

    float breakpoint_lower = (region_lower == 0) ? MINVAL : bins[region_lower - 1];
    float breakpoint_upper = (region_upper == max_cardinality - 1) ? MAXVAL : bins[region_upper];

    if (breakpoint_lower > value) {
        ts_type diff = breakpoint_lower - value;
        return factor * diff * diff;
    }
    if (breakpoint_upper < value) {
        ts_type diff = value - breakpoint_upper;
        return factor * diff * diff;
    }
    return 0.0f;
}

ts_type messi_minidist_range_raw(isax_index *index,
                                 float *paa_or_fft,
                                 sax_type *sax_min,
                                 sax_type *sax_max,
                                 sax_type *sax_cardinalities,
                                 float bsf) {
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->n_segments;

    ts_type distance = 0.0f;

    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        int i = 0;
        if (!index->settings->is_norm &&
            (index->settings->n_coefficients == 0 || index->coefficients[0] == 0)) {
            distance += messi_minidist_range_bins(index->bins[0], paa_or_fft[0],
                                                  sax_min[0], sax_max[0],
                                                  max_cardinality, 1.0f);
            if (distance > bsf) {
                return distance;
            }
            i = (index->settings->n_coefficients == 0) ? 2 : 1;
        }
        for (; i < number_of_segments; ++i) {
            distance += messi_minidist_range_bins(index->bins[i], paa_or_fft[i],
                                                  sax_min[i], sax_max[i],
                                                  max_cardinality, 2.0f);
            if (distance > bsf) {
                return distance;
            }
        }
        return distance;
    }

    if (index->settings->function_type == 5) {
        for (int i = 0; i < number_of_segments; i++) {
            distance += messi_minidist_range_bins(index->bins[i], paa_or_fft[i],
                                                  sax_min[i], sax_max[i],
                                                  max_cardinality, 1.0f);
            if (distance > bsf) {
                return distance;
            }
        }
        return distance;
    }

    int offset = ((max_cardinality - 1) * (max_cardinality - 2)) / 2;
    float ratio_sqrt = index->settings->mindist_sqrt;

    for (int i = 0; i < number_of_segments; i++) {
        sax_type c_c = sax_cardinalities[i];
        sax_type c_m = max_bit_cardinality;

        int shift = c_m - c_c;
        sax_type region_lower = (sax_min[i] >> shift) << shift;
        unsigned int low_mask = (shift == 0) ? 0U : ((1U << shift) - 1U);
        sax_type region_upper = (sax_type)(((sax_max[i] >> shift) << shift) | low_mask);

        float breakpoint_lower = (region_lower == 0)
                ? MINVAL
                : sax_breakpoints[offset + region_lower - 1];
        float breakpoint_upper = (region_upper == max_cardinality - 1)
                ? MAXVAL
                : sax_breakpoints[offset + region_upper];

        if (breakpoint_lower > paa_or_fft[i]) {
            float diff = breakpoint_lower - paa_or_fft[i];
            distance += diff * diff;
        } else if (breakpoint_upper < paa_or_fft[i]) {
            float diff = paa_or_fft[i] - breakpoint_upper;
            distance += diff * diff;
        }
        // if (distance > bsf) {
        //     return ratio_sqrt * distance;
        // }
    }

    return ratio_sqrt * distance;
}

static ts_type get_lb_distance(const ts_type *bins, float value, sax_type v, sax_type c_c,
                               sax_type c_m, int max_cardinality, float factor, int use_raw) {
    int shift = c_m - c_c;
    sax_type region_lower = use_raw ? (sax_type)((v >> shift) << shift) : (sax_type)(v << shift);
    sax_type region_upper = region_lower | ((shift == 0) ? 0U : ((1U << shift) - 1U));

    float breakpoint_lower = (region_lower == 0) ? MINVAL : bins[region_lower - 1];
    float breakpoint_upper = (region_upper == max_cardinality - 1) ? MAXVAL : bins[region_upper];

    if (breakpoint_lower > value) {
        ts_type diff = breakpoint_lower - value;
        return factor * diff * diff;
    }
    if (breakpoint_upper < value) {
        ts_type diff = value - breakpoint_upper;
        return factor * diff * diff;
    }
    return 0.0f;
}

float minidist_paa_to_isax(float *paa, sax_type *sax,
                                  sax_type *sax_cardinalities,
                                  const isax_index_settings *settings,
                                  int use_raw) {
    sax_type max_bit_cardinality = settings->sax_bit_cardinality;
    int max_cardinality = settings->sax_alphabet_cardinality;
    int number_of_segments = settings->n_segments;
    float ratio_sqrt = settings->mindist_sqrt;

    float distance = 0.0f;
    int offset = ((max_cardinality - 1) * (max_cardinality - 2)) / 2;
    const ts_type *bins = sax_breakpoints + offset;

    for (int i = 0; i < number_of_segments; i++) {
        sax_type c_c = sax_cardinalities[i];
        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        distance += get_lb_distance(bins, paa[i], v, c_c, c_m, max_cardinality, ratio_sqrt, use_raw);
    }

    return distance;
}

ts_type minidist_fft_to_sfa(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf) {
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->n_segments;

    ts_type distance = 0.0f;
    int start_index = 0;

    if (!index->settings->is_norm &&
        (index->settings->n_coefficients == 0 || index->coefficients[0] == 0)) {
        sax_type c_c = sax_cardinalities[0];
        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[0];
        distance += get_lb_distance(index->bins[0], fft[0], v, c_c, c_m, max_cardinality, 1.0f, 0);
        if (distance > bsf) {
            return distance;
        }
        start_index = (index->settings->n_coefficients == 0) ? 2 : 1;
    }

    for (int i = start_index; i < number_of_segments; i++) {
        sax_type c_c = sax_cardinalities[i];
        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        distance += get_lb_distance(index->bins[i], fft[i], v, c_c, c_m, max_cardinality, 2.0f, 0);
        if (distance > bsf) {
            return distance;
        }
    }
    return distance;
}

ts_type minidist_pca_to_spartan(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities, float bsf) {
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->n_segments;

    ts_type distance = 0.0f;
    for (int i = 0; i < number_of_segments; i++) {
        sax_type c_c = sax_cardinalities[i];
        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        distance += get_lb_distance(index->bins[i], pca[i], v, c_c, c_m, max_cardinality, 1.0f, 0);
        if (distance > bsf) {
            return distance;
        }
    }
    return distance;
}
