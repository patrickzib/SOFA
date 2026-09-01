#include "config.h"
#include "../../globals.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "math.h"

#include "ads/calc_utils.h"
#include "ads/lower_bound_simd.h"
#include "ads/sax/sax.h"
#include "ads/sfa/sfa.h"
#include "ads/spartan/spartan.h"
#include "ads/sax/sax_breakpoints.h"

double messi_monotonic_seconds(void) {
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    return (double) time.tv_sec + (double) time.tv_nsec / 1000000000.0;
}

ts_type messi_minidist_raw(isax_index *index, float *paa_or_fft, sax_type *sax,
                           sax_type *sax_cardinalities, float bsf) {
    if (index->settings->n_segments == 16 && sizeof(sax_type) == 1) {
        if (index->settings->function_type == 4)
            return messi_lower_bound_16(index, paa_or_fft, sax, sax_cardinalities, bsf, 2.0f);
        if (index->settings->function_type == 6)
            return messi_lower_bound_16(index, paa_or_fft, sax, sax_cardinalities, bsf, 1.0f);
        if (index->settings->function_type == 7)
            return messi_lower_bound_16(index, paa_or_fft, sax, sax_cardinalities, bsf, 1.0f);
        if (index->settings->function_type == 5)
            return messi_lower_bound_16(index, paa_or_fft, sax, sax_cardinalities, bsf, 1.0f);
        return minidist_paa_to_isax_raw_SIMD(paa_or_fft, sax, sax_cardinalities, index->settings);
    }
    if (index->settings->function_type == 4 || index->settings->function_type == 6 || index->settings->function_type == 7)
        return minidist_fft_to_sfa(index, paa_or_fft, sax, sax_cardinalities, bsf);
    if (index->settings->function_type == 5)
        return minidist_pca_to_spartan(index, paa_or_fft, sax, sax_cardinalities, bsf);
    return minidist_paa_to_isax(paa_or_fft, sax, sax_cardinalities, index->settings, 1);
}

ts_type messi_minidist(isax_index *index, float *paa_or_fft, sax_type *sax,
                       sax_type *sax_cardinalities, float bsf) {
    if (index->settings->function_type == 4 || index->settings->function_type == 6 || index->settings->function_type == 7)
        return minidist_fft_to_sfa(index, paa_or_fft, sax, sax_cardinalities, bsf);
    if (index->settings->function_type == 5)
        return minidist_pca_to_spartan(index, paa_or_fft, sax, sax_cardinalities, bsf);
    return minidist_paa_to_isax(paa_or_fft, sax, sax_cardinalities, index->settings, 0);
}


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

enum response configure_dynamic_bit_allocation(isax_index *index,
                                                 const double *variance,
                                                 int dimensions,
                                                 int budget,
                                                 int min_bit_val,
                                                 int max_bit_val) {
    isax_index_settings *settings = index->settings;
    if (!settings->dynamic_root_split_variance && settings->index_type != MESSI_INDEX_TRIE) {
        return SUCCESS;
    }
    if (variance == NULL || dimensions <= 0 || budget < 0 || min_bit_val < 0 || max_bit_val <= 0 ||
        min_bit_val > max_bit_val) {
        return FAILURE;
    }
    if (max_bit_val > 8) max_bit_val = 8;
    if (min_bit_val > max_bit_val || budget < min_bit_val * dimensions) return FAILURE;

    const int allocate_bits = settings->dynamic_root_split_variance ||
                               (settings->index_type == MESSI_INDEX_TRIE &&
                                settings->trie_dynamic_alphabet);
    sax_type *bits = allocate_bits
                         ? calloc((size_t) dimensions, sizeof(*bits)) : NULL;
    double *stored_variance = calloc((size_t) dimensions, sizeof(*stored_variance));
    if ((allocate_bits && bits == NULL) || stored_variance == NULL) {
        free(bits);
        free(stored_variance);
        return FAILURE;
    }
    memcpy(stored_variance, variance, sizeof(*stored_variance) * (size_t) dimensions);
    free(settings->symbolic_variances);
    settings->symbolic_variances = stored_variance;

    if (!allocate_bits) return SUCCESS;

    for (int i = 0; i < dimensions; ++i) bits[i] = (sax_type) min_bit_val;

    for (int assigned = min_bit_val * dimensions; assigned < budget; ++assigned) {
        int best = -1;
        double best_gain = -1.0;
        for (int i = 0; i < dimensions; ++i) {
            if (bits[i] >= max_bit_val) {
                continue;
            }
            double gain = ldexp(stored_variance[i], -(int) bits[i] - 1);
            if (gain > best_gain) {
                best_gain = gain;
                best = i;
            }
        }
        if (best < 0) {
            free(bits);
            free(settings->symbolic_variances);
            settings->symbolic_variances = NULL;
            return FAILURE;
        }
        ++bits[best];
    }

    free(settings->root_bit_cardinalities);
    settings->root_bit_cardinalities = bits;
    fprintf(stderr, ">>> variance root split (%d bits, min %d, max %d):",
            budget, min_bit_val, max_bit_val);
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

ts_type messi_minidist_range_raw_partitioned(isax_index *index,
                                             float *paa_or_fft,
                                             sax_type *sax_min,
                                             sax_type *sax_max,
                                             sax_type *sax_cardinalities,
                                             float bsf,
                                             int unselected_start_dimension,
                                             float *unselected_distance) {
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->n_segments;

    ts_type distance = 0.0f;
    ts_type unselected = 0.0f;

    if (index->settings->function_type == 6 || index->settings->function_type == 7) {
        for (int i = 0; i < number_of_segments; ++i) {
            const ts_type contribution = messi_minidist_range_bins(index->bins[i], paa_or_fft[i],
                                                                    sax_min[i], sax_max[i],
                                                                    max_cardinality, 1.0f);
            distance += contribution;
            if (i >= unselected_start_dimension) unselected += contribution;
            if (distance > bsf) {
                if (unselected_distance != NULL) *unselected_distance = unselected;
                return distance;
            }
        }
        if (unselected_distance != NULL) *unselected_distance = unselected;
        return distance;
    }

    if (index->settings->function_type == 4) {
        int i = 0;
        if (!index->settings->is_norm &&
            (index->settings->n_coefficients == 0 || index->coefficients[0] == 0)) {
            const ts_type contribution = messi_minidist_range_bins(index->bins[0], paa_or_fft[0],
                                                                    sax_min[0], sax_max[0],
                                                                    max_cardinality, 1.0f);
            distance += contribution;
            if (i >= unselected_start_dimension) unselected += contribution;
            if (distance > bsf) {
                if (unselected_distance != NULL) *unselected_distance = unselected;
                return distance;
            }
            i = (index->settings->n_coefficients == 0) ? 2 : 1;
        }
        for (; i < number_of_segments; ++i) {
            const ts_type contribution = messi_minidist_range_bins(index->bins[i], paa_or_fft[i],
                                                                    sax_min[i], sax_max[i],
                                                                    max_cardinality, 2.0f);
            distance += contribution;
            if (i >= unselected_start_dimension) unselected += contribution;
            if (distance > bsf) {
                if (unselected_distance != NULL) *unselected_distance = unselected;
                return distance;
            }
        }
        if (unselected_distance != NULL) *unselected_distance = unselected;
        return distance;
    }

    if (index->settings->function_type == 5) {
        for (int i = 0; i < number_of_segments; i++) {
            const ts_type contribution = messi_minidist_range_bins(index->bins[i], paa_or_fft[i],
                                                                    sax_min[i], sax_max[i],
                                                                    max_cardinality, 1.0f);
            distance += contribution;
            if (i >= unselected_start_dimension) unselected += contribution;
            if (distance > bsf) {
                if (unselected_distance != NULL) *unselected_distance = unselected;
                return distance;
            }
        }
        if (unselected_distance != NULL) *unselected_distance = unselected;
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

        ts_type contribution = 0.0f;
        if (breakpoint_lower > paa_or_fft[i]) {
            float diff = breakpoint_lower - paa_or_fft[i];
            contribution = diff * diff;
        } else if (breakpoint_upper < paa_or_fft[i]) {
            float diff = paa_or_fft[i] - breakpoint_upper;
            contribution = diff * diff;
        }
        distance += contribution;
        if (i >= unselected_start_dimension) unselected += contribution;
        // if (ratio_sqrt * distance > bsf) {
        //     return ratio_sqrt * distance;
        // }
    }

    if (unselected_distance != NULL) *unselected_distance = ratio_sqrt * unselected;
    return ratio_sqrt * distance;
}

ts_type messi_minidist_range_raw(isax_index *index,
                                 float *paa_or_fft,
                                 sax_type *sax_min,
                                 sax_type *sax_max,
                                 sax_type *sax_cardinalities,
                                 float bsf) {
    return messi_minidist_range_raw_partitioned(index, paa_or_fft, sax_min, sax_max,
                                                 sax_cardinalities, bsf,
                                                 index->settings->n_segments, NULL);
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

/* The query transform is fixed while scanning a leaf.  With full 8-bit
 * symbols, cache each dimension/symbol contribution once rather than
 * repeatedly reconstructing its bin interval for every record. */
int messi_build_record_lb_table(const isax_index *index,
                                const float *paa_or_fft,
                                int dimensions,
                                float table[MESSI_RECORD_LB_MAX_DIMENSIONS][256]) {
    if (index == NULL || index->settings == NULL || paa_or_fft == NULL || table == NULL ||
        dimensions < 16 || dimensions > MESSI_RECORD_LB_MAX_DIMENSIONS ||
        index->settings->sax_bit_cardinality != 8 ||
        index->settings->sax_alphabet_cardinality != 256) return 0;

    const isax_index_settings *settings = index->settings;
    const int function_type = settings->function_type;
    int sfa_start = 0;
    if (function_type == 4 && !settings->is_norm &&
        (settings->n_coefficients == 0 ||
         (index->coefficients != NULL && index->coefficients[0] == 0))) {
        sfa_start = settings->n_coefficients == 0 ? 2 : 1;
    }

    const ts_type *sax_bins = NULL;
    float sax_factor = 1.0f;
    if (function_type != 4 && function_type != 5 && function_type != 6 && function_type != 7) {
        const int offset = ((settings->sax_alphabet_cardinality - 1) *
                            (settings->sax_alphabet_cardinality - 2)) / 2;
        sax_bins = sax_breakpoints + offset;
        sax_factor = settings->mindist_sqrt;
    }

    for (int dimension = 0; dimension < dimensions; ++dimension) {
        const int ignored = function_type == 4 &&
                            dimension > 0 && dimension < sfa_start;
        const ts_type *bins = sax_bins;
        float factor = sax_factor;
        if (function_type == 4 || function_type == 5 || function_type == 6 || function_type == 7) {
            if (index->bins == NULL || index->bins[dimension] == NULL) return 0;
            bins = index->bins[dimension];
            factor = function_type == 4
                         ? (dimension == 0 && sfa_start > 0 ? 1.0f : 2.0f)
                         : 1.0f;
        }
        for (int symbol = 0; symbol < 256; ++symbol)
            table[dimension][symbol] = ignored ? 0.0f :
                get_lb_distance(bins, paa_or_fft[dimension], (sax_type) symbol, 8, 8,
                                256, factor, 1);
    }
    return 1;
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

    if (index->settings->function_type == 4 && !index->settings->is_norm &&
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
    const float factor = index->settings->function_type == 4 ? 2.0f : 1.0f;
        distance += get_lb_distance(index->bins[i], fft[i], v, c_c, c_m, max_cardinality, factor, 0);
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
