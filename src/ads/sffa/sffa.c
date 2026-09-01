#include "config.h"
#include "ads/sffa/sffa.h"

#include "ads/calc_utils.h"
#include "ads/sax/ts.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double variance;
    int coordinate;
} sffa_rank;

static int compare_float(const void *a, const void *b) {
    const ts_type x = *(const ts_type *) a, y = *(const ts_type *) b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static int compare_rank(const void *a, const void *b) {
    const sffa_rank *x = a, *y = b;
    if (x->variance > y->variance) return -1;
    if (x->variance < y->variance) return 1;
    return x->coordinate - y->coordinate;
}

static uint64_t next_random(uint64_t *state) {
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static unsigned int symbol_from_value(const ts_type *boundaries, int count, ts_type value) {
    int low = 0, high = count;
    while (low < high) {
        const int middle = low + (high - low) / 2;
        if (value < boundaries[middle]) high = middle;
        else low = middle + 1;
    }
    return (unsigned int) low;
}

enum response sffa_bins_init(isax_index *index) {
    if (index == NULL || index->settings == NULL) return FAILURE;
    const int dimensions = index->settings->n_segments;
    const int boundaries = index->settings->sax_alphabet_cardinality - 1;
    if (dimensions <= 0 || boundaries <= 0) return FAILURE;
    index->bins = calloc((size_t) dimensions, sizeof(*index->bins));
    index->binsv = calloc((size_t) dimensions * boundaries, sizeof(*index->binsv));
    index->sffa_coordinates = calloc((size_t) dimensions, sizeof(*index->sffa_coordinates));
    if (index->bins == NULL || index->binsv == NULL || index->sffa_coordinates == NULL) {
        sffa_free_bins(index);
        return FAILURE;
    }
    for (int i = 0; i < dimensions; ++i) {
        index->bins[i] = index->binsv + (size_t) i * boundaries;
    }
    index->sffa_scale = 1.0f;
    return SUCCESS;
}

void sffa_free_bins(isax_index *index) {
    if (index == NULL) return;
    free(index->bins);
    free(index->binsv);
    free(index->sffa_coordinates);
    index->bins = NULL;
    index->binsv = NULL;
    index->sffa_coordinates = NULL;
    index->sffa_scale = 0.0f;
}

enum response sffa_project(const isax_index *index, const ts_type *ts,
                           ts_type *out, sffa_workspace *ws) {
    if (index == NULL || index->settings == NULL || out == NULL ||
        index->sffa_coordinates == NULL || sffa_fractional_transform(index, ts, ws) != SUCCESS)
        return FAILURE;
    const int dimensions = index->settings->n_segments;
    const int n = index->settings->timeseries_size;
    for (int i = 0; i < dimensions; ++i) {
        const sffa_coordinate coordinate = index->sffa_coordinates[i];
        if (coordinate.coefficient < 0 || coordinate.coefficient >= n) return FAILURE;
        out[i] = index->sffa_scale * ws->output[coordinate.coefficient][coordinate.imaginary ? 1 : 0];
    }
    return SUCCESS;
}

void sffa_from_values(const isax_index *index, const ts_type *values, sax_type *word) {
    const int count = index->settings->sax_alphabet_cardinality - 1;
    for (int i = 0; i < index->settings->n_segments; ++i)
        word[i] = (sax_type) symbol_from_value(index->bins[i], count, values[i]);
}

enum response sffa_from_ts(const isax_index *index, const ts_type *ts,
                           sax_type *word, sffa_workspace *ws, ts_type *values) {
    if (word == NULL || sffa_project(index, ts, values, ws) != SUCCESS) return FAILURE;
    sffa_from_values(index, values, word);
    return SUCCESS;
}

/* A conservative real-operator norm for selected Re/Im rows.  We calculate
 * G=B B^T exactly enough for float input and use max row sum of |G| as a
 * certified upper bound on ||B||_2^2.  The tiny margin absorbs roundoff. */
static enum response sffa_choose_scale(isax_index *index, sffa_workspace *ws) {
    const int d = index->settings->n_segments, n = index->settings->timeseries_size;
    ts_type *rows = calloc((size_t) d * n, sizeof(*rows));
    ts_type *basis = calloc((size_t) n, sizeof(*basis));
    if (rows == NULL || basis == NULL) { free(rows); free(basis); return FAILURE; }
    for (int col = 0; col < n; ++col) {
        basis[col] = 1.0f;
        if (sffa_fractional_transform(index, basis, ws) != SUCCESS) {
            free(rows); free(basis); return FAILURE;
        }
        for (int row = 0; row < d; ++row) {
            const sffa_coordinate c = index->sffa_coordinates[row];
            rows[(size_t) row * n + col] = ws->output[c.coefficient][c.imaginary ? 1 : 0];
        }
        basis[col] = 0.0f;
    }
    double max_row_sum = 0.0;
    for (int i = 0; i < d; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < d; ++j) {
            double dot = 0.0;
            for (int k = 0; k < n; ++k)
                dot += (double) rows[(size_t) i * n + k] * rows[(size_t) j * n + k];
            row_sum += fabs(dot);
        }
        if (row_sum > max_row_sum) max_row_sum = row_sum;
    }
    free(rows); free(basis);
    if (max_row_sum <= 0.0) return FAILURE;
    index->sffa_scale = (ts_type) (0.9999 / sqrt(max_row_sum));
    return SUCCESS;
}

static enum response sffa_fit_from_sample(isax_index *index, const ts_type *sample,
                                          unsigned int sample_size, sffa_workspace *ws) {
    const int n = index->settings->timeseries_size;
    const int d = index->settings->n_segments;
    const int pool_size = 2 * n;
    if (sample == NULL || sample_size < 2 || (size_t) pool_size > SIZE_MAX / sample_size)
        return FAILURE;
    ts_type *sampled = calloc((size_t) pool_size * sample_size, sizeof(*sampled));
    sffa_rank *rank = calloc((size_t) pool_size, sizeof(*rank));
    if (sampled == NULL || rank == NULL) goto failed;
    for (unsigned int row = 0; row < sample_size; ++row) {
        if (sffa_fractional_transform(index, sample + (size_t) row * n, ws) != SUCCESS) goto failed;
        for (int coefficient = 0; coefficient < n; ++coefficient) {
            sampled[(size_t) coefficient * sample_size + row] = ws->output[coefficient][0];
            sampled[(size_t) (n + coefficient) * sample_size + row] = ws->output[coefficient][1];
        }
    }
    for (int coordinate = 0; coordinate < pool_size; ++coordinate) {
        double mean = 0.0, variance = 0.0;
        const ts_type *values = sampled + (size_t) coordinate * sample_size;
        for (unsigned int row = 0; row < sample_size; ++row) mean += values[row];
        mean /= sample_size;
        for (unsigned int row = 0; row < sample_size; ++row) {
            const double delta = values[row] - mean;
            variance += delta * delta;
        }
        rank[coordinate].variance = variance / sample_size;
        rank[coordinate].coordinate = coordinate;
    }
    qsort(rank, (size_t) pool_size, sizeof(*rank), compare_rank);
    for (int i = 0; i < d; ++i) {
        const int coordinate = rank[i].coordinate;
        index->sffa_coordinates[i].coefficient = coordinate % n;
        index->sffa_coordinates[i].imaginary = coordinate >= n;
    }
    if (sffa_choose_scale(index, ws) != SUCCESS) goto failed;
    free(index->settings->symbolic_variances);
    index->settings->symbolic_variances = calloc((size_t) d, sizeof(*index->settings->symbolic_variances));
    if (index->settings->symbolic_variances == NULL) goto failed;
    for (int i = 0; i < d; ++i) {
        const int coordinate = rank[i].coordinate;
        ts_type *values = sampled + (size_t) coordinate * sample_size;
        for (unsigned int row = 0; row < sample_size; ++row) values[row] *= index->sffa_scale;
        index->settings->symbolic_variances[i] = rank[i].variance * index->sffa_scale * index->sffa_scale;
        qsort(values, sample_size, sizeof(*values), compare_float);
        const int alphabet = index->settings->sax_alphabet_cardinality;
        for (int bin = 0; bin < alphabet - 1; ++bin) {
            if (index->settings->histogram_type == 1) {
                unsigned int at = (unsigned int) (((unsigned long long) (bin + 1) * sample_size) / alphabet);
                if (at >= sample_size) at = sample_size - 1;
                index->bins[i][bin] = values[at];
            } else {
                index->bins[i][bin] = values[0] +
                    (values[sample_size - 1] - values[0]) * (ts_type) (bin + 1) / alphabet;
            }
        }
    }
    if (configure_dynamic_bit_allocation(index, index->settings->symbolic_variances, d,
                                         d, 0, index->settings->sax_bit_cardinality) != SUCCESS)
        goto failed;
    free(rank); free(sampled);
    return SUCCESS;
failed:
    free(rank); free(sampled);
    return FAILURE;
}

static double sffa_score_holdout(isax_index *index, const ts_type *validation,
                                 unsigned int validation_size, sffa_workspace *ws) {
    const int n = index->settings->timeseries_size;
    const int d = index->settings->n_segments;
    const int alphabet = index->settings->sax_alphabet_cardinality;
    ts_type *values = calloc((size_t) validation_size * d, sizeof(*values));
    sax_type *words = calloc((size_t) validation_size * d, sizeof(*words));
    if (values == NULL || words == NULL) { free(values); free(words); return -1.0; }
    for (unsigned int row = 0; row < validation_size; ++row) {
        ts_type *row_values = values + (size_t) row * d;
        if (sffa_project(index, validation + (size_t) row * n, row_values, ws) != SUCCESS) {
            free(values); free(words); return -1.0;
        }
        for (int dimension = 0; dimension < d; ++dimension)
            words[(size_t) row * d + dimension] =
                (sax_type) symbol_from_value(index->bins[dimension], alphabet - 1, row_values[dimension]);
    }

    double bound_total = 0.0, exact_total = 0.0;
    unsigned int pairs = 0;
    for (unsigned int left = 0; left < validation_size && pairs < 4096U; ++left) {
        for (unsigned int right = left + 1; right < validation_size && pairs < 4096U; ++right) {
            double exact = 0.0, bound = 0.0;
            const ts_type *left_ts = validation + (size_t) left * n;
            const ts_type *right_ts = validation + (size_t) right * n;
            const ts_type *left_values = values + (size_t) left * d;
            const sax_type *right_word = words + (size_t) right * d;
            for (int i = 0; i < n; ++i) {
                const double delta = left_ts[i] - right_ts[i];
                exact += delta * delta;
            }
            if (exact <= 1e-12) continue;
            for (int dimension = 0; dimension < d; ++dimension) {
                const int symbol = right_word[dimension];
                const ts_type lower = symbol == 0 ? -FLT_MAX : index->bins[dimension][symbol - 1];
                const ts_type upper = symbol == alphabet - 1 ? FLT_MAX : index->bins[dimension][symbol];
                double delta = 0.0;
                if (left_values[dimension] < lower) delta = lower - left_values[dimension];
                else if (left_values[dimension] > upper) delta = left_values[dimension] - upper;
                bound += delta * delta;
            }
            /* The order-selection objective must itself remain a valid
             * lower bound on the raw squared distance. */
            if (bound > exact + 1e-5 * fmax(1.0, exact)) {
                free(values); free(words); return -1.0;
            }
            bound_total += bound;
            exact_total += exact;
            ++pairs;
        }
    }
    free(values); free(words);
    return pairs == 0 || exact_total <= 0.0 ? -1.0 : bound_total / exact_total;
}

static int sffa_score_is_better(double score, double order,
                                double best_score, double best_order) {
    const double tolerance = 1e-12;
    if (score > best_score + tolerance) return 1;
    if (fabs(score - best_score) <= tolerance && fabs(order - 1.0) < fabs(best_order - 1.0)) return 1;
    return 0;
}

static enum response sffa_select_order(isax_index *index, const ts_type *sample,
                                       unsigned int sample_size, sffa_workspace *ws) {
    const unsigned int fit_size = (sample_size * 3U) / 4U;
    const unsigned int validation_size = sample_size - fit_size;
    double best_order = 1.0, best_score = -1.0;
    if (fit_size < 2 || validation_size < 2) return FAILURE;
    for (int tenths = 2; tenths <= 10; ++tenths) {
        const double order = tenths / 10.0;
        index->settings->sffa_order = order;
        if (sffa_fit_from_sample(index, sample, fit_size, ws) != SUCCESS) return FAILURE;
        const double score = sffa_score_holdout(index, sample + (size_t) fit_size * index->settings->timeseries_size,
                                                validation_size, ws);
        if (score < 0.0) return FAILURE;
        if (sffa_score_is_better(score, order, best_score, best_order)) {
            best_order = order;
            best_score = score;
        }
    }
    for (int direction = -1; direction <= 1; direction += 2) {
        const double order = best_order + direction * 0.05;
        if (order < 0.2 - 1e-12 || order > 1.0 + 1e-12) continue;
        index->settings->sffa_order = order;
        if (sffa_fit_from_sample(index, sample, fit_size, ws) != SUCCESS) return FAILURE;
        const double score = sffa_score_holdout(index, sample + (size_t) fit_size * index->settings->timeseries_size,
                                                validation_size, ws);
        if (score < 0.0) return FAILURE;
        if (sffa_score_is_better(score, order, best_score, best_order)) {
            best_order = order;
            best_score = score;
        }
    }
    index->settings->sffa_order = best_order;
    if (sffa_fit_from_sample(index, sample, sample_size, ws) != SUCCESS) return FAILURE;
    fprintf(stderr, ">>> SFFA: auto selected order %.6g (held-out symbolic tightness %.6g)\n",
            best_order, best_score);
    return SUCCESS;
}

enum response sffa_set_bins(isax_index *index, const char *filename,
                            long ts_num, int maxquerythread,
                            int filetype_int, int apply_znorm) {
    (void) maxquerythread;
    if (index == NULL || index->settings == NULL || filename == NULL ||
        index->bins == NULL || index->sffa_coordinates == NULL || ts_num <= 0)
        return FAILURE;
    const int n = index->settings->timeseries_size;
    const int d = index->settings->n_segments;
    const unsigned int sample_size = index->settings->sample_size;
    const int pool_size = 2 * n;
    if (n <= 0 || d <= 0 || d > pool_size || sample_size == 0 || sample_size > (unsigned long) ts_num)
        return FAILURE;
    if (index->settings->sffa_auto_order && sample_size < 8) {
        fprintf(stderr, "error: SFFA automatic order selection needs at least 8 training samples.\n");
        return FAILURE;
    }
    if (sample_size < (unsigned int) index->settings->sax_alphabet_cardinality * 8U)
        fprintf(stderr, "warning: SFFA sample size %u is small for alphabet %d; use at least 8 samples/bin when possible.\n",
                sample_size, index->settings->sax_alphabet_cardinality);
    if ((size_t) n > SIZE_MAX / sample_size) return FAILURE;

    ts_type *sample = calloc((size_t) n * sample_size, sizeof(*sample));
    ts_type *swap_row = calloc((size_t) n, sizeof(*swap_row));
    file_type *raw = filetype_int ? calloc((size_t) n, sizeof(*raw)) : NULL;
    sffa_workspace ws = {0};
    FILE *input = NULL;
    if (sample == NULL || swap_row == NULL || (filetype_int && raw == NULL)) goto failed;
    input = fopen(filename, "rb");
    if (input == NULL) goto failed;
    sffa_workspace_init(&ws, (unsigned long) n);
    if (ws.input == NULL || ws.forward_plan == NULL || ws.inverse_plan == NULL) goto failed;
    uint64_t rng = ((uint64_t) index->settings->sampling_seed << 32) | 1U;
    for (unsigned int row = 0; row < sample_size; ++row) {
        unsigned long record = row;
        ts_type *destination = sample + (size_t) row * n;
        if (index->settings->sample_type == 2)
            record = ((unsigned long long) row * (unsigned long) ts_num) / sample_size;
        else if (index->settings->sample_type == 3)
            record = (unsigned long) (next_random(&rng) % (uint64_t) ts_num);
        const size_t item_size = filetype_int ? sizeof(*raw) : sizeof(*destination);
        if (fseek(input, (long) (record * (unsigned long) n * item_size), SEEK_SET) != 0) goto failed;
        if (filetype_int) {
            if (fread(raw, sizeof(*raw), (size_t) n, input) != (size_t) n) goto failed;
            for (int i = 0; i < n; ++i) destination[i] = (ts_type) raw[i];
        } else if (fread(destination, sizeof(*destination), (size_t) n, input) != (size_t) n) goto failed;
        if (apply_znorm) znorm(destination, n);
    }
    /* Keep fit and validation independent even for sequential/stratified sampling. */
    rng ^= UINT64_C(0xd1b54a32d192ed03);
    for (unsigned int row = sample_size - 1; row > 0; --row) {
        const unsigned int other = (unsigned int) (next_random(&rng) % (uint64_t) (row + 1));
        if (other == row) continue;
        memcpy(swap_row, sample + (size_t) row * n, (size_t) n * sizeof(*sample));
        memcpy(sample + (size_t) row * n, sample + (size_t) other * n, (size_t) n * sizeof(*sample));
        memcpy(sample + (size_t) other * n, swap_row, (size_t) n * sizeof(*sample));
    }
    if (index->settings->sffa_auto_order) {
        if (sffa_select_order(index, sample, sample_size, &ws) != SUCCESS) goto failed;
    } else if (sffa_fit_from_sample(index, sample, sample_size, &ws) != SUCCESS) {
        goto failed;
    }
    fprintf(stderr, ">>> SFFA: order %.6g, selected %d scalar Re/Im coordinates, scale %.6g\n",
            index->settings->sffa_order, d, index->sffa_scale);
    fclose(input); sffa_workspace_destroy(&ws); free(raw); free(swap_row); free(sample);
    return SUCCESS;
failed:
    if (input) fclose(input);
    sffa_workspace_destroy(&ws); free(raw); free(swap_row); free(sample);
    return FAILURE;
}
