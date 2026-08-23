#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "config.h"
#include "globals.h"
#include "ads/calc_utils.h"
#include "ads/lower_bound_simd.h"
#include "ads/spartan/pca.h"
#include "ads/sax/sax.h"
#include "ads/sax/ts.h"
#include "ads/spartan/spartan.h"
#include "ads/sfa/sfa.h"
typedef struct spartan_bins_data {
    isax_index *index;
    ts_type **coeff_mem_array;
    long int start_number;
    long int stop_number;
} spartan_bins_data;

static uint64_t spartan_rng_next(uint64_t *state) {
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static long spartan_random_at_most(uint64_t *state, long max) {
    if (max <= 0) return 0;
    const uint64_t range = (uint64_t) max + 1;
    const uint64_t limit = UINT64_MAX - (UINT64_MAX % range);
    uint64_t value;
    do {
        value = spartan_rng_next(state);
    } while (value >= limit);
    return (long) (value % range);
}

static int spartan_compare_positions(const void *left, const void *right) {
    const long int a = *(const long int *) left;
    const long int b = *(const long int *) right;
    return (a > b) - (a < b);
}

enum response spartan_bins_init(isax_index *index) {
    int num_symbols = index->settings->sax_alphabet_cardinality;
    int n_segments = index->settings->n_segments;

    index->bins = NULL;
    index->bins = (ts_type **) calloc(n_segments, sizeof(ts_type *));
    index->binsv = (ts_type *) calloc(n_segments * (num_symbols - 1), sizeof(ts_type));

    // allocate num_symbols-1 memory slots for each word
    for (int i = 0; i < n_segments; ++i) {
        index->bins[i] = calloc(num_symbols - 1, sizeof(ts_type));
        for (int j = 0; j < num_symbols - 1; ++j) {
            index->bins[i][j] = FLT_MAX;
        }
    }
    for (int j = 0; j < n_segments * (num_symbols - 1); ++j) {
        index->binsv[j] = FLT_MAX;
    }
    fprintf(stderr, ">>> SPARTAN: Initialized bins[%d][%d] \n", n_segments, num_symbols - 1);

    return SUCCESS;
}

void spartan_free_bins(isax_index *index) {
    if (index == NULL || index->bins == NULL) {
        pca_free(index);
        return;
    }
    for (int i = 0; i < index->settings->n_segments; ++i) {
        free(index->bins[i]);
    }
    free(index->bins);
    pca_free(index);
}

static enum response spartan_collect_samples(isax_index *index, const char *ifilename,
                                             long int ts_num, int filetype_int,
                                             int apply_znorm, ts_type *samples,
                                             unsigned int sample_size) {
    if (sample_size == 0) {
        return FAILURE;
    }
    FILE *ifile = fopen(ifilename, "rb");
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        return FAILURE;
    }

    unsigned long ts_length = index->settings->timeseries_size;
    file_type *ts_orig1 = NULL;
    if (filetype_int) {
        ts_orig1 = (file_type *) calloc(ts_length, sizeof(file_type));
    }

    unsigned int records = sample_size;
    if ((long int) records > ts_num) {
        records = (unsigned int) ts_num;
    }

    long int *positions = NULL;
    if (index->settings->sample_type == 3) {
        uint64_t rng_state = (uint64_t) index->settings->sampling_seed;
        positions = malloc(sizeof(long int) * records);
        for (unsigned int i = 0; i < records; ++i) {
            positions[i] = (long int) i;
        }
        for (long int i = (long int) records; i < ts_num; ++i) {
            long int j = spartan_random_at_most(&rng_state, i);
            if (j < (long int) records) {
                positions[j] = i;
            }
        }
        /* Sample order does not affect PCA or binning.  Sorting lets the
         * collection pass use sequential I/O instead of one seek per record. */
        qsort(positions, records, sizeof(*positions), spartan_compare_positions);
    }

    if (index->settings->sample_type == 2 || index->settings->sample_type == 3) {
        /* A seek after every selected record turns a one-million-record sample
         * into one million tiny I/O requests.  Read contiguous blocks instead
         * and copy the uniform or reservoir-selected records within each one. */
        const unsigned long block_records = 8192;
        const size_t block_values = (size_t) block_records * ts_length;
        ts_type *float_block = filetype_int ? NULL : malloc(block_values * sizeof(*float_block));
        file_type *int_block = filetype_int ? malloc(block_values * sizeof(*int_block)) : NULL;
        if ((!filetype_int && float_block == NULL) || (filetype_int && int_block == NULL)) {
            free(float_block); free(int_block); free(ts_orig1); free(positions); fclose(ifile);
            return FAILURE;
        }

        unsigned int sample = 0;
        for (unsigned long first = 0; first < (unsigned long) ts_num; first += block_records) {
            const unsigned long count = ((unsigned long) ts_num - first) < block_records
                                            ? ((unsigned long) ts_num - first) : block_records;
            const size_t values = (size_t) count * ts_length;
            const size_t got = filetype_int
                                   ? fread(int_block, sizeof(*int_block), values, ifile)
                                   : fread(float_block, sizeof(*float_block), values, ifile);
            if (got != values) {
                free(float_block); free(int_block); free(ts_orig1); free(positions); fclose(ifile);
                return FAILURE;
            }
            while (sample < records) {
                const unsigned long position = index->settings->sample_type == 2
                    ? (unsigned long) (((unsigned long long) sample *
                                        (unsigned long long) ts_num) / records)
                    : (unsigned long) positions[sample];
                if (position < first) { ++sample; continue; }
                if (position >= first + count) break;
                ts_type *destination = samples + (size_t) sample * ts_length;
                const size_t offset = (size_t) (position - first) * ts_length;
                if (filetype_int) {
                    for (unsigned long j = 0; j < ts_length; ++j)
                        destination[j] = (ts_type) int_block[offset + j];
                } else {
                    memcpy(destination, float_block + offset, ts_length * sizeof(*destination));
                }
                ++sample;
            }
        }
        free(float_block);
        free(int_block);
        if (sample != records) {
            free(ts_orig1); free(positions); fclose(ifile);
            return FAILURE;
        }
    } else {
        for (unsigned int i = 0; i < records; ++i) {
            ts_type *ts = samples + (i * ts_length);
            if (filetype_int) {
                if (fread(ts_orig1, sizeof(file_type), ts_length, ifile) != ts_length) {
                    free(ts_orig1); free(positions); fclose(ifile); return FAILURE;
                }
                for (unsigned long j = 0; j < ts_length; ++j) ts[j] = (ts_type) ts_orig1[j];
            } else if (fread(ts, sizeof(ts_type), ts_length, ifile) != ts_length) {
                free(ts_orig1); free(positions); fclose(ifile); return FAILURE;
            }
        }
    }

    if (apply_znorm) {
        #pragma omp parallel for schedule(static)
        for (unsigned int i = 0; i < records; ++i) {
            znorm(samples + (i * ts_length), ts_length);
        }
    }

    fclose(ifile);
    free(ts_orig1);
    free(positions);

    return SUCCESS;
}

static void *spartan_order_divide_worker(void *transferdata) {
    spartan_bins_data *bins_data = (spartan_bins_data *) transferdata;
    ts_type **coeff_mem_array = bins_data->coeff_mem_array;

    isax_index *index = bins_data->index;
    long int start_number = bins_data->start_number;
    long int stop_number = bins_data->stop_number;

    unsigned int sample_size = index->settings->sample_size;
    int n_segments = index->settings->n_segments;
    ts_type *cur_coeff_line;

    for (int j = start_number; j < stop_number; ++j) {
        cur_coeff_line = (ts_type *) coeff_mem_array[j];
        qsort(cur_coeff_line, sample_size, sizeof(ts_type), &compare_ts_type);
    }

    // equi-depth splitting
    if (index->settings->histogram_type == 1) {
        int num_symbols = index->settings->sax_alphabet_cardinality;
        ts_type depth = (ts_type) sample_size / num_symbols;

        for (int i = start_number; i < stop_number; ++i) {
            float bin_index = 0.0;
            cur_coeff_line = coeff_mem_array[i];
            for (int j = 0; j < num_symbols - 1; ++j) {
                bin_index += depth;
                index->bins[i][j] = cur_coeff_line[(int) bin_index];
            }
        }
    }
    // equi-width splitting
    else if (index->settings->histogram_type == 2) {
        int num_symbols = index->settings->sax_alphabet_cardinality;

        for (int i = start_number; i < stop_number; ++i) {
            cur_coeff_line = coeff_mem_array[i];
            ts_type first = cur_coeff_line[0];
            ts_type last = cur_coeff_line[sample_size - 1];
            ts_type interval_width = (last - first) / (ts_type) num_symbols;
            for (int j = 0; j < num_symbols - 1; ++j) {
                index->bins[i][j] = interval_width * (j + 1) + first;
            }
        }
    }

    if (n_segments == 0) {
        fprintf(stderr, "warning: SPARTAN has zero segments.\n");
    }
    return NULL;
}

enum response spartan_set_bins(isax_index *index, const char *ifilename, long int ts_num,
                               int maxquerythread, int filetype_int, int apply_znorm) {
    if (index == NULL || index->settings == NULL) {
        return FAILURE;
    }
    int dim = index->settings->n_segments;
    int ts_length = index->settings->timeseries_size;
    unsigned int sample_size = index->settings->sample_size;
    if (dim <= 0 || sample_size == 0 || ts_num <= 0 || sample_size > (unsigned long) ts_num) {
        fprintf(stderr, "error: invalid SPARTAN sampling settings.\n");
        return FAILURE;
    }

    int worker_threads = maxquerythread;
    if (worker_threads < 1) {
        worker_threads = 1;
    }
    if (worker_threads > dim) {
        worker_threads = dim;
    }

    fprintf(stderr, ">>> SPARTAN binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START
    double binning_start = messi_monotonic_seconds();

    ts_type *samples = calloc((size_t) sample_size * ts_length, sizeof(ts_type));
    if (samples == NULL) {
        fprintf(stderr, "error: failed to allocate SPARTAN sample buffer.\n");
        return FAILURE;
    }

    if (spartan_collect_samples(index, ifilename, ts_num, filetype_int, apply_znorm, samples, sample_size) != SUCCESS) {
        free(samples);
        return FAILURE;
    }
    double sample_end = messi_monotonic_seconds();

    pca_free(index);
    if (pca_fit(index, samples, sample_size, ts_length) != SUCCESS) {
        free(samples);
        return FAILURE;
    }
    double pca_end = messi_monotonic_seconds();

    ts_type **coeff_mem_array = (ts_type **) calloc(dim, sizeof(ts_type *));
    if (coeff_mem_array == NULL) {
        free(samples);
        pca_free(index);
        return FAILURE;
    }
    for (int k = 0; k < dim; ++k) {
        coeff_mem_array[k] = (ts_type *) calloc(sample_size, sizeof(ts_type));
        if (coeff_mem_array[k] == NULL) {
            for (int j = 0; j < k; ++j) {
                free(coeff_mem_array[j]);
            }
            free(coeff_mem_array);
            free(samples);
            pca_free(index);
            return FAILURE;
        }
    }

    ts_type *projection = calloc(dim, sizeof(ts_type));
    if (projection == NULL) {
        for (int k = 0; k < dim; ++k) {
            free(coeff_mem_array[k]);
        }
        free(coeff_mem_array);
        free(samples);
        pca_free(index);
        return FAILURE;
    }
    for (unsigned int i = 0; i < sample_size; ++i) {
        const ts_type *row = samples + (i * ts_length);
        if (pca_from_ts(index, row, projection) != SUCCESS) {
            free(projection);
            for (int k = 0; k < dim; ++k) {
                free(coeff_mem_array[k]);
            }
            free(coeff_mem_array);
            free(samples);
            pca_free(index);
            return FAILURE;
        }
        for (int k = 0; k < dim; ++k) {
            coeff_mem_array[k][i] = projection[k];
        }
    }
    free(projection);
    free(samples);
    double projection_end = messi_monotonic_seconds();

    if (configure_dynamic_bit_allocation(index, index->pca_explained_variance,
                                           index->settings->n_segments,
                                           index->settings->index_type == MESSI_INDEX_TRIE &&
                                                   index->settings->trie_dynamic_alphabet
                                               ? index->settings->trie_alphabet_budget_bits * index->settings->n_segments
                                               : (index->settings->n_segments < (int) (sizeof(root_mask_type) * 8)
                                                      ? index->settings->n_segments
                                                      : (int) (sizeof(root_mask_type) * 8)),
                                           index->settings->index_type == MESSI_INDEX_TRIE &&
                                                   index->settings->trie_dynamic_alphabet
                                               ? index->settings->trie_min_bits : 0,
                                           index->settings->index_type == MESSI_INDEX_TRIE &&
                                                   index->settings->trie_dynamic_alphabet
                                               ? index->settings->trie_max_bits
                                               : index->settings->sax_bit_cardinality) != SUCCESS) {
        for (int k = 0; k < dim; ++k) {
            free(coeff_mem_array[k]);
        }
        free(coeff_mem_array);
        pca_free(index);
        return FAILURE;
    }

    pthread_t threadid[worker_threads];
    spartan_bins_data *input_data = malloc(sizeof(spartan_bins_data) * (size_t) worker_threads);
    if (input_data == NULL) {
        for (int k = 0; k < dim; ++k) {
            free(coeff_mem_array[k]);
        }
        free(coeff_mem_array);
        pca_free(index);
        return FAILURE;
    }

    for (int i = 0; i < worker_threads; i++) {
        input_data[i].index = index;
        input_data[i].coeff_mem_array = coeff_mem_array;
        input_data[i].start_number = i * (dim / worker_threads);
        input_data[i].stop_number = (i + 1) * (dim / worker_threads);
    }

    input_data[worker_threads - 1].start_number = (worker_threads - 1) * (dim / worker_threads);
    input_data[worker_threads - 1].stop_number = dim;

    for (int i = 0; i < worker_threads; i++) {
        pthread_create(&(threadid[i]), NULL, spartan_order_divide_worker, (void *) &(input_data[i]));
    }

    for (int i = 0; i < worker_threads; i++) {
        pthread_join(threadid[i], NULL);
    }
    double bins_end = messi_monotonic_seconds();

    for (int k = 0; k < dim; ++k) {
        free(coeff_mem_array[k]);
    }
    free(coeff_mem_array);
    free(input_data);

    COUNT_BINNING_TIME_END

    spartan_print_bins(index);
    fprintf(stderr,
            ">>> SPARTAN binning timing: sample=%.3fs PCA=%.3fs projection=%.3fs "
            "bin-sort=%.3fs total=%.3fs\n",
            sample_end - binning_start, pca_end - sample_end,
            projection_end - pca_end, bins_end - projection_end,
            bins_end - binning_start);
    fprintf(stderr, ">>> Finished SPARTAN binning\n");
    return SUCCESS;
}

/* Return the first boundary strictly above value.  This is equivalent to the
 * former linear scan, including its treatment of values equal to a boundary,
 * but reduces 256-symbol quantization from up to 255 comparisons to eight. */
static unsigned int spartan_symbol_from_value(const ts_type *boundaries, int count,
                                              ts_type value) {
    int low = 0, high = count;
    while (low < high) {
        const int middle = low + (high - low) / 2;
        if (value < boundaries[middle]) high = middle;
        else low = middle + 1;
    }
    return (unsigned int) low;
}

void spartan_from_pca(isax_index *index, const ts_type *coeffs, sax_type *sax_out) {
    const int n_segments = index->settings->n_segments;
    const int boundary_count = index->settings->sax_alphabet_cardinality - 1;
    for (int k = 0; k < n_segments; ++k)
        sax_out[k] = (sax_type) spartan_symbol_from_value(index->bins[k], boundary_count, coeffs[k]);
}

enum response spartan_from_ts(isax_index *index, const ts_type *ts, sax_type *sax_out,
                              ts_type *coeff_scratch) {
    if (index == NULL || ts == NULL || sax_out == NULL || coeff_scratch == NULL) {
        return FAILURE;
    }

    if (pca_from_ts(index, ts, coeff_scratch) != SUCCESS) {
        return FAILURE;
    }
    spartan_from_pca(index, coeff_scratch, sax_out);
    return SUCCESS;
}

enum response pca_from_ts(const isax_index *index, const ts_type *ts, ts_type *out) {
    if (index == NULL || ts == NULL || out == NULL) {
        return FAILURE;
    }
    int input_dim = index->pca_dim;
    int output_dim = index->settings->n_segments;
    int components = index->pca_components_count;
    if (components <= 0 || index->pca_components == NULL || index->pca_mean == NULL || input_dim <= 0) {
        for (int i = 0; i < output_dim; ++i) {
            out[i] = 0.0f;
        }
        return FAILURE;
    }

    for (int k = 0; k < components; ++k) {
        double acc = 0.0;
        const ts_type *component = index->pca_components + (k * input_dim);
        for (int i = 0; i < input_dim; ++i) {
            // Temporary: disable PCA mean-centering.
            double centered = (double) ts[i] - index->pca_mean[i];
            acc += centered * component[i];

            // acc += (double) ts[i] * component[i];
        }
        out[k] = (ts_type) acc;
    }
    for (int k = components; k < output_dim; ++k) {
        out[k] = 0.0f;
    }
    return SUCCESS;
}



void spartan_print_bins(isax_index *index) {
    fprintf(stderr, ">>> SPARTAN: Sample size %u\n", index->settings->sample_size);
    if (index->settings->histogram_type == 1) {
        fprintf(stderr, ">>> SPARTAN: Using Equi-depth histograms\n");
    } else if (index->settings->histogram_type == 2) {
        fprintf(stderr, ">>> SPARTAN: Using Equi-width histograms\n");
    }
}
