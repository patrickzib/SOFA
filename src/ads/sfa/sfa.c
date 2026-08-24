//  
//  sfa.c
//  sfa  C version for MESSI
//
//  Based on sfa_trie code by Karima Echihabi on 18/11/2017
//  Copyright 2017 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#include "globals.h"
#include <stdio.h>
#include <pthread.h>
#ifdef VALUES

#include <values.h>

#endif

#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>

#include "ads/calc_utils.h"
#include "ads/lower_bound_simd.h"
#include "ads/isax_index.h"
#include "ads/sfa/dft.h"

#include "ads/sfa/sfa.h"

/*
  This functions allocates a two dimensional array
  of num_words rows and (num_symbols-1) columns
  The array will contain the discretization
  intervals
*/
enum response sfa_bins_init(isax_index *index) {
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
    fprintf(stderr, ">>> SFA: Initialized bins[%d][%d] \n", n_segments,
            num_symbols - 1);

    if (index->settings->n_coefficients != 0) {
        index->coefficients = calloc(n_segments / 2, sizeof(int));
    }

    return SUCCESS;
}

/*
  This functions frees the allocated bins-array
*/
void sfa_free_bins(isax_index *index) {
    for (int i = 0; i < index->settings->n_segments; ++i) {
        free(index->bins[i]);
    }
    free(index->bins);
}


/*
  In this function, the intervals are calculated (multiple coeff. binning).
  The coefficients with the highest variance or the first ones are chosen and
  these values are saved to bins
*/
enum response sfa_set_bins(
    isax_index *index, const char *ifilename,
    long int ts_num, int maxquerythread,
    int filetype_int, int apply_znorm) {
    const int n_segments = index->settings->n_segments;
    const int ts_length = index->settings->timeseries_size;
    const unsigned int sample_size = index->settings->sample_size;
    const int use_variance = index->settings->n_coefficients > 0;
    const int n_coefficients = use_variance
                                   ? index->settings->n_coefficients
                                   : n_segments;
    if (sample_size == 0 || ts_num <= 0 || sample_size > (unsigned long) ts_num) {
        fprintf(stderr, "error: SFA sample size must be between 1 and the dataset size.\n");
        return FAILURE;
    }

    int worker_threads = maxquerythread;
    if (worker_threads < 1) {
        worker_threads = 1;
    }
    if ((unsigned int) worker_threads > sample_size) {
        worker_threads = (int) sample_size;
    }

    fprintf(stderr, ">>> Binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START
    double binning_start = messi_monotonic_seconds();

    ts_type **dft_mem_array = calloc(n_coefficients, sizeof(*dft_mem_array));
    if (dft_mem_array == NULL) {
        return FAILURE;
    }
    for (int i = 0; i < n_coefficients; ++i) {
        dft_mem_array[i] = calloc(sample_size, sizeof(**dft_mem_array));
        if (dft_mem_array[i] == NULL) {
            free_dft_memory(index, i, dft_mem_array);
            return FAILURE;
        }
    }

    pthread_t threadid[worker_threads];
    bins_data_inmemory *input_data =
            malloc(worker_threads * sizeof(*input_data));
    if (input_data == NULL) {
        free_dft_memory(index, n_coefficients, dft_mem_array);
        return FAILURE;
    }

    /*
     * Phase 1:
     * sample time series and calculate DFT coefficients
     */
    const long sample_chunk = sample_size / worker_threads;
    const long ts_chunk = ts_num / worker_threads;

    for (int i = 0; i < worker_threads; ++i) {
        fftw_workspace fftw = {0};
        fftw_workspace_init(&fftw, ts_length);

        bins_data_inmemory *data = &input_data[i];

        data->index = index;
        data->dft_mem_array = dft_mem_array;
        data->filename = ifilename;
        data->workernumber = i;

        data->records = sample_chunk;
        data->records_offset = sample_chunk;

        switch (index->settings->sample_type) {
            case 1: /* first-n sampling */
                data->start_number = i * sample_chunk;
                data->stop_number = (i + 1) * sample_chunk;
                break;

            case 2: /* uniform sampling */
                data->start_number = i * ts_chunk;
                data->stop_number = (i + 1) * ts_chunk;
                break;

            case 3: /* random sampling */
                data->start_number = i * ts_chunk;
                data->stop_number = (i + 1) * ts_chunk;
                break;
        }

        data->filetype_int = filetype_int;
        data->apply_znorm = apply_znorm;
        data->status = SUCCESS;
        data->fftw = fftw;
    }

    /* Give any remainder to the last sampling worker. */
    input_data[worker_threads - 1].records =
            sample_size - (worker_threads - 1) * sample_chunk;

    if (index->settings->sample_type == 1) {
        input_data[worker_threads - 1].stop_number = sample_size;
    } else if (index->settings->sample_type == 2 || index->settings->sample_type == 3) {
        input_data[worker_threads - 1].stop_number = ts_num;
    }

    for (int i = 0; i < worker_threads; ++i) {
        pthread_create(
            &threadid[i],
            NULL,
            set_bins_worker_dft,
            &input_data[i]);
    }

    for (int i = 0; i < worker_threads; ++i) {
        pthread_join(threadid[i], NULL);
    }

    for (int i = 0; i < worker_threads; ++i) {
        if (input_data[i].status != SUCCESS) {
            for (int j = 0; j < worker_threads; ++j) {
                fftw_workspace_destroy(&input_data[j].fftw);
            }
            free(input_data);
            free_dft_memory(index, n_coefficients, dft_mem_array);
            return FAILURE;
        }
    }
    double sampling_end = messi_monotonic_seconds();

    /*
     * Optional variance-based coefficient selection
     */
    ts_type **dft_mem_array_coeff = NULL;

    if (use_variance) {
        dft_mem_array_coeff =
                calculate_variance_coeff(index, dft_mem_array);

        free_dft_memory(index, n_coefficients, dft_mem_array);
        if (dft_mem_array_coeff == NULL) {
            free(input_data);
            return FAILURE;
        }
    }

    ts_type **root_split_coefficients =
            use_variance ? dft_mem_array_coeff : dft_mem_array;
    double variance_mean = 0.0;
    for (int i = 0; i < n_segments; ++i) variance_mean += index->settings->symbolic_variances[i];
    variance_mean /= (double) n_segments;
    if (variance_mean > 0.0) {
        for (int i = 0; i < n_segments; ++i)
            index->settings->symbolic_variances[i] /= variance_mean;
    }
    const int root_budget = index->settings->index_type == MESSI_INDEX_TRIE &&
                                    index->settings->trie_dynamic_alphabet
                                ? index->settings->trie_alphabet_budget_bits * n_segments
                                : (n_segments < (int) (sizeof(root_mask_type) * 8)
                                       ? n_segments : (int) (sizeof(root_mask_type) * 8));
    if (configure_dynamic_bit_allocation(index, index->settings->symbolic_variances,
                                           n_segments, root_budget,
                                           index->settings->index_type == MESSI_INDEX_TRIE &&
                                                   index->settings->trie_dynamic_alphabet
                                               ? index->settings->trie_min_bits : 0,
                                           index->settings->index_type == MESSI_INDEX_TRIE &&
                                                   index->settings->trie_dynamic_alphabet
                                               ? index->settings->trie_max_bits
                                               : index->settings->sax_bit_cardinality) != SUCCESS) {
        free(input_data);
        free_dft_memory(index, n_segments, root_split_coefficients);
        return FAILURE;
    }
    double selection_end = messi_monotonic_seconds();

    /*
     * Phase 2:
     * divide coefficient rows evenly among workers
     */
    const int base = n_segments / worker_threads;
    const int remainder = n_segments % worker_threads;

    for (int i = 0; i < worker_threads; ++i) {
        bins_data_inmemory *data = &input_data[i];

        fftw_workspace_destroy(&data->fftw);

        const int extra_before = i < remainder ? i : remainder;
        const int work = base + (i < remainder ? 1 : 0);

        data->start_number = i * base + extra_before;
        data->stop_number = data->start_number + work;

        if (use_variance) {
            data->dft_mem_array = dft_mem_array_coeff;
        }
    }

    for (int i = 0; i < worker_threads; ++i) {
        pthread_create(
            &threadid[i],
            NULL,
            order_divide_worker,
            &input_data[i]);
    }

    for (int i = 0; i < worker_threads; ++i) {
        pthread_join(threadid[i], NULL);
    }
    double bins_end = messi_monotonic_seconds();

    free(input_data);

    free_dft_memory(index, n_segments, root_split_coefficients);

    COUNT_BINNING_TIME_END

    sfa_print_bins(index);
    fprintf(stderr, ">>> SFA binning timing\n");
    fprintf(stderr, "    sample + DFT          : %.3f s\n",
            sampling_end - binning_start);
    fprintf(stderr, "    variance + root split : %.3f s\n",
            selection_end - sampling_end);
    fprintf(stderr, "    bin sorting           : %.3f s\n",
            bins_end - selection_end);
    fprintf(stderr, "    total                 : %.3f s\n",
            bins_end - binning_start);
    fprintf(stderr, ">>> Finished binning\n");
    return SUCCESS;
}

/*
  This function calculates the variance for each coefficient in dft_mem_array.
  It returns a trimmed dft_mem_array with only the highest-variance coeff.
*/
ts_type **calculate_variance_coeff(isax_index *index, ts_type **dft_mem_array) {
    int n_coefficients = index->settings->n_coefficients;
    int n_segments = index->settings->n_segments;
    unsigned int sample_size = index->settings->sample_size;

    struct variance_coeff_index var_coeff_index[n_coefficients / 2];

    for (int i = 0; i < n_coefficients / 2; ++i) {
        double mean_real = 0.0;
        double mean_imag = 0.0;
        double var_real = 0.0;
        double var_imag = 0.0;

        for (int j = 0; j < sample_size; ++j) {
            mean_real += dft_mem_array[i * 2][j];
            mean_imag += dft_mem_array[i * 2 + 1][j];
        }
        mean_real = mean_real / (double) sample_size;
        mean_imag = mean_imag / (double) sample_size;

        for (int j = 0; j < sample_size; ++j) {
            var_real += (dft_mem_array[i * 2][j] - mean_real) * (
                dft_mem_array[i * 2][j] - mean_real);
            var_imag += (dft_mem_array[i * 2 + 1][j] - mean_imag) * (
                dft_mem_array[i * 2 + 1][j] - mean_imag);
        }
        var_real = var_real / (double) sample_size;
        var_imag = var_imag / (double) sample_size;

        double total_var = var_real + var_imag;

        var_coeff_index[i].variance = total_var;
        var_coeff_index[i].variance_real = var_real;
        var_coeff_index[i].variance_imag = var_imag;
        var_coeff_index[i].coeff_index = i;
    }

    /*
    fprintf(stderr, "Variance: ");
    for (int i = 0; i < n_coefficients; ++i) {
        // fprintf(stderr, "%.3f\tposition %d\n", var_coeff_index[i].variance, var_coeff_index[i].coeff_index);
        fprintf(stderr, "%.4f, ", var_coeff_index[i].variance);
    }
    fprintf(stderr, "\n");
    */

    qsort(var_coeff_index, n_coefficients / 2, sizeof(var_coeff_index[0]), compare_var);

    const int candidate_complex = n_coefficients / 2;
    const int print_candidates = candidate_complex < 8 ? candidate_complex : 8;
    fprintf(stderr, ">>> SFA: variance ranking (top %d/%d):\n",
            print_candidates, candidate_complex);
    for (int i = 0; i < print_candidates; ++i) {
        fprintf(stderr, "%d, (%.4f) ", var_coeff_index[i].coeff_index,
                var_coeff_index[i].variance);
    }
    if (print_candidates < candidate_complex) {
        fprintf(stderr, "... (+%d more)", candidate_complex - print_candidates);
    }
    fprintf(stderr, "\n");

    for (int i = 0; i < n_segments / 2; ++i) {
        index->coefficients[i] = var_coeff_index[i].coeff_index;
    }

    /*
     * iSAX uses all transformed dimensions as its bound, so retain the
     * historical coefficient-index order there.  Trie, however, keeps extra
     * dimensions for MBRs/splitting and uses only a prefix for record bounds.
     * Put the strongest coefficients at the front of that prefix; otherwise
     * truncating the numerically sorted list silently discards the variance
     * ranking.
     */
    const int selected_complex = n_segments / 2;
    const int record_complex =
            (index->settings->index_type == MESSI_INDEX_TRIE &&
             index->settings->trie_bound_dimensions > 0)
                ? index->settings->trie_bound_dimensions / 2
                : selected_complex;

    if (index->settings->index_type != MESSI_INDEX_TRIE) {
        qsort(index->coefficients, selected_complex, sizeof(int), compare_int);
    } else {
        /* The ranking array is already descending by variance.  The first
         * record_complex entries become the trie record-bound prefix and the
         * remaining selected entries are retained for MBR dimensions. */
        int *ordered = calloc((size_t) selected_complex, sizeof(*ordered));
        if (ordered == NULL) {
            return NULL;
        }
        for (int i = 0; i < selected_complex; ++i) {
            ordered[i] = var_coeff_index[i].coeff_index;
        }
        memcpy(index->coefficients, ordered,
               sizeof(*ordered) * (size_t) selected_complex);
        free(ordered);
    }

    fprintf(stderr, ">>> SFA: Coefficients Used: ");
    for (int i = 0; i < selected_complex; ++i) {
        fprintf(stderr, "%d, ", index->coefficients[i]);
    }
    fprintf(stderr, "\n");

    if (index->settings->index_type == MESSI_INDEX_TRIE) {
        fprintf(stderr, ">>> SFA: Trie record-bound coeffs: ");
        for (int i = 0; i < record_complex && i < selected_complex; ++i) {
            fprintf(stderr, "%d, ", index->coefficients[i]);
        }
        fprintf(stderr, "\n");
    }

    /* Keep the training variances aligned with the selected representation.
     * Dynamic alphabet allocation consumes these values directly. */
    double *selected_variance = calloc((size_t) n_segments, sizeof(*selected_variance));
    if (selected_variance == NULL) return NULL;
    for (int i = 0; i < n_segments / 2; ++i) {
        selected_variance[i * 2] = var_coeff_index[i].variance_real;
        selected_variance[i * 2 + 1] = var_coeff_index[i].variance_imag;
    }
    free(index->settings->symbolic_variances);
    index->settings->symbolic_variances = selected_variance;

    ts_type **dft_mem_array_coeff = (ts_type **) calloc(n_segments, sizeof(ts_type *));
    for (int k = 0; k < n_segments; ++k) {
        dft_mem_array_coeff[k] = (ts_type *) calloc(sample_size, sizeof(ts_type));
    }

    for (int i = 0; i < n_segments / 2; ++i) {
        int coeff = index->coefficients[i];

        memcpy(dft_mem_array_coeff[i * 2],
               dft_mem_array[coeff * 2],
               sizeof(ts_type) * sample_size);
        memcpy(dft_mem_array_coeff[i * 2 + 1],
               dft_mem_array[coeff * 2 + 1],
               sizeof(ts_type) * sample_size);
    }

    return dft_mem_array_coeff;
}

/*
    Worker method for sampling values, calculating FFT coefficients (the first n_coefficients coefficients) and saving them to dft_mem_array
*/
void *set_bins_worker_dft(void *transferdata) {
    struct bins_data_inmemory *bins_data = (bins_data_inmemory *) transferdata;

    ts_type **dft_mem_array = bins_data->dft_mem_array;

    isax_index *index = ((bins_data_inmemory *) transferdata)->index;
    unsigned long start_number = bins_data->start_number;
    unsigned long stop_number = bins_data->stop_number;
    uint64_t rng_state = ((uint64_t) index->settings->sampling_seed << 32) ^
                         (uint64_t) (bins_data->workernumber + 1);

    unsigned long ts_length = index->settings->timeseries_size;
    long records = bins_data->records;
    if (records <= 0) {
        bins_data->status = FAILURE;
        return NULL;
    }

    int n_coefficients = 0;
    // Variance based coefficients
    if (index->settings->n_coefficients > 0) {
        n_coefficients = index->settings->n_coefficients;
    }
    // first coefficients
    else {
        n_coefficients = index->settings->n_segments;
    }
    if (index->settings->function_type == 6) {
        n_coefficients = (int) ts_length;
    }

    unsigned long skip_elements = 0;
    if (index->settings->sample_type == 2) {
        unsigned long span = stop_number - start_number;
        if (stop_number <= start_number || (unsigned long) records > span) {
            bins_data->status = FAILURE;
            return NULL;
        }
        skip_elements = (span / (unsigned long) records) - 1;
    }
    if (index->settings->sample_type == 3 && stop_number <= start_number) {
        bins_data->status = FAILURE;
        return NULL;
    }

    unsigned long start_index = start_number * ts_length * sizeof(ts_type);
    int filetype_int = bins_data->filetype_int;
    int apply_znorm = bins_data->apply_znorm;

    FILE *ifile;
    ifile = fopen(bins_data->filename, "rb");
    if (ifile == NULL) {
        bins_data->status = FAILURE;
        return NULL;
    }
    fseek(ifile, start_index, SEEK_SET);

    unsigned long position_count = start_number;


    ts_type *ts = bins_data->fftw.ts;
    fftw_workspace *fftw = &bins_data->fftw;

    file_type *ts_orig1 = NULL;
    ts_type *ts_orig2 = NULL;

    if (filetype_int) {
        ts_orig1 = (file_type *) calloc(index->settings->timeseries_size,
                                        sizeof(file_type));
    } else {
        ts_orig2 = (ts_type *)
                calloc(index->settings->timeseries_size, sizeof(ts_type));
    }
    if ((filetype_int && ts_orig1 == NULL) || (!filetype_int && ts_orig2 == NULL)) {
        bins_data->status = FAILURE;
        free(ts_orig1);
        free(ts_orig2);
        fclose(ifile);
        return NULL;
    }

    for (int i = 0; i < records; ++i) {
        //choose random position for random sampling
        if (index->settings->sample_type == 3) {
            unsigned long span = stop_number - start_number;
            unsigned long position = start_number +
                                     (unsigned long) random_at_most_seed(&rng_state,
                                                                         (long int) span - 1);
            fseek(ifile, (position * ts_length * sizeof(ts_type)), SEEK_SET);
        }

        if (filetype_int) {
            if (fread(ts_orig1, sizeof(file_type), ts_length, ifile) != ts_length) {
                bins_data->status = FAILURE;
                break;
            }
            for (int j = 0; j < ts_length; ++j) {
                ts[j] = (ts_type) ts_orig1[j];
            }
        } else {
            if (fread(ts_orig2, sizeof(ts_type), ts_length, ifile) != ts_length) {
                bins_data->status = FAILURE;
                break;
            }
            for (int j = 0; j < ts_length; ++j) {
                ts[j] = ts_orig2[j];
            }
        }
        // apply z-normalization
        if (apply_znorm) {
            znorm(ts, ts_length);
        }

        enum response transform_status;
        if (index->settings->function_type == 6) {
            transform_status = fft_full_real_from_ts(index, fftw);
        } else {
            int transform_size = index->settings->n_coefficients != 0
                                     ? index->settings->n_coefficients
                                     : index->settings->n_segments;
            transform_status = fft_from_ts(index, transform_size, 0, fftw);
        }
        if (transform_status != SUCCESS) {
            bins_data->status = FAILURE;
            free(ts_orig1);
            free(ts_orig2);
            fclose(ifile);
            return NULL;
        }

        for (int j = 0; j < n_coefficients; ++j) {
            ts_type value = fftw->transform[j];
            dft_mem_array[j][i + (bins_data->workernumber * bins_data->records_offset)]
                    = value;
        }

        // skip elements for uniform sampling
        if (index->settings->sample_type == 2) {
            fseek(ifile, skip_elements * ts_length * sizeof(ts_type), SEEK_CUR);
            position_count += (1 + skip_elements);
        }

        /*
            // TODO which one is correct? above or this?
            // skip elements for uniform sampling
            if (index->settings->sample_type == 2) {
                fseek(ifile, skip_elements, SEEK_CUR);
            }
         */
    }

    free(ts_orig1);
    free(ts_orig2);
    fclose(ifile);
    return NULL;
}

/*
    Worker method for coefficient-wise splitting and saving the interval to bins
*/
void *order_divide_worker(void *transferdata) {
    struct bins_data_inmemory *bins_data = (bins_data_inmemory *) transferdata;

    ts_type **dft_mem_array = bins_data->dft_mem_array;

    isax_index *index = ((bins_data_inmemory *) transferdata)->index;
    unsigned long start_number = bins_data->start_number;
    unsigned long stop_number = bins_data->stop_number;

    unsigned int sample_size = index->settings->sample_size;
    int n_segments = index->settings->n_segments;
    ts_type *cur_coeff_line;

    for (int j = start_number; j < stop_number; ++j) {
        cur_coeff_line = (ts_type *) dft_mem_array[j];
        qsort(cur_coeff_line, sample_size, sizeof(ts_type), &compare_ts_type);
    }

    // equi-depth splitting
    if (index->settings->histogram_type == 1) {
        int num_symbols = index->settings->sax_alphabet_cardinality;
        ts_type depth = (ts_type) sample_size / num_symbols;

        for (int i = start_number; i < stop_number; ++i) {
            float bin_index = 0.0;
            cur_coeff_line = dft_mem_array[i];
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
            cur_coeff_line = dft_mem_array[i];
            ts_type first = cur_coeff_line[0];
            ts_type last = cur_coeff_line[sample_size - 1];
            ts_type interval_width = (last - first) / (ts_type) num_symbols;
            for (int j = 0; j < num_symbols - 1; ++j) {
                index->bins[i][j] = interval_width * (j + 1) + first;
            }
        }
    }

    if (n_segments == 0) {
        fprintf(stderr, "warning: SFA has zero segments.\n");
    }
    return NULL;
}


/*
    Method for printing the discretization intervals to the console
*/
void sfa_print_bins(isax_index *index) {
    fprintf(stderr, ">>> SFA: Sample size %u\n", index->settings->sample_size);
    if (index->settings->histogram_type == 1) {
        fprintf(stderr, ">>> SFA: Using Equi-depth histograms\n");
    } else if (index->settings->histogram_type == 2) {
        fprintf(stderr, ">>> SFA: Using Equi-width histograms\n");
    }

    /*
    int n_segments = index->settings->n_segments;
    fprintf(stderr,"[\n");
    for (int i = 0; i < n_segments; ++i)
    {
        fprintf(stderr,"-Inf\t");
        for (int j=0; j < index->settings->sax_alphabet_cardinality-1; ++j)
        {
            ts_type value = roundf(index->bins[i][j]*100.0)/100.0;
            if (value == FLT_MAX)
	            fprintf(stderr,",Inf\n");
	        else
	            fprintf(stderr,",%g",value);
        }
        fprintf(stderr,";\n");
    }
    fprintf(stderr,"]\n");
    */
}

void free_dft_memory(isax_index *index, int n_coefficients, ts_type **dft_mem_array) {
    // int n_segments = index->settings->n_segments;
    for (int k = 0; k < n_coefficients; ++k) {
        free(dft_mem_array[k]);
    }
    free(dft_mem_array);
}

//compare-functions for qsort
int compare_ts_type(const void *a, const void *b) {
    ts_type ts_a = *((ts_type *) a);
    ts_type ts_b = *((ts_type *) b);

    if (ts_a < ts_b)
        return -1;
    else return ts_a > ts_b;
}

int compare_var(const void *a, const void *b) {
    struct variance_coeff_index *a1 = (struct variance_coeff_index *) a;
    struct variance_coeff_index *a2 = (struct variance_coeff_index *) b;
    if ((*a1).variance > (*a2).variance)
        return -1;
    else if ((*a1).variance<(*a2).variance)
        return 1;
    else
        return 0;
}

int compare_int(const void *a, const void *b) {
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;
    return *ia - *ib;
}

void sfa_printbin(unsigned long long n, int size) {
    char *b = malloc(sizeof(char) * (size + 1));
    int i;

    for (i = 0; i < size; i++) {
        b[i] = '0';
    }

    for (i = 0; i < size; i++, n = n / 2)
        if (n % 2) b[size - 1 - i] = '1';

    b[size] = '\0';
    printf("%s\n", b);
    free(b);
}

void sfa_print(sax_type *sax, int segments, int cardinality) {
    int i;
    for (i = 0; i < segments; i++) {
        printf("%d:\t\n", i);
        sfa_printbin(sax[i], cardinality);
    }
    printf("\n");
}

void fft_print(ts_type *fft, int segments) {
    int i;
    for (i = 0; i < segments; i++) {
        printf("%d:\t%.3f\n", i, fft[i]);
    }
    printf("\n");
}

/*
    This function calculates random numbers between 0 and max
*/
long random_at_most(long max) {
    unsigned long
            num_bins = (unsigned long) max + 1,
            num_rand = (unsigned long) RAND_MAX + 1,
            bin_size = num_rand / num_bins,
            defect = num_rand % num_bins;

    long x;
    do {
        x = random();
    } while (num_rand - defect <= (unsigned long) x);

    return x / bin_size;
}

/* Deterministic, worker-local generator used by SFA/PISA bin sampling. */
static uint64_t sfa_rng_next(uint64_t *state) {
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

long random_at_most_seed(uint64_t *state, long max) {
    if (max <= 0) return 0;
    const uint64_t range = (uint64_t) max + 1;
    const uint64_t limit = UINT64_MAX - (UINT64_MAX % range);
    uint64_t value;
    do {
        value = sfa_rng_next(state);
    } while (value >= limit);
    return (long) (value % range);
}
