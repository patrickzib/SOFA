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
#if ADS_HAVE_AVX2
#include "immintrin.h"
#endif

#ifdef VALUES

#include <values.h>

#endif

#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include <unistd.h>

#include "ads/calc_utils.h"
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
    index->bins = (ts_type **) calloc(n_segments, sizeof(ts_type * ));
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
    fprintf(stderr, ">>> SFA: Initialized bins[%d][%d] \n", n_segments, num_symbols - 1);

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
  In this function, the intervals are caluclated (multiple coeff. binning).
  The coefficients with the highest variance or the first ones are chosen and
  these values are saved to bins
*/
void sfa_set_bins(
        isax_index *index, const char *ifilename,
        long int ts_num, int maxquerythread,
        int filetype_int, int apply_znorm) {

    int n_segments = index->settings->n_segments;
    int n_coefficients = 0;
    int ts_length = index->settings->timeseries_size;
    unsigned int sample_size = index->settings->sample_size;
    int use_variance = index->settings->n_coefficients > 0;

    // a) select best coefficients based no variance
    if (use_variance) {
        n_coefficients = index->settings->n_coefficients;
    }
        // b) use the first coefficients
    else {
        n_coefficients = index->settings->n_segments;
    }

    fprintf(stderr, ">>> Binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START

    ts_type **dft_mem_array = (ts_type **) calloc(n_coefficients, sizeof(ts_type * ));
    for (int k = 0; k < n_coefficients; ++k) {
        dft_mem_array[k] = (ts_type *) calloc(sample_size, sizeof(ts_type));
    }

    //build the bins out of a sample of the data
    //read whole sample in memory
    pthread_t threadid[maxquerythread];
    bins_data_inmemory *input_data = malloc(sizeof(bins_data_inmemory) * (maxquerythread));

    fftw_workspace fftw = {0};

    for (int i = 0; i < maxquerythread; i++) {
        //create FFTW-objects for the threads
        fftw_workspace_init(&fftw, ts_length);

        input_data[i].index = index;
        input_data[i].dft_mem_array = dft_mem_array;
        input_data[i].filename = ifilename;
        input_data[i].workernumber = i;
        input_data[i].records = sample_size / maxquerythread;
        input_data[i].records_offset = sample_size / maxquerythread;

        // first-n-sampling
        if (index->settings->sample_type == 1) {
            input_data[i].start_number = i * (sample_size / maxquerythread);
            input_data[i].stop_number = (i + 1) * (sample_size / maxquerythread);
        }
        // uniform sampling
        else if (index->settings->sample_type == 2) {
            input_data[i].start_number = i * (ts_num / maxquerythread);
            input_data[i].stop_number = (i + 1) * (ts_num / maxquerythread);
        }
        // random sampling
        else if (index->settings->sample_type == 3) {
            input_data[i].start_number = 0;
            input_data[i].stop_number = ts_num;
        }

        input_data[i].filetype_int = filetype_int;
        input_data[i].apply_znorm = apply_znorm;

        input_data[i].fftw = fftw;
    }

    // reset values for last worker to keep lost segments at the end
    input_data[maxquerythread - 1].records = sample_size - (maxquerythread - 1) * (sample_size / maxquerythread);

    if (index->settings->sample_type == 1) {
        input_data[maxquerythread - 1].stop_number = sample_size;
    } else if (index->settings->sample_type == 2) {
        input_data[maxquerythread - 1].stop_number = ts_num;
    }

    // initiate worker threads for sampling values and calculating FFTs
    for (int i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, set_bins_worker_dft, (void *) &(input_data[i]));
    }

    // wait for the finish of other threads
    for (int i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }

    ts_type **dft_mem_array_coeff;
    if (use_variance) {
        // calculate coefficient-wise variance
        dft_mem_array_coeff = calculate_variance_coeff(index, dft_mem_array);

        free_dft_memory(index, n_coefficients, dft_mem_array);
    }

    for (int i = 0; i < maxquerythread; i++) {
        fftw_workspace_destroy(&input_data[i].fftw);

        input_data[i].start_number = i * (n_segments / maxquerythread);
        input_data[i].stop_number = (i + 1) * (n_segments / maxquerythread);

        if (use_variance) {
            // replace with new ordering
            input_data[i].dft_mem_array = dft_mem_array_coeff;
        }
    }

    input_data[maxquerythread - 1].start_number = (maxquerythread - 1) * (n_segments / maxquerythread);
    input_data[maxquerythread - 1].stop_number = n_segments;

    //initiate worker threads for splitting coefficients into intervals
    for (int i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, order_divide_worker, (void *) &(input_data[i]));
    }
    //wait for the finish of other threads
    for (int i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }

    free(input_data);

    if (use_variance) {
        free_dft_memory(index, index->settings->n_segments, dft_mem_array_coeff);
    } else {
        free_dft_memory(index, index->settings->n_segments, dft_mem_array);
    }

    COUNT_BINNING_TIME_END

    sfa_print_bins(index);
    fprintf(stderr, ">>> Finished binning\n");
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
            var_real += (dft_mem_array[i * 2][j] - mean_real) * (dft_mem_array[i * 2][j] - mean_real);
            var_imag += (dft_mem_array[i * 2 + 1][j] - mean_imag) * (dft_mem_array[i * 2 + 1][j] - mean_imag);
        }
        var_real = var_real / (double) sample_size;
        var_imag = var_imag / (double) sample_size;

        double total_var = var_real + var_imag;

        var_coeff_index[i].variance = total_var;
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

    fprintf(stderr, ">>> SFA: Best Indices Sorted:\n");
    for (int i = 0; i < n_coefficients / 2; ++i) {
        fprintf(stderr, "%d, (%.4f) ", var_coeff_index[i].coeff_index, var_coeff_index[i].variance);
    }
    fprintf(stderr, "\n");

    for (int i = 0; i < n_segments / 2; ++i) {
        index->coefficients[i] = var_coeff_index[i].coeff_index;
    }

    // sorting needed?
    qsort(index->coefficients, n_segments / 2, sizeof(int), compare_int);
    fprintf(stderr, ">>> SFA: Hightest Variance Coeffs Sorted: ");
    for (int i = 0; i < n_segments / 2; ++i) {
        fprintf(stderr, "%d, ", index->coefficients[i]);
    }
    fprintf(stderr, "\n");

    ts_type **dft_mem_array_coeff = (ts_type **) calloc(n_segments, sizeof(ts_type * ));
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

    unsigned long ts_length = index->settings->timeseries_size;

    int n_coefficients = 0;
    // Variance based coefficients
    if (index->settings->n_coefficients > 0) {
        n_coefficients = index->settings->n_coefficients;
    }
        // first coefficients
    else {
        n_coefficients = index->settings->n_segments;
    }

    unsigned long start_index = start_number * ts_length * sizeof(ts_type);
    int filetype_int = bins_data->filetype_int;
    int apply_znorm = bins_data->apply_znorm;

    FILE *ifile;
    ifile = fopen(bins_data->filename, "rb");
    fseek(ifile, start_index, SEEK_SET);

    unsigned long skip_elements;
    long records = bins_data->records;

    unsigned long position_count = start_number;

    //set number of elements to skip for uniform sampling
    if (index->settings->sample_type == 2) {
        skip_elements = (((stop_number - start_number) / records) - 1);

        // TODO skip_elements = (((stop_number - start_number) / records) - 1) * ts_length * sizeof(ts_type);
    }

    ts_type *ts = bins_data->fftw.ts;
    fftw_workspace *fftw = &bins_data->fftw;

    file_type *ts_orig1 = NULL;
    ts_type *ts_orig2 = NULL;

    if (filetype_int) {
        ts_orig1 = (file_type *) calloc(index->settings->timeseries_size, sizeof(file_type));
    } else {
        ts_orig2 = (ts_type *) calloc(index->settings->timeseries_size, sizeof(ts_type));
    }

    for (int i = 0; i < records; ++i) {
        //choose random position for random sampling
        if (index->settings->sample_type == 3) {
            long int position = start_number + random_at_most(records);
            fseek(ifile, (position * ts_length * sizeof(ts_type)), SEEK_SET);
        }

        if (filetype_int) {
            fread(ts_orig1, sizeof(file_type), ts_length, ifile);
            for (int j = 0; j < ts_length; ++j) {
                ts[j] = (ts_type) ts_orig1[j];
            }
        } else {
            fread(ts_orig2, sizeof(ts_type), ts_length, ifile);
            for (int j = 0; j < ts_length; ++j) {
                ts[j] = ts_orig2[j];
            }
        }
        // apply z-normalization
        if (apply_znorm) {
            znorm(ts, ts_length);
        }

        int use_best = index->settings->n_coefficients != 0;
        if (use_best) {
            fft_from_ts(index, index->settings->n_coefficients, 0, fftw);
        } else {
            fft_from_ts(index, index->settings->n_segments, 0, fftw);
        }

        for (int j = 0; j < n_coefficients; ++j) {
            ts_type value = fftw->transform[j];
            dft_mem_array[j][i + (bins_data->workernumber * bins_data->records_offset)] = value;
        }

        // skip elements for uniform sampling
        if (index->settings->sample_type == 2) {
            fseek(ifile, skip_elements * ts_length * sizeof(ts_type), SEEK_CUR);
            position_count += (1 + skip_elements);
            if (position_count >= stop_number) {
                fprintf(stderr, "pos %lu; stop_number %lu\n", position_count, stop_number);
            }
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
    else if ((*a1).variance < (*a2).variance)
        return 1;
    else
        return 0;
}

int compare_int(const void *a, const void *b) {
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;
    return *ia - *ib;
}

ts_type
get_lb_distance(const ts_type *bins, const float fft, const sax_type v, const sax_type c_c,
                sax_type c_m, int max_cardinality, float factor) {
    sax_type region_lower = (v << (c_m - c_c));
    sax_type region_upper = (~((int) MAXFLOAT << (c_m - c_c)) | region_lower);

    ts_type distance = 0.0;
    float breakpoint_lower = 0.0;
    float breakpoint_upper = 0.0;

    if (region_lower == 0) {
        breakpoint_lower = MINVAL;
    } else {
        breakpoint_lower = bins[region_lower - 1];
    }
    if (region_upper == max_cardinality - 1) {
        breakpoint_upper = MAXVAL;
    } else {
        breakpoint_upper = bins[region_upper];
    }

    if (breakpoint_lower > fft) {
        ts_type value = breakpoint_lower - fft;
        distance += factor * value * value;
    } else if (breakpoint_upper < fft) {
        ts_type value = (fft - breakpoint_upper);
        distance += factor * value * value;
    }
    return distance;
}

/*
    This function calculates a mindist (lower bounding dist.) between a query (FFT coeff.) and a SFA representation
*/
ts_type minidist_fft_to_sfa(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf) {
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->n_segments;

    ts_type distance = 0.0;
    int i = 0;

    // Special case: for non-normalized series, treat the first coefficient specially.
    if (!index->settings->is_norm &&
        (index->settings->n_coefficients == 0 || index->coefficients[0] == 0)) {
        distance += get_lb_distance(index->bins[0], fft[0], sax[0], max_cardinality,
                                    max_bit_cardinality, max_cardinality, 1.0);
        if (distance > bsf) {
            return distance;
        }
        // Skip the imaginary part of the first coefficient when no variance-based selection is used.
        i = (index->settings->n_coefficients == 0) ? 2 : 1;
    }

    for (; i < number_of_segments; ++i) {
        distance += get_lb_distance(
                index->bins[i], fft[i], sax[i], max_cardinality,
                max_bit_cardinality, max_cardinality, 2.0);

        if (distance > bsf) {
            return distance;
        }
    }

    return distance;
}

/*
    This function calculates a mindist (lower bounding dist.) between a query (FFT coeff.) and a SFA representation
*/
ts_type minidist_fft_to_sfa_raw(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf) {
    return minidist_fft_to_sfa(index, fft, sax, sax_cardinalities, bsf);
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


#if ADS_HAVE_AVX2
ts_type
minidist_fft_to_sfa_rawe_SIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf) {

    int region_upper[16], region_lower[16];
    float distancef[8], distancef2[8];
    int offset = 0;
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;refac

    __m256i vectorsignbit = _mm256_set1_epi32(0xffffffff);

    __m128i sax_cardinalitiesv8 = _mm_lddqu_si128((const void *) index->settings->max_sax_cardinalities);
    __m256i sax_cardinalitiesv16 = _mm256_cvtepu8_epi16(sax_cardinalitiesv8);
    __m128i sax_cardinalitiesv16_0 = _mm256_extractf128_si256(sax_cardinalitiesv16, 0);
    __m256i c_cv_0 = _mm256_cvtepu16_epi32(sax_cardinalitiesv16_0);

    __m128i saxv8 = _mm_lddqu_si128((const void *) sax);
    __m256i saxv16 = _mm256_cvtepu8_epi16(saxv8);
    __m128i saxv16_0 = _mm256_extractf128_si256(saxv16, 0);

    __m256i v_0 = _mm256_cvtepu16_epi32(saxv16_0);


    __m256i c_m = _mm256_set1_epi32(max_bit_cardinality);
    __m256i cm_ccv_0 = _mm256_sub_epi32(c_m, c_cv_0);

    __m256i region_lowerv_0 = _mm256_srlv_epi32(v_0, cm_ccv_0);

    region_lowerv_0 = _mm256_sllv_epi32(region_lowerv_0, cm_ccv_0);


    __m256i v1 = _mm256_andnot_si256(_mm256_setzero_si256(), vectorsignbit);

    __m256i region_upperv_0 = _mm256_sllv_epi32(v1, cm_ccv_0);

    region_upperv_0 = _mm256_andnot_si256(region_upperv_0, vectorsignbit);
    region_upperv_0 = _mm256_or_si256(region_upperv_0, region_lowerv_0);

    //lower
    __m256i lower_juge_zerov_0 = _mm256_cmpeq_epi32(region_lowerv_0, _mm256_setzero_si256());


    __m256i lower_juge_nzerov_0 = _mm256_andnot_si256(lower_juge_zerov_0, vectorsignbit);

    __m256 minvalv = _mm256_set1_ps(MINVAL);
    __m256i bitsizev = _mm256_set1_epi16((short) index->settings->sax_alphabet_cardinality - 1);
    __m256i bit1v = _mm256_set1_epi32(1);
    __m256i offsetvs = _mm256_set_epi16(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

    __m256i vsssssoffsetvs2 = _mm256_mullo_epi16(bitsizev, offsetvs);
    __m128i offsetv0s = _mm256_extractf128_si256(vsssssoffsetvs2, 0);

    __m128i offsetv1s = _mm256_extractf128_si256(vsssssoffsetvs2, 1);
    __m256i offsetv0 = _mm256_cvtepu16_epi32(offsetv0s);

    __m256i region_lowerbinv0 = _mm256_add_epi32(offsetv0, region_lowerv_0);

    region_lowerbinv0 = _mm256_sub_epi32(region_lowerbinv0, bit1v);

    __m256 lsax_breakpoints_shiftv_0 = _mm256_i32gather_ps(index->binsv, region_lowerbinv0, 4);

    __m256 breakpoint_lowerv_0 = (__m256) _mm256_or_si256(_mm256_and_si256(lower_juge_zerov_0, (__m256i) minvalv),
                                                          _mm256_and_si256(lower_juge_nzerov_0,
                                                                           (__m256i) lsax_breakpoints_shiftv_0));


    //upper
    __m256i region_upperbinv0 = _mm256_add_epi32(offsetv0, region_upperv_0);

    __m256 usax_breakpoints_shiftv_0 = _mm256_i32gather_ps(index->binsv, region_upperbinv0, 4);

    __m256i upper_juge_maxv_0 = _mm256_cmpeq_epi32(region_upperv_0, _mm256_set1_epi32(index->settings->sax_alphabet_cardinality - 1));

    __m256i upper_juge_nmaxv_0 = _mm256_andnot_si256(upper_juge_maxv_0, vectorsignbit);

    __m256 breakpoint_upperv_0 = (__m256) _mm256_or_si256(
            _mm256_and_si256(upper_juge_maxv_0, (__m256i) _mm256_set1_ps(MAXVAL)),
            _mm256_and_si256(upper_juge_nmaxv_0, (__m256i) usax_breakpoints_shiftv_0));


    //dis
    __m256 paav_0, paav_1;

    paav_0 = _mm256_loadu_ps(fft);

    __m256 dis_juge_upv_0 = _mm256_cmp_ps(breakpoint_lowerv_0, paav_0, _CMP_GT_OS);

    __m256 dis_juge_lov_0 = (__m256) _mm256_and_si256((__m256i) _mm256_cmp_ps(breakpoint_lowerv_0, paav_0, _CMP_NGT_US),
                                                      (__m256i) _mm256_cmp_ps(breakpoint_upperv_0, paav_0, _CMP_LT_OS));

    __m256 dis_juge_elv_0 = (__m256) _mm256_andnot_si256(
            _mm256_or_si256((__m256i) dis_juge_upv_0, (__m256i) dis_juge_lov_0), vectorsignbit);

    __m256 dis_lowv_0 = _mm256_mul_ps(_mm256_sub_ps(breakpoint_lowerv_0, paav_0),
                                      _mm256_sub_ps(breakpoint_lowerv_0, paav_0));
    __m256 dis_uppv_0 = _mm256_mul_ps(_mm256_sub_ps(breakpoint_upperv_0, paav_0),
                                      _mm256_sub_ps(breakpoint_upperv_0, paav_0));


    __m256 distancev_0 = (__m256) _mm256_or_si256(
            _mm256_or_si256(_mm256_and_si256((__m256i) dis_juge_upv_0, (__m256i) dis_lowv_0),
                            _mm256_and_si256((__m256i) dis_juge_lov_0, (__m256i) dis_uppv_0)),
            _mm256_and_si256((__m256i) dis_juge_elv_0, (__m256i) _mm256_set1_ps(0.0)));

    __m256 distancev2 = _mm256_hadd_ps(distancev_0, distancev_0);
    __m256 distancevf = _mm256_hadd_ps(distancev2, distancev2);

    _mm256_storeu_ps(distancef, distancevf);
    if ((distancef[0] + distancef[4]) * 2 > bsf) {
        return (distancef[0] + distancef[4]) * 2;
    }

    __m128i sax_cardinalitiesv16_1 = _mm256_extractf128_si256(sax_cardinalitiesv16, 1);
    __m256i c_cv_1 = _mm256_cvtepu16_epi32(sax_cardinalitiesv16_1);
    __m128i saxv16_1 = _mm256_extractf128_si256(saxv16, 1);
    __m256i v_1 = _mm256_cvtepu16_epi32(saxv16_1);
    __m256i cm_ccv_1 = _mm256_sub_epi32(c_m, c_cv_1);
    __m256i region_lowerv_1 = _mm256_srlv_epi32(v_1, cm_ccv_1);
    region_lowerv_1 = _mm256_sllv_epi32(region_lowerv_1, cm_ccv_1);
    __m256i region_upperv_1 = _mm256_sllv_epi32(v1, cm_ccv_1);
    region_upperv_1 = _mm256_andnot_si256(region_upperv_1, vectorsignbit);
    region_upperv_1 = _mm256_or_si256(region_upperv_1, region_lowerv_1);
    __m256i lower_juge_zerov_1 = _mm256_cmpeq_epi32(region_lowerv_1, _mm256_setzero_si256());
    __m256i lower_juge_nzerov_1 = _mm256_andnot_si256(lower_juge_zerov_1, vectorsignbit);
    __m256i offsetv1 = _mm256_cvtepu16_epi32(offsetv1s);
    __m256i region_lowerbinv1 = _mm256_add_epi32(offsetv1, region_lowerv_1);
    region_lowerbinv1 = _mm256_sub_epi32(region_lowerbinv1, bit1v);
    __m256 lsax_breakpoints_shiftv_1 = _mm256_i32gather_ps(index->binsv, region_lowerbinv1, 4);
    __m256 breakpoint_lowerv_1 = (__m256) _mm256_or_si256(_mm256_and_si256(lower_juge_zerov_1, (__m256i) minvalv),
                                                          _mm256_and_si256(lower_juge_nzerov_1,
                                                                           (__m256i) lsax_breakpoints_shiftv_1));

    __m256i region_upperbinv1 = _mm256_add_epi32(offsetv1, region_upperv_1);
    __m256 usax_breakpoints_shiftv_1 = _mm256_i32gather_ps(index->binsv, region_upperbinv1, 4);
    __m256i upper_juge_maxv_1 = _mm256_cmpeq_epi32(region_upperv_1, _mm256_set1_epi32(index->settings->sax_alphabet_cardinality - 1));
    __m256i upper_juge_nmaxv_1 = _mm256_andnot_si256(upper_juge_maxv_1, vectorsignbit);
    __m256 breakpoint_upperv_1 = (__m256) _mm256_or_si256(
            _mm256_and_si256(upper_juge_maxv_1, (__m256i) _mm256_set1_ps(MAXVAL)),
            _mm256_and_si256(upper_juge_nmaxv_1, (__m256i) usax_breakpoints_shiftv_1));
    paav_1 = _mm256_loadu_ps(&(fft[8]));
    __m256 dis_juge_upv_1 = _mm256_cmp_ps(breakpoint_lowerv_1, paav_1, _CMP_GT_OS);


    __m256 dis_juge_lov_1 = (__m256) _mm256_and_si256((__m256i) _mm256_cmp_ps(breakpoint_lowerv_1, paav_1, _CMP_NGT_US),
                                                      (__m256i) _mm256_cmp_ps(breakpoint_upperv_1, paav_1, _CMP_LT_OS));

    __m256 dis_juge_elv_1 = (__m256) _mm256_andnot_si256(
            _mm256_or_si256((__m256i) dis_juge_upv_1, (__m256i) dis_juge_lov_1), vectorsignbit);

    __m256 dis_lowv_1 = _mm256_mul_ps(_mm256_sub_ps(breakpoint_lowerv_1, paav_1),
                                      _mm256_sub_ps(breakpoint_lowerv_1, paav_1));

    __m256 dis_uppv_1 = _mm256_mul_ps(_mm256_sub_ps(breakpoint_upperv_1, paav_1),
                                      _mm256_sub_ps(breakpoint_upperv_1, paav_1));

    __m256 distancev_1 = (__m256) _mm256_or_si256(
            _mm256_or_si256(_mm256_and_si256((__m256i) dis_juge_upv_1, (__m256i) dis_lowv_1),
                            _mm256_and_si256((__m256i) dis_juge_lov_1, (__m256i) dis_uppv_1)),
            _mm256_and_si256((__m256i) dis_juge_elv_1, (__m256i) _mm256_set1_ps(0.0)));


    distancev2 = _mm256_hadd_ps(distancev_1, distancev_1);
    distancevf = _mm256_hadd_ps(distancev2, distancev2);
    _mm256_storeu_ps(distancef2, distancevf);

    return (distancef[0] + distancef[4] + distancef2[0] + distancef2[4]) * 2;

}
#endif
