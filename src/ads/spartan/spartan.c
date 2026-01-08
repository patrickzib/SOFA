#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "globals.h"
#include "ads/calc_utils.h"
#include "ads/spartan/pca.h"
#include "ads/sax/sax.h"
#include "ads/sax/ts.h"
#include "ads/spartan/spartan.h"
#include "ads/sfa/sfa.h"
#if ADS_HAVE_AVX2
#include "immintrin.h"
#endif

typedef struct spartan_bins_data {
    isax_index *index;
    ts_type **coeff_mem_array;
    long int start_number;
    long int stop_number;
} spartan_bins_data;

static long spartan_random_at_most(long max) {
    unsigned long num_bins = (unsigned long) max + 1;
    unsigned long num_rand = (unsigned long) RAND_MAX + 1;
    unsigned long bin_size = num_rand / num_bins;
    unsigned long defect = num_rand % num_bins;

    long x;
    do {
        x = random();
    } while (num_rand - defect <= (unsigned long) x);

    return x / bin_size;
}

enum response spartan_bins_init(isax_index *index) {
    int num_symbols = index->settings->sax_alphabet_cardinality;
    int n_segments = index->settings->n_segments;

    index->bins = NULL;
    index->bins = (ts_type **) calloc(n_segments, sizeof(ts_type *));
    index->binsv = (ts_type *) calloc(n_segments * (num_symbols - 1), sizeof(ts_type));

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
    ts_type *ts = malloc(sizeof(ts_type) * ts_length);
    file_type *ts_orig1 = NULL;
    ts_type *ts_orig2 = NULL;
    if (filetype_int) {
        ts_orig1 = (file_type *) calloc(ts_length, sizeof(file_type));
    } else {
        ts_orig2 = (ts_type *) calloc(ts_length, sizeof(ts_type));
    }

    unsigned int records = sample_size;
    if ((long int) records > ts_num) {
        records = (unsigned int) ts_num;
    }

    unsigned long start_number = 0;
    unsigned long stop_number = ts_num;
    if (index->settings->sample_type == 1) {
        stop_number = records;
    }

    unsigned long skip_elements = 0;
    if (index->settings->sample_type == 2 && records > 0) {
        skip_elements = ((stop_number - start_number) / records);
        if (skip_elements > 0) {
            skip_elements -= 1;
        }
    }

    for (unsigned int i = 0; i < records; ++i) {
        if (index->settings->sample_type == 3) {
            long int position = spartan_random_at_most((long) ts_num - 1);
            fseek(ifile, (position * ts_length * sizeof(ts_type)), SEEK_SET);
        }

        if (filetype_int) {
            fread(ts_orig1, sizeof(file_type), ts_length, ifile);
            for (unsigned long j = 0; j < ts_length; ++j) {
                ts[j] = (ts_type) ts_orig1[j];
            }
        } else {
            fread(ts_orig2, sizeof(ts_type), ts_length, ifile);
            memcpy(ts, ts_orig2, sizeof(ts_type) * ts_length);
        }

        if (apply_znorm) {
            znorm(ts, ts_length);
        }

        memcpy(samples + (i * ts_length), ts, sizeof(ts_type) * ts_length);

        if (index->settings->sample_type == 2 && skip_elements > 0) {
            fseek(ifile, skip_elements * ts_length * sizeof(ts_type), SEEK_CUR);
        }
    }

    fclose(ifile);
    free(ts);
    free(ts_orig1);
    free(ts_orig2);

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
    } else if (index->settings->histogram_type == 2) {
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

void spartan_set_bins(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread,
                      int filetype_int, int apply_znorm) {
    int dim = index->settings->n_segments;
    int ts_length = index->settings->timeseries_size;
    unsigned int sample_size = index->settings->sample_size;

    fprintf(stderr, ">>> SPARTAN binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START

    ts_type *samples = calloc((size_t) sample_size * ts_length, sizeof(ts_type));
    if (samples == NULL) {
        fprintf(stderr, "error: failed to allocate SPARTAN sample buffer.\n");
        return;
    }

    if (spartan_collect_samples(index, ifilename, ts_num, filetype_int, apply_znorm, samples, sample_size) != SUCCESS) {
        free(samples);
        return;
    }

    if (pca_fit(index, samples, sample_size, ts_length) != SUCCESS) {
        free(samples);
        return;
    }

    ts_type **coeff_mem_array = (ts_type **) calloc(dim, sizeof(ts_type *));
    for (int k = 0; k < dim; ++k) {
        coeff_mem_array[k] = (ts_type *) calloc(sample_size, sizeof(ts_type));
    }

    ts_type *projection = calloc(dim, sizeof(ts_type));
    for (unsigned int i = 0; i < sample_size; ++i) {
        const ts_type *row = samples + (i * ts_length);
        pca_from_ts(index, row, projection);
        for (int k = 0; k < dim; ++k) {
            coeff_mem_array[k][i] = projection[k];
        }
    }
    free(projection);
    free(samples);

    pthread_t threadid[maxquerythread];
    spartan_bins_data *input_data = malloc(sizeof(spartan_bins_data) * (size_t) maxquerythread);

    for (int i = 0; i < maxquerythread; i++) {
        input_data[i].index = index;
        input_data[i].coeff_mem_array = coeff_mem_array;
        input_data[i].start_number = i * (dim / maxquerythread);
        input_data[i].stop_number = (i + 1) * (dim / maxquerythread);
    }

    input_data[maxquerythread - 1].start_number = (maxquerythread - 1) * (dim / maxquerythread);
    input_data[maxquerythread - 1].stop_number = dim;

    for (int i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, spartan_order_divide_worker, (void *) &(input_data[i]));
    }

    for (int i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }

    for (int k = 0; k < dim; ++k) {
        free(coeff_mem_array[k]);
    }
    free(coeff_mem_array);
    free(input_data);

    COUNT_BINNING_TIME_END

    spartan_print_bins(index);
    fprintf(stderr, ">>> Finished SPARTAN binning\n");
}

void spartan_from_pca(isax_index *index, const ts_type *coeffs, sax_type *sax_out) {
    int n_segments = index->settings->n_segments;
    for (int k = 0; k < n_segments; ++k) {
        unsigned int c;
        for (c = 0; c < index->settings->sax_alphabet_cardinality - 1; c++) {
            if (coeffs[k] < index->bins[k][c]) {
                break;
            }
        }
        sax_out[k] = (unsigned char) (c);
    }
}

enum response spartan_from_ts(isax_index *index, const ts_type *ts, sax_type *sax_out) {
    int dim = index->settings->n_segments;
    ts_type *coeffs = calloc(dim, sizeof(ts_type));
    if (coeffs == NULL) {
        free(coeffs);
        return FAILURE;
    }

    if (pca_from_ts(index, ts, coeffs) != SUCCESS) {
        free(coeffs);
        return FAILURE;
    }
    spartan_from_pca(index, coeffs, sax_out);

    free(coeffs);

    if (sax_out != NULL) {
        return SUCCESS;
    }
    fprintf(stderr, "SPARTAN error\n");
    return FAILURE;
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
            double centered = (double) ts[i] - index->pca_mean[i];
            acc += centered * component[i];
        }
        out[k] = (ts_type) acc;
    }
    for (int k = components; k < output_dim; ++k) {
        out[k] = 0.0f;
    }
    return SUCCESS;
}

ts_type minidist_pca_to_spartan(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities, float bsf) {
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->n_segments;

    ts_type distance = 0.0;
    for (int i = 0; i < number_of_segments; i++) {
        distance += get_lb_distance(
                index->bins[i], pca[i], sax[i], sax_cardinalities[i],
                max_bit_cardinality, max_cardinality, 1.0);

        if (distance > bsf) {
            return distance;
        }
    }
    return distance;
}

ts_type minidist_pca_to_spartan_raw(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities,
                                    float bsf) {
    return minidist_pca_to_spartan(index, pca, sax, sax_cardinalities, bsf);
}

#if ADS_HAVE_AVX2
ts_type
minidist_pca_to_spartan_rawe_SIMD(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities, float bsf) {

    int region_upper[16], region_lower[16];
    float distancef[8], distancef2[8];
    int offset = 0;
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;

    __m256i vectorsignbit = _mm256_set1_epi32(0xffffffff);

    __m128i sax_cardinalitiesv8 = _mm_lddqu_si128((const void *) sax_cardinalities);
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

    paav_0 = _mm256_loadu_ps(pca);

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
    if (distancef[0] + distancef[4] > bsf) {
        return distancef[0] + distancef[4];
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
    paav_1 = _mm256_loadu_ps(&(pca[8]));
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

    return distancef[0] + distancef[4] + distancef2[0] + distancef2[4];
}
#endif

void spartan_print_bins(isax_index *index) {
    fprintf(stderr, ">>> SPARTAN: Sample size %u\n", index->settings->sample_size);
    if (index->settings->histogram_type == 1) {
        fprintf(stderr, ">>> SPARTAN: Using Equi-depth histograms\n");
    } else if (index->settings->histogram_type == 2) {
        fprintf(stderr, ">>> SPARTAN: Using Equi-width histograms\n");
    }
}
