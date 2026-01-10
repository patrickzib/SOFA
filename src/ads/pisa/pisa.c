#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "globals.h"
#include "ads/pisa/pisa.h"
#include "ads/sfa/sfa.h"
#include "ads/spartan/pca.h"
#include "ads/spartan/spartan.h"

enum response pisa_bins_init(isax_index *index) {
    return sfa_bins_init(index);
}

void pisa_free_bins(isax_index *index) {
    if (index == NULL) {
        return;
    }
    sfa_free_bins(index);
    pca_free(index);
}

static int pisa_fft_dim(const isax_index *index) {
    if (index == NULL || index->settings == NULL) {
        return 0;
    }
    if (index->settings->n_coefficients > 0) {
        return index->settings->n_coefficients;
    }
    return index->settings->n_segments;
}

enum response pisa_pca_from_ts(isax_index *index, const ts_type *ts, ts_type *out, fftw_workspace *fftw) {
    if (index == NULL || ts == NULL || out == NULL || fftw == NULL) {
        return FAILURE;
    }
    int fft_dim = pisa_fft_dim(index);
    if (fft_dim <= 0) {
        return FAILURE;
    }
    memcpy(fftw->ts, ts, sizeof(ts_type) * index->settings->timeseries_size);
    fft_from_ts(index, fft_dim, 0, fftw);
    return pca_from_ts(index, fftw->transform, out);
}

enum response pisa_from_ts(isax_index *index, const ts_type *ts, sax_type *sax_out, fftw_workspace *fftw) {
    if (index == NULL || ts == NULL || sax_out == NULL || fftw == NULL) {
        return FAILURE;
    }
    int n_segments = index->settings->n_segments;
    ts_type *coeffs = calloc(n_segments, sizeof(ts_type));
    if (coeffs == NULL) {
        return FAILURE;
    }
    if (pisa_pca_from_ts(index, ts, coeffs, fftw) != SUCCESS) {
        free(coeffs);
        return FAILURE;
    }
    sfa_from_fft(index, coeffs, sax_out);
    free(coeffs);
    return SUCCESS;
}

void pisa_set_bins(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread,
                   int filetype_int, int apply_znorm) {
    if (index == NULL || index->settings == NULL) {
        return;
    }
    int n_segments = index->settings->n_segments;
    int ts_length = index->settings->timeseries_size;
    unsigned int sample_size = index->settings->sample_size;
    int fft_dim = pisa_fft_dim(index);
    if (fft_dim <= 0 || n_segments <= 0 || sample_size == 0) {
        fprintf(stderr, "warning: invalid PISA settings.\n");
        return;
    }

    fprintf(stderr, ">>> PISA binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START

    ts_type **dft_mem_array = (ts_type **) calloc(fft_dim, sizeof(ts_type *));
    if (dft_mem_array == NULL) {
        fprintf(stderr, "error: failed to allocate PISA FFT samples.\n");
        return;
    }
    for (int k = 0; k < fft_dim; ++k) {
        dft_mem_array[k] = (ts_type *) calloc(sample_size, sizeof(ts_type));
        if (dft_mem_array[k] == NULL) {
            fprintf(stderr, "error: failed to allocate PISA FFT sample buffer.\n");
            free_dft_memory(index, k, dft_mem_array);
            return;
        }
    }

    pthread_t threadid[maxquerythread];
    bins_data_inmemory *input_data = malloc(sizeof(bins_data_inmemory) * (size_t) maxquerythread);
    if (input_data == NULL) {
        free_dft_memory(index, fft_dim, dft_mem_array);
        return;
    }

    fftw_workspace fftw = {0};
    for (int i = 0; i < maxquerythread; i++) {
        fftw_workspace_init(&fftw, ts_length);

        input_data[i].index = index;
        input_data[i].dft_mem_array = dft_mem_array;
        input_data[i].filename = ifilename;
        input_data[i].workernumber = i;
        input_data[i].records = sample_size / maxquerythread;
        input_data[i].records_offset = sample_size / maxquerythread;

        if (index->settings->sample_type == 1) {
            input_data[i].start_number = i * (sample_size / maxquerythread);
            input_data[i].stop_number = (i + 1) * (sample_size / maxquerythread);
        } else if (index->settings->sample_type == 2) {
            input_data[i].start_number = i * (ts_num / maxquerythread);
            input_data[i].stop_number = (i + 1) * (ts_num / maxquerythread);
        } else if (index->settings->sample_type == 3) {
            input_data[i].start_number = 0;
            input_data[i].stop_number = ts_num;
        }

        input_data[i].filetype_int = filetype_int;
        input_data[i].apply_znorm = apply_znorm;
        input_data[i].fftw = fftw;
    }

    input_data[maxquerythread - 1].records =
        sample_size - (maxquerythread - 1) * (sample_size / maxquerythread);

    if (index->settings->sample_type == 1) {
        input_data[maxquerythread - 1].stop_number = sample_size;
    } else if (index->settings->sample_type == 2) {
        input_data[maxquerythread - 1].stop_number = ts_num;
    }

    for (int i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, set_bins_worker_dft, (void *) &(input_data[i]));
    }

    for (int i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }

    for (int i = 0; i < maxquerythread; i++) {
        fftw_workspace_destroy(&input_data[i].fftw);
    }

    ts_type *samples = calloc((size_t) sample_size * (size_t) fft_dim, sizeof(ts_type));
    if (samples == NULL) {
        fprintf(stderr, "error: failed to allocate PISA PCA samples.\n");
        free(input_data);
        free_dft_memory(index, fft_dim, dft_mem_array);
        return;
    }
    for (unsigned int i = 0; i < sample_size; ++i) {
        for (int k = 0; k < fft_dim; ++k) {
            samples[i * fft_dim + k] = dft_mem_array[k][i];
        }
    }

    free_dft_memory(index, fft_dim, dft_mem_array);

    pca_free(index);
    if (pca_fit(index, samples, sample_size, fft_dim) != SUCCESS) {
        free(samples);
        free(input_data);
        return;
    }

    ts_type **coeff_mem_array = (ts_type **) calloc(n_segments, sizeof(ts_type *));
    if (coeff_mem_array == NULL) {
        free(samples);
        free(input_data);
        return;
    }
    for (int k = 0; k < n_segments; ++k) {
        coeff_mem_array[k] = (ts_type *) calloc(sample_size, sizeof(ts_type));
        if (coeff_mem_array[k] == NULL) {
            fprintf(stderr, "error: failed to allocate PISA projection buffer.\n");
            for (int j = 0; j < k; ++j) {
                free(coeff_mem_array[j]);
            }
            free(coeff_mem_array);
            free(samples);
            free(input_data);
            return;
        }
    }

    ts_type *projection = calloc(n_segments, sizeof(ts_type));
    if (projection == NULL) {
        for (int j = 0; j < n_segments; ++j) {
            free(coeff_mem_array[j]);
        }
        free(coeff_mem_array);
        free(samples);
        free(input_data);
        return;
    }

    for (unsigned int i = 0; i < sample_size; ++i) {
        if (pca_from_ts(index, samples + (i * fft_dim), projection) != SUCCESS) {
            continue;
        }
        for (int k = 0; k < n_segments; ++k) {
            coeff_mem_array[k][i] = projection[k];
        }
    }

    free(projection);
    free(samples);

    for (int i = 0; i < maxquerythread; i++) {
        input_data[i].index = index;
        input_data[i].dft_mem_array = coeff_mem_array;
        input_data[i].start_number = i * (n_segments / maxquerythread);
        input_data[i].stop_number = (i + 1) * (n_segments / maxquerythread);
    }

    input_data[maxquerythread - 1].start_number = (maxquerythread - 1) * (n_segments / maxquerythread);
    input_data[maxquerythread - 1].stop_number = n_segments;

    for (int i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, order_divide_worker, (void *) &(input_data[i]));
    }
    for (int i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }

    for (int k = 0; k < n_segments; ++k) {
        free(coeff_mem_array[k]);
    }
    free(coeff_mem_array);
    free(input_data);

    COUNT_BINNING_TIME_END

    sfa_print_bins(index);
    fprintf(stderr, ">>> Finished PISA binning\n");
}
