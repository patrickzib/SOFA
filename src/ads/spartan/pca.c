#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "globals.h"
#include "ads/build_progress.h"
#include "ads/spartan/pca.h"

#if HAVE_CBLAS
#include CBLAS_HEADER
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

extern void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                   double *w, double *work, int *lwork, int *info);

static void pca_free_model(isax_index *index) {
    if (index->pca_mean != NULL) {
        free(index->pca_mean);
        index->pca_mean = NULL;
    }
    if (index->pca_components != NULL) {
        free(index->pca_components);
        index->pca_components = NULL;
    }
    free(index->pca_bias);
    index->pca_bias = NULL;
    index->pca_components_count = 0;
    index->pca_dim = 0;
    free(index->pca_explained_variance);
    index->pca_explained_variance = NULL;
}

void pca_free(isax_index *index) {
    if (index == NULL) {
        return;
    }
    pca_free_model(index);
}

static int pca_lapack_eigen(double *matrix, int n, double *eigvals) {
    char jobz = 'V';
    char uplo = 'U';
    int lda = n;
    int info = 0;
    int lwork = -1;
    double wkopt = 0.0;

    dsyev_(&jobz, &uplo, &n, matrix, &lda, eigvals, &wkopt, &lwork, &info);
    if (info != 0) {
        return info;
    }

    lwork = (int) wkopt;
    if (lwork < 1) {
        return -1;
    }

    double *work = malloc(sizeof(double) * (size_t) lwork);
    if (work == NULL) {
        return -1;
    }

    dsyev_(&jobz, &uplo, &n, matrix, &lda, eigvals, work, &lwork, &info);
    free(work);
    return info;
}

static int pca_compare_variance(const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    if (da[1] > db[1]) {
        return -1;
    }
    if (da[1] < db[1]) {
        return 1;
    }
    return 0;
}

enum response pca_fit(
    isax_index *index,
    const ts_type *samples,
    unsigned int sample_size,
    int dim) {
    if (index == NULL ||
        index->settings == NULL ||
        samples == NULL ||
        sample_size == 0 ||
        dim <= 0) {
        fprintf(stderr, "error: invalid PCA input.\n");
        return FAILURE;
    }

    int components = index->settings->n_segments;
    if (components <= 0) {
        fprintf(stderr, "error: PCA requires at least one component.\n");
        return FAILURE;
    }

    if (components > dim) {
        components = dim;
    }

    ts_type *mean = calloc((size_t) dim, sizeof(ts_type));
    double *cov = calloc((size_t) dim * dim, sizeof(double));
    double *eigvecs = calloc((size_t) dim * dim, sizeof(double));
    double *eigvals = calloc((size_t) dim, sizeof(double));
    double *ranked = calloc((size_t) dim * 2, sizeof(double));
    ts_type *components_matrix = calloc((size_t) components * dim, sizeof(ts_type));
    double *bias = calloc((size_t) components, sizeof(double));

    if (mean == NULL || cov == NULL || eigvecs == NULL || eigvals == NULL || ranked ==
        NULL || components_matrix == NULL || bias == NULL) {
        free(mean);
        free(cov);
        free(eigvecs);
        free(eigvals);
        free(ranked);
        free(components_matrix);
        free(bias);
        fprintf(stderr, "error: failed to allocate PCA buffers.\n");
        return FAILURE;
    }

    for (unsigned int i = 0; i < sample_size; ++i) {
        const ts_type *row = samples + (i * dim);
        for (int j = 0; j < dim; ++j) {
            mean[j] += row[j];
        }
    }

    for (int j = 0; j < dim; ++j) {
        mean[j] /= (ts_type) sample_size;
    }

    for (unsigned int i = 0; i < sample_size; ++i) {
        const ts_type *row = samples + (i * dim);
        for (int a = 0; a < dim; ++a) {
            double va = (double) row[a] - mean[a];
            for (int b = a; b < dim; ++b) {
                double vb = (double) row[b] - mean[b];
                cov[a * dim + b] += va * vb;
            }
        }
    }

    double denom = (sample_size > 1) ? (double) (sample_size - 1) : 1.0;
    for (int a = 0; a < dim; ++a) {
        for (int b = a; b < dim; ++b) {
            double value = cov[a * dim + b] / denom;
            cov[a * dim + b] = value;
            cov[b * dim + a] = value;
        }
    }

    memcpy(eigvecs, cov, sizeof(double) * (size_t) dim * (size_t) dim);
    int info = pca_lapack_eigen(eigvecs, dim, eigvals);
    if (info != 0) {
        fprintf(stderr, "error: LAPACK dsyev failed (%d).\n", info);
        free(mean);
        free(cov);
        free(eigvecs);
        free(eigvals);
        free(ranked);
        free(components_matrix);
        free(bias);
        return FAILURE;
    }

    for (int i = 0; i < dim; ++i) {
        ranked[i * 2] = (double) i;
        ranked[i * 2 + 1] = eigvals[i];
    }
    qsort(ranked, dim, sizeof(double) * 2, pca_compare_variance);

    fprintf(stderr, ">>> PCA eigenvalue ranking (top %d/%d):\n",
            components < 8 ? components : 8, components);
    const int print_components = components < 8 ? components : 8;
    for (int i = 0; i < print_components; ++i) {
        fprintf(stderr, "%d, (%.4f) ", (int) ranked[i * 2], ranked[i * 2 + 1]);
    }
    if (print_components < components) {
        fprintf(stderr, "... (+%d more)", components - print_components);
    }
    fprintf(stderr, "\n");

    for (int k = 0; k < components; ++k) {
        int idx = (int) ranked[k * 2];
        for (int i = 0; i < dim; ++i) {
            components_matrix[k * dim + i] = (ts_type) eigvecs[i + idx * dim];
            bias[k] -= (double) mean[i] * components_matrix[k * dim + i];
        }
    }

    index->pca_mean = mean;
    index->pca_components = components_matrix;
    index->pca_bias = bias;
    index->pca_explained_variance = calloc((size_t) index->settings->n_segments,
                                           sizeof(*index->pca_explained_variance));
    if (index->pca_explained_variance == NULL) {
        pca_free_model(index);
        free(cov);
        free(eigvecs);
        free(eigvals);
        free(ranked);
        return FAILURE;
    }
    double variance_sum = 0.0;
    for (int k = 0; k < components; ++k) variance_sum += ranked[k * 2 + 1];
    for (int k = 0; k < components; ++k) {
        index->pca_explained_variance[k] = variance_sum > 0.0
            ? ranked[k * 2 + 1] / variance_sum : 0.0;
    }
    index->pca_components_count = components;
    index->pca_dim = dim;

    free(cov);
    free(eigvecs);
    free(eigvals);
    free(ranked);

    return SUCCESS;
}

void pca_report_projection_backend(int workers, unsigned int rows) {
    static int reported_projection_mode = 0;
    if (reported_projection_mode) return;

    messi_build_progress_begin_diagnostic();
#if HAVE_CBLAS
#if HAVE_OPENBLAS_CONFIG
    const char *openblas_config = openblas_get_config();
    const int parallel_sgemm = openblas_config != NULL &&
        strstr(openblas_config, "SINGLE_THREADED") != NULL && workers > 1 && rows > 1;
    if (parallel_sgemm) {
        fprintf(stderr, ">>> PCA projection: CBLAS SGEMM; OpenBLAS is single-threaded, using %d OpenMP workers\n",
                workers < (int) rows ? workers : (int) rows);
    } else if (openblas_config != NULL) {
        fprintf(stderr, ">>> PCA projection: CBLAS SGEMM; OpenBLAS threaded (%s)\n",
                openblas_config);
    } else {
        fprintf(stderr, ">>> PCA projection: CBLAS SGEMM\n");
    }
#else
    fprintf(stderr, ">>> PCA projection: CBLAS SGEMM\n");
#endif
#else
    fprintf(stderr, ">>> PCA projection: scalar fallback (CBLAS not compiled)\n");
#endif
    messi_build_progress_end_diagnostic();
    reported_projection_mode = 1;
}

enum response pca_project_batch(const isax_index *index, const ts_type *input,
                                unsigned int rows, ts_type *output, int workers) {
    if (index == NULL || input == NULL || output == NULL || rows == 0 ||
        index->settings == NULL || index->pca_components == NULL || index->pca_bias == NULL ||
        index->pca_dim <= 0 || index->pca_components_count <= 0) {
        return FAILURE;
    }

    const int input_dim = index->pca_dim;
    const int output_dim = index->settings->n_segments;
    const int components = index->pca_components_count;
    if (components > output_dim) return FAILURE;

    pca_report_projection_backend(workers, rows);

#if HAVE_CBLAS
    memset(output, 0, sizeof(*output) * (size_t) rows * output_dim);
    int parallel_sgemm = 0;
#if HAVE_OPENBLAS_CONFIG
    const char *openblas_config = openblas_get_config();
#ifdef _OPENMP
    parallel_sgemm = openblas_config != NULL &&
        strstr(openblas_config, "SINGLE_THREADED") != NULL && workers > 1 && rows > 1;
#endif
#endif
    if (parallel_sgemm) {
        const int worker_count = workers < (int) rows ? workers : (int) rows;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
        for (int worker = 0; worker < worker_count; ++worker) {
            const unsigned int first = (unsigned int) (((unsigned long long) rows * worker) / worker_count);
            const unsigned int last = (unsigned int) (((unsigned long long) rows * (worker + 1)) / worker_count);
            const unsigned int count = last - first;
            if (count == 0) continue;
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        (int) count, components, input_dim, 1.0f,
                        input + (size_t) first * input_dim, input_dim,
                        index->pca_components, input_dim, 0.0f,
                        output + (size_t) first * output_dim, output_dim);
        }
    } else {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    (int) rows, components, input_dim, 1.0f,
                    input, input_dim, index->pca_components, input_dim,
                    0.0f, output, output_dim);
    }
    for (unsigned int row = 0; row < rows; ++row) {
        ts_type *projected = output + (size_t) row * output_dim;
        for (int k = 0; k < components; ++k) projected[k] += (ts_type) index->pca_bias[k];
    }
    return SUCCESS;
#else
    for (unsigned int row = 0; row < rows; ++row) {
        if (pca_from_ts(index, input + (size_t) row * input_dim,
                        output + (size_t) row * output_dim) != SUCCESS)
            return FAILURE;
    }
    return SUCCESS;
#endif
}
