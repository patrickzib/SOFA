#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "globals.h"
#include "ads/spartan/pca.h"

static void pca_free_model(isax_index *index) {
    if (index->pca_mean != NULL) {
        free(index->pca_mean);
        index->pca_mean = NULL;
    }
    if (index->pca_components != NULL) {
        free(index->pca_components);
        index->pca_components = NULL;
    }
    index->pca_components_count = 0;
    index->pca_dim = 0;
}

void pca_free(isax_index *index) {
    if (index == NULL) {
        return;
    }
    pca_free_model(index);
}

static void pca_identity(double *matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

static void pca_jacobi_eigen(double *matrix, int n, double *eigvals, double *eigvecs) {
    const double eps = 1e-10;
    const int max_iter = 50 * n * n;

    pca_identity(eigvecs, n);

    for (int iter = 0; iter < max_iter; ++iter) {
        int p = 0;
        int q = 1;
        double max_off = 0.0;

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double value = fabs(matrix[i * n + j]);
                if (value > max_off) {
                    max_off = value;
                    p = i;
                    q = j;
                }
            }
        }

        if (max_off < eps) {
            break;
        }

        double app = matrix[p * n + p];
        double aqq = matrix[q * n + q];
        double apq = matrix[p * n + q];
        double phi = 0.5 * atan2(2.0 * apq, (aqq - app));
        double c = cos(phi);
        double s = sin(phi);

        for (int i = 0; i < n; ++i) {
            double aip = matrix[i * n + p];
            double aiq = matrix[i * n + q];
            matrix[i * n + p] = c * aip - s * aiq;
            matrix[i * n + q] = s * aip + c * aiq;
        }

        for (int i = 0; i < n; ++i) {
            double api = matrix[p * n + i];
            double aqi = matrix[q * n + i];
            matrix[p * n + i] = c * api - s * aqi;
            matrix[q * n + i] = s * api + c * aqi;
        }

        matrix[p * n + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        matrix[q * n + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        matrix[p * n + q] = 0.0;
        matrix[q * n + p] = 0.0;

        for (int i = 0; i < n; ++i) {
            double vip = eigvecs[i * n + p];
            double viq = eigvecs[i * n + q];
            eigvecs[i * n + p] = c * vip - s * viq;
            eigvecs[i * n + q] = s * vip + c * viq;
        }
    }

    for (int i = 0; i < n; ++i) {
        eigvals[i] = matrix[i * n + i];
    }
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

enum response pca_fit(isax_index *index, const ts_type *samples, unsigned int sample_size, int dim) {
    int components = index->settings->n_segments;
    if (components > dim) {
        components = dim;
    }

    ts_type *mean = calloc((size_t) dim, sizeof(ts_type));
    double *cov = calloc((size_t) dim * dim, sizeof(double));
    double *eigvecs = calloc((size_t) dim * dim, sizeof(double));
    double *eigvals = calloc((size_t) dim, sizeof(double));
    double *ranked = calloc((size_t) dim * 2, sizeof(double));
    ts_type *components_matrix = calloc((size_t) components * dim, sizeof(ts_type));

    if (mean == NULL || cov == NULL || eigvecs == NULL || eigvals == NULL || ranked == NULL || components_matrix == NULL) {
        free(mean);
        free(cov);
        free(eigvecs);
        free(eigvals);
        free(ranked);
        free(components_matrix);
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

    pca_jacobi_eigen(cov, dim, eigvals, eigvecs);

    for (int i = 0; i < dim; ++i) {
        ranked[i * 2] = (double) i;
        ranked[i * 2 + 1] = eigvals[i];
    }
    qsort(ranked, dim, sizeof(double) * 2, pca_compare_variance);

    for (int k = 0; k < components; ++k) {
        int idx = (int) ranked[k * 2];
        for (int i = 0; i < dim; ++i) {
            components_matrix[k * dim + i] = (ts_type) eigvecs[i * dim + idx];
        }
    }

    index->pca_mean = mean;
    index->pca_components = components_matrix;
    index->pca_components_count = components;
    index->pca_dim = dim;

    free(cov);
    free(eigvecs);
    free(eigvals);
    free(ranked);

    return SUCCESS;
}
