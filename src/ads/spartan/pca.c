#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "globals.h"
#include "ads/spartan/pca.h"

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

typedef struct {
    int piece;
    int component;
    int vector_index;
    double variance;
} pca_piece_candidate;

static int pca_compare_piece_candidate(const void *a, const void *b) {
    const pca_piece_candidate *left = (const pca_piece_candidate *) a;
    const pca_piece_candidate *right = (const pca_piece_candidate *) b;
    if (left->variance > right->variance) return -1;
    if (left->variance < right->variance) return 1;
    if (left->piece != right->piece) return left->piece - right->piece;
    return left->component - right->component;
}

/* Piecewise PCA remains an orthonormal transform: each local component is
 * zero outside its own disjoint time window.  Thus the existing squared
 * coefficient lower bounds stay valid without a separate representation. */
static enum response pca_fit_piecewise(isax_index *index, const ts_type *samples,
                                       unsigned int sample_size, int dim, int pieces) {
    const int components = index->settings->n_segments < dim
                               ? index->settings->n_segments : dim;
    if (pieces < 2 || pieces > dim || pieces > components) {
        fprintf(stderr, "error: SPARTAN PCA pieces must be between 1 and the retained component count (%d).\n",
                components);
        return FAILURE;
    }

    ts_type *mean = calloc((size_t) dim, sizeof(*mean));
    ts_type *candidate_vectors = calloc((size_t) dim * dim, sizeof(*candidate_vectors));
    pca_piece_candidate *candidates = calloc((size_t) dim, sizeof(*candidates));
    unsigned char *selected = calloc((size_t) dim, sizeof(*selected));
    ts_type *components_matrix = calloc((size_t) components * dim, sizeof(*components_matrix));
    double *explained_variance = calloc((size_t) components, sizeof(*explained_variance));
    if (mean == NULL || candidate_vectors == NULL || candidates == NULL || selected == NULL ||
        components_matrix == NULL || explained_variance == NULL) {
        free(mean); free(candidate_vectors); free(candidates); free(selected);
        free(components_matrix); free(explained_variance);
        fprintf(stderr, "error: failed to allocate piecewise PCA buffers.\n");
        return FAILURE;
    }

    const int base_width = dim / pieces;
    const int extra = dim % pieces;
    int candidate_count = 0;
    int start = 0;
    for (int piece = 0; piece < pieces; ++piece) {
        const int width = base_width + (piece < extra ? 1 : 0);
        double *cov = calloc((size_t) width * width, sizeof(*cov));
        double *eigvecs = calloc((size_t) width * width, sizeof(*eigvecs));
        double *eigvals = calloc((size_t) width, sizeof(*eigvals));
        if (cov == NULL || eigvecs == NULL || eigvals == NULL) {
            free(cov); free(eigvecs); free(eigvals);
            free(mean); free(candidate_vectors); free(candidates); free(selected);
            free(components_matrix); free(explained_variance);
            fprintf(stderr, "error: failed to allocate local PCA buffers.\n");
            return FAILURE;
        }
        for (unsigned int sample = 0; sample < sample_size; ++sample) {
            const ts_type *row = samples + (size_t) sample * dim + start;
            for (int d = 0; d < width; ++d) mean[start + d] += row[d];
        }
        for (int d = 0; d < width; ++d) mean[start + d] /= (ts_type) sample_size;
        for (unsigned int sample = 0; sample < sample_size; ++sample) {
            const ts_type *row = samples + (size_t) sample * dim + start;
            for (int a = 0; a < width; ++a) {
                const double va = (double) row[a] - mean[start + a];
                for (int b = a; b < width; ++b) {
                    const double vb = (double) row[b] - mean[start + b];
                    cov[a * width + b] += va * vb;
                }
            }
        }
        const double denom = sample_size > 1 ? (double) sample_size - 1.0 : 1.0;
        for (int a = 0; a < width; ++a) for (int b = a; b < width; ++b) {
            const double value = cov[a * width + b] / denom;
            cov[a * width + b] = value;
            cov[b * width + a] = value;
        }
        memcpy(eigvecs, cov, sizeof(*eigvecs) * (size_t) width * width);
        const int info = pca_lapack_eigen(eigvecs, width, eigvals);
        free(cov);
        if (info != 0) {
            free(eigvecs); free(eigvals);
            free(mean); free(candidate_vectors); free(candidates); free(selected);
            free(components_matrix); free(explained_variance);
            fprintf(stderr, "error: LAPACK local PCA failed (%d).\n", info);
            return FAILURE;
        }
        for (int local = 0; local < width; ++local) {
            candidates[candidate_count] = (pca_piece_candidate) {
                piece, local, candidate_count, eigvals[local]
            };
            ts_type *target = candidate_vectors + (size_t) candidate_count * dim + start;
            for (int d = 0; d < width; ++d) target[d] = (ts_type) eigvecs[d + local * width];
            ++candidate_count;
        }
        free(eigvecs); free(eigvals);
        start += width;
    }

    qsort(candidates, (size_t) candidate_count, sizeof(*candidates), pca_compare_piece_candidate);
    for (int piece = 0; piece < pieces; ++piece) {
        for (int candidate = 0; candidate < candidate_count; ++candidate) {
            if (candidates[candidate].piece == piece) { selected[candidate] = 1; break; }
        }
    }
    int selected_count = pieces;
    for (int candidate = 0; candidate < candidate_count && selected_count < components; ++candidate) {
        if (!selected[candidate]) { selected[candidate] = 1; ++selected_count; }
    }

    double variance_sum = 0.0;
    int output = 0;
    for (int candidate = 0; candidate < candidate_count && output < components; ++candidate) {
        if (!selected[candidate]) continue;
        memcpy(components_matrix + (size_t) output * dim,
               candidate_vectors + (size_t) candidates[candidate].vector_index * dim,
               sizeof(*components_matrix) * dim);
        explained_variance[output] = candidates[candidate].variance;
        variance_sum += candidates[candidate].variance;
        ++output;
    }
    for (int component = 0; component < components; ++component)
        explained_variance[component] = variance_sum > 0.0
            ? explained_variance[component] / variance_sum : 0.0;

    fprintf(stderr, ">>> SPARTAN: Piecewise PCA: %d pieces (window lengths %d--%d); "
                    "one component per piece, then global variance ranking\n",
            pieces, base_width, base_width + (extra > 0 ? 1 : 0));
    fprintf(stderr, ">>> SPARTAN: PCA eigenvalue ranking (top %d/%d):\n",
            components < 8 ? components : 8, components);
    int printed = 0;
    for (int candidate = 0; candidate < candidate_count && printed < 8; ++candidate) {
        if (!selected[candidate]) continue;
        fprintf(stderr, "piece %d, local %d (%.4f) ", candidates[candidate].piece,
                candidates[candidate].component, candidates[candidate].variance);
        ++printed;
    }
    if (components > 8) fprintf(stderr, "... (+%d more)", components - 8);
    fprintf(stderr, "\n");

    index->pca_mean = mean;
    index->pca_components = components_matrix;
    index->pca_explained_variance = explained_variance;
    index->pca_components_count = components;
    index->pca_dim = dim;
    free(candidate_vectors); free(candidates); free(selected);
    return SUCCESS;
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

    const int pieces = index->settings->spartan_pca_pieces;
    if (pieces > 1) return pca_fit_piecewise(index, samples, sample_size, dim, pieces);

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

    if (mean == NULL || cov == NULL || eigvecs == NULL || eigvals == NULL || ranked ==
        NULL || components_matrix == NULL) {
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
        return FAILURE;
    }

    for (int i = 0; i < dim; ++i) {
        ranked[i * 2] = (double) i;
        ranked[i * 2 + 1] = eigvals[i];
    }
    qsort(ranked, dim, sizeof(double) * 2, pca_compare_variance);

    fprintf(stderr, ">>> SPARTAN: PCA eigenvalue ranking (top %d/%d):\n",
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
        }
    }

    index->pca_mean = mean;
    index->pca_components = components_matrix;
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
