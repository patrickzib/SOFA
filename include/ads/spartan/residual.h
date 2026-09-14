#ifndef MESSI_SPARTAN_RESIDUAL_H
#define MESSI_SPARTAN_RESIDUAL_H
#include <math.h>
#include "ads/isax_index.h"
/* Assumes orthonormal PCA; no floating-point error correction.
 * Suffix and residual overlap and must be combined using max. */
typedef struct { int enabled; } spartan_residual_model;
typedef struct { float radius; } spartan_residual_value;
static inline spartan_residual_model spartan_residual_init(const isax_index *index) {
    spartan_residual_model model = {0};
    if (index->pca_dim <= 0 || index->pca_components_count <= 0 ||
        !index->pca_mean || !index->pca_components || !index->pca_bias) return model;
    for (int d = 0; d < index->pca_dim; ++d)
        if (!isfinite(index->pca_mean[d])) return model;
    for (int i = 0; i < index->pca_components_count; ++i) {
        if (!isfinite(index->pca_bias[i])) return model;
        for (int d = 0; d < index->pca_dim; ++d)
            if (!isfinite(index->pca_components[(size_t)i*index->pca_dim+d])) return model;
    }
    model.enabled = 1;
    return model;
}
static inline spartan_residual_value spartan_residual_encode(
        const isax_index *index, const spartan_residual_model *model,
        const float *raw, const float *projection, int k) {
    spartan_residual_value value = {NAN};
    if (!model->enabled || k < 0 || k > index->pca_components_count) return value;
    double norm = 0.0, prefix = 0.0;
    for (int d = 0; d < index->pca_dim; ++d) {
        if (!isfinite(raw[d])) return value;
        double centered = (double)raw[d]-index->pca_mean[d];
        norm += centered*centered;
    }
    for (int i = 0; i < k; ++i) {
        if (!isfinite(projection[i])) return value;
        prefix += (double)projection[i]*projection[i];
    }
    value.radius = (float)sqrt(fmax(0.0, norm-prefix));
    return value;
}
static inline double spartan_residual_gap(double q, double low, double high) {
    if (!isfinite(q) || !isfinite(low) || !isfinite(high)) return 0.0;
    return fmax(0.0, fmax(low-q, q-high));
}
#endif
