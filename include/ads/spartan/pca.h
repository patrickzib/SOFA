#ifndef MESSI_PCA_H
#define MESSI_PCA_H

#include "config.h"
#include "../../../globals.h"
#include "ads/isax_index.h"

enum response pca_fit(isax_index *index, const ts_type *samples, unsigned int sample_size, int dim);
void pca_free(isax_index *index);
enum response pca_from_ts(const isax_index *index, const ts_type *ts, ts_type *out);
/* Project contiguous input rows into contiguous output rows.  Uses CBLAS
 * SGEMM when configured, otherwise the scalar projection path. */
enum response pca_project_batch(const isax_index *index, const ts_type *input,
                                unsigned int rows, ts_type *output);

#endif /* MESSI_PCA_H */
