#ifndef MESSI_PCA_H
#define MESSI_PCA_H

#include "config.h"
#include "../../../globals.h"
#include "ads/isax_index.h"

enum response pca_fit(isax_index *index, const ts_type *samples, unsigned int sample_size, int dim);
void pca_free(isax_index *index);

#endif /* MESSI_PCA_H */
