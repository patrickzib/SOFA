#ifndef MESSI_SPARTAN_H
#define MESSI_SPARTAN_H

#include "config.h"
#include "../../../globals.h"
#include "ads/isax_index.h"

enum response spartan_bins_init(isax_index *index);
void spartan_free_bins(isax_index *index);
enum response spartan_set_bins(isax_index *index, const char *ifilename, long int ts_num,
                               int maxquerythread, int filetype_int, int apply_znorm);
enum response collect_binning_samples(isax_index *index, const char *ifilename,
                                      long int ts_num, int filetype_int, int apply_znorm,
                                      ts_type *samples, unsigned int sample_size);
void spartan_print_bins(isax_index *index);

enum response spartan_from_ts(isax_index *index, const ts_type *ts, sax_type *sax_out,
                              ts_type *coeff_scratch);
void spartan_from_pca(isax_index *index, const ts_type *coeffs, sax_type *sax_out);
enum response pca_from_ts(const isax_index *index, const ts_type *ts, ts_type *out);

ts_type minidist_pca_to_spartan(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities, float bsf);

#endif /* MESSI_SPARTAN_H */
