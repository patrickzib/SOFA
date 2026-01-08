#ifndef MESSI_SPARTAN_H
#define MESSI_SPARTAN_H

#include "config.h"
#include "../../../globals.h"
#include "ads/isax_index.h"

enum response spartan_bins_init(isax_index *index);
void spartan_free_bins(isax_index *index);
void spartan_set_bins(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread,
                      int filetype_int, int apply_znorm);
void spartan_print_bins(isax_index *index);

enum response spartan_from_ts(isax_index *index, const ts_type *ts, sax_type *sax_out);
void spartan_from_pca(isax_index *index, const ts_type *coeffs, sax_type *sax_out);
enum response pca_from_ts(const isax_index *index, const ts_type *ts, ts_type *out);

ts_type minidist_pca_to_spartan(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities, float bsf);
ts_type minidist_pca_to_spartan_raw(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities, float bsf);
ts_type minidist_pca_to_spartan_rawe_SIMD(isax_index *index, float *pca, sax_type *sax, sax_type *sax_cardinalities,
                                          float bsf);

#if !ADS_HAVE_AVX2
#define minidist_pca_to_spartan_rawe_SIMD minidist_pca_to_spartan_raw
#endif

#endif /* MESSI_SPARTAN_H */
