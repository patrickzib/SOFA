#ifndef MESSI_PISA_H
#define MESSI_PISA_H

#include "config.h"
#include "../../../globals.h"
#include "ads/isax_index.h"
#include "ads/sfa/dft.h"

enum response pisa_bins_init(isax_index *index);
void pisa_free_bins(isax_index *index);
void pisa_set_bins(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread,
                   int filetype_int, int apply_znorm);

enum response pisa_pca_from_ts(isax_index *index, const ts_type *ts, ts_type *out, fftw_workspace *fftw);
enum response pisa_from_ts(isax_index *index, const ts_type *ts, sax_type *sax_out, fftw_workspace *fftw);

#endif /* MESSI_PISA_H */
