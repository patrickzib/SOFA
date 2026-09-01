#ifndef MESSI_SFFA_H
#define MESSI_SFFA_H

#include "ads/isax_index.h"
#include <fftw3.h>

/* A scalar coordinate in the complex fractional-Fourier representation. */
typedef struct sffa_coordinate {
    int coefficient;
    unsigned char imaginary;
} sffa_coordinate;

typedef struct sffa_workspace {
    fftwf_complex *input;
    fftwf_complex *forward;
    fftwf_complex *inverse;
    fftwf_complex *output;
    fftwf_plan forward_plan;
    fftwf_plan inverse_plan;
} sffa_workspace;

void sffa_workspace_init(sffa_workspace *ws, unsigned long length);
void sffa_workspace_destroy(sffa_workspace *ws);

/* Compute the unitary fractional DFT.  It is the spectral fractional power
 * of the normalized DFT: order 0 is identity, order 1 is the normalized DFT.
 */
enum response sffa_fractional_transform(const isax_index *index,
                                        const ts_type *ts,
                                        sffa_workspace *ws);
enum response sffa_project(const isax_index *index, const ts_type *ts,
                           ts_type *out, sffa_workspace *ws);
void sffa_from_values(const isax_index *index, const ts_type *values,
                      sax_type *word);
enum response sffa_from_ts(const isax_index *index, const ts_type *ts,
                           sax_type *word, sffa_workspace *ws,
                           ts_type *values);

enum response sffa_bins_init(isax_index *index);
void sffa_free_bins(isax_index *index);
enum response sffa_set_bins(isax_index *index, const char *filename,
                            long ts_num, int maxquerythread,
                            int filetype_int, int apply_znorm);

#endif /* MESSI_SFFA_H */
