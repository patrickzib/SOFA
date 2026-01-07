#ifndef sfalib_dft_h
#define sfalib_dft_h

#include "../../../config.h"
#include "../../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <fftw3.h>
#include <sys/types.h>

typedef struct fftw_workspace {
    ts_type *ts;
    fftwf_complex *ts_out;
    fftwf_plan plan_forward;
    ts_type *transform;
} fftw_workspace;

static inline void fftw_workspace_init(fftw_workspace *ws, unsigned long ts_length) {
    ws->ts = fftwf_malloc(sizeof(ts_type) * ts_length);
    ws->ts_out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * (ts_length / 2 + 1));
    ws->plan_forward = fftwf_plan_dft_r2c_1d(ts_length, ws->ts, ws->ts_out, FFTW_ESTIMATE);
    ws->transform = fftwf_malloc(sizeof(ts_type) * ts_length);
}

static inline void fftw_workspace_destroy(fftw_workspace *ws) {
    if (ws->plan_forward) {
        fftwf_destroy_plan(ws->plan_forward);
    }
    if (ws->ts) {
        fftwf_free(ws->ts);
    }
    if (ws->ts_out) {
        fftwf_free(ws->ts_out);
    }
    if (ws->transform) {
        fftwf_free(ws->transform);
    }
    ws->ts = NULL;
    ws->ts_out = NULL;
    ws->plan_forward = NULL;
    ws->transform = NULL;
}

void fft_from_ts(isax_index *index, int coeff_number, int best_only, fftw_workspace *fftw);
void sfa_from_fft(isax_index *index, ts_type * cur_transform, unsigned char * cur_sfa_word);

enum response sfa_from_ts(isax_index *index, sax_type *sax_out, fftw_workspace *fftw);

#endif
