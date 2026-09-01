#include "config.h"
#include "ads/sffa/sffa.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void sffa_workspace_init(sffa_workspace *ws, unsigned long length) {
    if (ws == NULL || length == 0) return;
    memset(ws, 0, sizeof(*ws));
    ws->input = fftwf_malloc(sizeof(*ws->input) * length);
    ws->forward = fftwf_malloc(sizeof(*ws->forward) * length);
    ws->inverse = fftwf_malloc(sizeof(*ws->inverse) * length);
    ws->output = fftwf_malloc(sizeof(*ws->output) * length);
    if (ws->input == NULL || ws->forward == NULL || ws->inverse == NULL || ws->output == NULL) return;
    ws->forward_plan = fftwf_plan_dft_1d((int) length, ws->input, ws->forward,
                                         FFTW_FORWARD, FFTW_ESTIMATE);
    ws->inverse_plan = fftwf_plan_dft_1d((int) length, ws->input, ws->inverse,
                                         FFTW_BACKWARD, FFTW_ESTIMATE);
}

void sffa_workspace_destroy(sffa_workspace *ws) {
    if (ws == NULL) return;
    if (ws->forward_plan) fftwf_destroy_plan(ws->forward_plan);
    if (ws->inverse_plan) fftwf_destroy_plan(ws->inverse_plan);
    if (ws->input) fftwf_free(ws->input);
    if (ws->forward) fftwf_free(ws->forward);
    if (ws->inverse) fftwf_free(ws->inverse);
    if (ws->output) fftwf_free(ws->output);
    memset(ws, 0, sizeof(*ws));
}

/* The normalized DFT F has F^4=I.  Its four orthogonal eigenspace
 * projectors yield F^p = sum_m c_m F^m, with
 * c_m = 1/4 sum_k exp(i*pi*k*(p-m)/2).  This avoids an approximate FrFT:
 * the transform used by the index is unitary to FFTW roundoff. */
enum response sffa_fractional_transform(const isax_index *index,
                                        const ts_type *ts,
                                        sffa_workspace *ws) {
    if (index == NULL || index->settings == NULL || ts == NULL || ws == NULL ||
        ws->input == NULL || ws->output == NULL || ws->forward_plan == NULL ||
        ws->inverse_plan == NULL) return FAILURE;
    const int n = index->settings->timeseries_size;
    if (n <= 0) return FAILURE;
    const float inv_sqrt_n = 1.0f / sqrtf((float) n);
    /* FFTW's backward DFT is unnormalized.  F^3 is the *normalized* inverse
     * DFT, so it needs the same 1/sqrt(N) scale as the forward branch. */
    const float inverse_scale = inv_sqrt_n;
    float cr[4] = {0}, ci[4] = {0};
    const double order = index->settings->sffa_order;
    for (int m = 0; m < 4; ++m) {
        for (int k = 0; k < 4; ++k) {
            const double phase = M_PI * (double) k * (order - (double) m) / 2.0;
            cr[m] += (float) (0.25 * cos(phase));
            ci[m] += (float) (0.25 * sin(phase));
        }
    }
    for (int i = 0; i < n; ++i) {
        ws->input[i][0] = ts[i];
        ws->input[i][1] = 0.0f;
    }
    fftwf_execute(ws->forward_plan);
    /* The inverse plan shares its input, so restore the real input first. */
    for (int i = 0; i < n; ++i) {
        ws->input[i][0] = ts[i];
        ws->input[i][1] = 0.0f;
    }
    fftwf_execute(ws->inverse_plan);
    for (int i = 0; i < n; ++i) {
        const int reverse = i == 0 ? 0 : n - i;
        float fr[4] = {ts[i], ws->forward[i][0] * inv_sqrt_n,
                       ts[reverse], ws->inverse[i][0] * inverse_scale};
        float fi[4] = {0.0f, ws->forward[i][1] * inv_sqrt_n,
                       0.0f, ws->inverse[i][1] * inverse_scale};
        float real = 0.0f, imag = 0.0f;
        for (int m = 0; m < 4; ++m) {
            real += cr[m] * fr[m] - ci[m] * fi[m];
            imag += cr[m] * fi[m] + ci[m] * fr[m];
        }
        ws->output[i][0] = real;
        ws->output[i][1] = imag;
    }
    return SUCCESS;
}
