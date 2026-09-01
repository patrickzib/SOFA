#include "ads/sffa/sffa.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int close_enough(float actual, float expected) {
    return fabsf(actual - expected) < 2e-4f;
}

int main(void) {
    enum { N = 8 };
    isax_index_settings settings;
    isax_index index;
    ts_type input[N] = {1.0f, -2.0f, 0.5f, 3.0f, -1.0f, 0.25f, 2.0f, -0.75f};
    memset(&settings, 0, sizeof(settings));
    memset(&index, 0, sizeof(index));
    settings.timeseries_size = N;
    index.settings = &settings;
    sffa_workspace ws = {0};
    sffa_workspace_init(&ws, N);
    if (ws.forward_plan == NULL || ws.inverse_plan == NULL) return 1;

    settings.sffa_order = 0.0;
    if (sffa_fractional_transform(&index, input, &ws) != SUCCESS) return 1;
    for (int i = 0; i < N; ++i) {
        if (!close_enough(ws.output[i][0], input[i]) || !close_enough(ws.output[i][1], 0.0f)) {
            fprintf(stderr, "SFFA p=0 identity mismatch at %d\n", i);
            return 1;
        }
    }

    settings.sffa_order = 1.0;
    if (sffa_fractional_transform(&index, input, &ws) != SUCCESS) return 1;
    for (int k = 0; k < N; ++k) {
        float real = 0.0f, imag = 0.0f;
        for (int n = 0; n < N; ++n) {
            const float phase = -2.0f * (float) M_PI * k * n / N;
            real += input[n] * cosf(phase) / sqrtf(N);
            imag += input[n] * sinf(phase) / sqrtf(N);
        }
        if (!close_enough(ws.output[k][0], real) || !close_enough(ws.output[k][1], imag)) {
            fprintf(stderr, "SFFA p=1 DFT mismatch at %d\n", k);
            return 1;
        }
    }

    settings.sffa_order = 3.0;
    if (sffa_fractional_transform(&index, input, &ws) != SUCCESS) return 1;
    for (int k = 0; k < N; ++k) {
        float real = 0.0f, imag = 0.0f;
        for (int n = 0; n < N; ++n) {
            const float phase = 2.0f * (float) M_PI * k * n / N;
            real += input[n] * cosf(phase) / sqrtf(N);
            imag += input[n] * sinf(phase) / sqrtf(N);
        }
        if (!close_enough(ws.output[k][0], real) || !close_enough(ws.output[k][1], imag)) {
            fprintf(stderr, "SFFA p=3 inverse-DFT mismatch at %d\n", k);
            return 1;
        }
    }

    settings.sffa_order = 0.5;
    if (sffa_fractional_transform(&index, input, &ws) != SUCCESS) return 1;
    float input_energy = 0.0f, output_energy = 0.0f;
    for (int i = 0; i < N; ++i) {
        input_energy += input[i] * input[i];
        output_energy += ws.output[i][0] * ws.output[i][0] + ws.output[i][1] * ws.output[i][1];
    }
    sffa_workspace_destroy(&ws);
    if (!close_enough(input_energy, output_energy)) {
        fprintf(stderr, "SFFA fractional transform is not unitary: %g != %g\n", input_energy, output_energy);
        return 1;
    }
    return 0;
}
