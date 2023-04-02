//
// Created by wcy on 3/31/23.
//

#include <stdio.h>
#include "fft.h"

int main() {
    int N = 32;
    complex_f8 *time_domain = avx_malloc(sizeof(complex_f8) * N);
    complex_f8 *freq_domain = avx_malloc(sizeof(complex_f8) * N);

    for (int i = 0; i < N; ++i) {
        time_domain[i].re = _mm256_setr_ps(1 + i, 2 + i, 3 + i, 4 + i, 5 + i, 6 + i, 7 + i, 8 + i);
        time_domain[i].im = _mm256_setr_ps(0, 0, 0, 0, 0, 0, 0, 0);
    }

    AVX_FFT *fft = avx_fft_new(&time_domain, &freq_domain, N);
    avx_fft(fft, Forward);

    for (int i = 0; i < N; ++i) {
        printf("%f %f\n", freq_domain[i].re[0], freq_domain[i].im[0]);
    }

    avx_fft(fft, Inverse);

    for (int i = 0; i < N; ++i) {
        printf("%f %f\n", time_domain[i].re[0] / N, time_domain[i].im[0] / N);
    }

    avx_fft_free(fft);
    avx_free(time_domain);
    avx_free(freq_domain);

    return 0;
}