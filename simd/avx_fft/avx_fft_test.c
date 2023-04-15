//
// Created by wcy on 2023/3/31.
//

#include <time.h>
#include <stdio.h>
#include "avx_fft.h"


int test_1() {
    int N = 17;
    complex_f8 *time_domain = avx_malloc(sizeof(complex_f8) * N);
    complex_f8 *freq_domain = avx_malloc(sizeof(complex_f8) * N);

    for (int i = 0; i < N; ++i) {
        time_domain[i].re = _mm256_setr_ps(1 + i * 8, 2 + i * 8, 3 + i * 8, 4 + i * 8, 5 + i * 8, 6 + i * 8, 7 + i * 8,
                                           8 + i * 8);
        time_domain[i].im = _mm256_setr_ps(0, 0, 0, 0, 0, 0, 0, 0);
    }

    AVX_FFT *fft = avx_fft_new(&time_domain, &freq_domain, N);
    avx_fft(fft, Forward);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 8; ++j) {
            printf("%f %f\n", freq_domain[i].re[j], freq_domain[i].im[j]);
        }
    }

    avx_fft(fft, Inverse);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 8; ++j) {
            printf("%f %f\n", time_domain[i].re[j] / (N * 8), time_domain[i].im[j] / (N * 8));
        }
    }

    avx_fft_free(fft);
    avx_free(time_domain);
    avx_free(freq_domain);

    return 0;
}

int test_2() {
    int N = 1 << 4;
    complex_f8 *time_domain = avx_malloc(sizeof(complex_f8) * N);
    complex_f8 *freq_domain = avx_malloc(sizeof(complex_f8) * N);


    AVX_FFT *fft = avx_fft_new(&time_domain, &freq_domain, N);

    int count = 1;
    struct timespec start, end;
    long long sum = 0;
    clock_gettime(CLOCK_REALTIME, &start);

    for (int i = 0; i < count; ++i) {
        for (int i = 0; i < N; ++i) {
            time_domain[i].re = _mm256_setr_ps(1 + i * 8,
                                               2 + i * 8,
                                               3 + i * 8,
                                               4 + i * 8,
                                               5 + i * 8,
                                               6 + i * 8,
                                               7 + i * 8,
                                               8 + i * 8);
            time_domain[i].im = _mm256_setr_ps(0, 0, 0, 0, 0, 0, 0, 0);
        }

        clock_gettime(CLOCK_REALTIME, &start);

        avx_fft(fft, Forward);
        avx_fft(fft, Inverse);

        clock_gettime(CLOCK_REALTIME, &end);


        sum += (end.tv_sec - start.tv_sec) * 1000000000 + (end.tv_nsec - start.tv_nsec);
    }

    printf("avx fft time: %f ms\n", sum / (double) count / 1000000);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 8; ++j) {
            printf("%f + %fi\n",
                   time_domain[i].re[j] / (N * 8),
                   time_domain[i].im[j] / (N * 8));
        }
    }

    avx_fft_free(fft);
    avx_free(time_domain);
    avx_free(freq_domain);

    return 0;
}

int main() {
//    test_1();
    test_2();
    return 0;
}