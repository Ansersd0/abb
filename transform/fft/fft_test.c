//
// Created by wcy on 2023/3/14.
//

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "fft2.h"
#include "real_fft.h"

void fft_2_test();

void fft_2_test_r();

void fft_2_test_s_1();

void fft_2_test_r_s_1();

int main() {
    fft_2_test();
    fft_2_test_r();
//    fft_2_test_s_1();
//    fft_2_test_r_s_1();
    return 0;
}


void fft_2_test() {
    int N = 1 << 14;
    f_complex *in, *out;
    in = malloc(sizeof(f_complex) * N);
    out = malloc(sizeof(f_complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i] = csetf(1 + i % 49, 2 + i % 94);
    }

    FFT_2_HEADER *fft = fft_2_new(&in, &out, N);

    struct timespec start, end;
    int count = 1000;
    long long sum = 0;

    for (int i = 0; i < count; ++i) {
//        for (int i = 0; i < N; ++i) {
//            in[i] = csetf(1 + i % 49, 2 + i % 94);
////        time_domain[i] = 1 + i;
//        }

        if (i > 0) {
            for (int i = 0; i < N; ++i) {
                in[i] = cdivf(in[i], csetf(N, 0));
//        time_domain[i] = 1 + i;
            }
        }

        clock_gettime(CLOCK_REALTIME, &start);

        fft_2(fft, Forward);
        fft_2(fft, Inverse);

        clock_gettime(CLOCK_REALTIME, &end);


        sum += (end.tv_sec - start.tv_sec) * 1000000000 + (end.tv_nsec - start.tv_nsec);
    }

//
    for (int i = 0; i < N % 51; ++i) {
        printf("%d: %f + %f * i\n", i, in[i].re / N, in[i].im / N);
    }

//    printf("Time taken %f us\n", ((double) (end - start)) / CLOCKS_PER_SEC * 1000000);
    printf("complex test Time taken %lf ns\n", (double) sum / count);

    fft_2_free(fft);
    free(in);
    free(out);


}


void fft_2_test_r() {
    int N = 1 <<15;
    float *in, *out;
    in = malloc(sizeof(float) * N);
    out = malloc(sizeof(float) * N);

//    for (int i = 0; i < N; ++i) {
//        in[i] = csetf(1 + i % 49, 2 + i % 94);
//    }

    PACK_HEADER *fft = pack_new(&in, &out, N);

    struct timespec start, end;

    int count = 1000;
    long long sum = 0;

    for (int i = 0; i < count; ++i) {
        for (int i = 0; i < N; ++i) {
            in[i] = 1 + i % 49;
//        time_domain[i] = 1 + i;
        }

        clock_gettime(CLOCK_REALTIME, &start);

        pack_fft(fft, Forward);
        pack_fft(fft, Inverse);

        clock_gettime(CLOCK_REALTIME, &end);

        sum += (end.tv_sec - start.tv_sec) * 1000000000 + (end.tv_nsec - start.tv_nsec);
    }

//
    for (int i = 0; i < N % 51; ++i) {
        printf("%d: %f\n", i, in[i] / N);
    }

//    printf("Time taken %f us\n", ((double) (end - start)) / CLOCKS_PER_SEC * 1000000);
    printf("real test Time taken %lf ns\n", (double) sum / count);

    pack_free(fft);
    free(in);
    free(out);

}

void fft_2_test_r_s_1() {
    int N = 1 << 25;
    float *in = malloc(sizeof(float) * N);
    float *out = malloc(sizeof(float) * N);

    PACK_HEADER *fft = pack_new(&in, &out, N);

    for (int i = 0; i < N; ++i) {
        in[i] = 1 + i % 49;
    }

    struct timespec start, end;

    clock_gettime(CLOCK_REALTIME, &start);
    pack_fft(fft, Forward);
    pack_fft(fft, Inverse);
    clock_gettime(CLOCK_REALTIME, &end);

    for (int i = 0; i < N % 51; ++i) {
        printf("%d: %f\n", i, in[i] / N);
    }

    printf("Time taken %lf ms\n",
           (double) (end.tv_sec - start.tv_sec) * 1000
           + (end.tv_nsec - start.tv_nsec) / 1000000);
}


void fft_2_test_s_1() {
    int N = 1 << 25;
    f_complex *in = malloc(sizeof(f_complex) * N);
    f_complex *out = malloc(sizeof(f_complex) * N);

    FFT_2_HEADER *fft = fft_2_new(&in, &out, N);

    for (int i = 0; i < N; ++i) {
        in[i] = csetf(1 + i % 49, 2 + i % 94);
    }

    struct timespec start, end;

    clock_gettime(CLOCK_REALTIME, &start);
    fft_2(fft, Forward);
    fft_2(fft, Inverse);
    clock_gettime(CLOCK_REALTIME, &end);

    for (int i = 0; i < N % 51; ++i) {
        printf("%d: %f + %f * i\n", i, in[i].re / N, in[i].im / N);
    }

    printf("Time taken %lf ms\n",
           (double) (end.tv_sec - start.tv_sec) * 1000
           + (end.tv_nsec - start.tv_nsec) / 1000000);
}