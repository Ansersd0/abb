//
// Created by wcy on 2023/3/15.
//

#ifndef ALGORITHM_FFT1_H
#define ALGORITHM_FFT1_H

#if defined(__GNUC__) || defined(__clang__)

#include <complex.h>
#include "common.h"

typedef struct {
    int N;
    int P[15];
    int P_n;

    float _Complex **W;
    float _Complex **time_domain;
    float _Complex **freq_domain;
} FFT_2_HEADER;

FFT_2_HEADER *fft_1_new(float _Complex **in, float _Complex **out, int N);

void fft_1_free(FFT_2_HEADER *fft);

void fft_1(FFT_2_HEADER *fft, fft_type type);

void my_fft(double _Complex *in, double _Complex *out, int N, fft_type type);

#endif //defined(__GNUC__) || defined(__clang__)
#endif //ALGORITHM_FFT1_H
