//
// Created by wcy on 2023/2/20.
//

#ifndef ALGORITHM_FFT2_H
#define ALGORITHM_FFT2_H


#include "common.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


typedef struct {
    float re;
    float im;
} f_complex;

#define _cplx_I ((f_complex) {0, 1})

#define caddf(a, b) ((f_complex) {a.re + b.re, a.im + b.im})
#define csubf(a, b) ((f_complex) {a.re - b.re, a.im - b.im})
#define cmulf(a, b) ((f_complex) {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re})
#define cmulf_I(a) ((f_complex) {-a.im, a.re})
#define smulf(a) ((f_complex) {-a, 0})
#define smulf_I(a) ((f_complex) {0, a})
#define cdivf(a, b) ((f_complex) {(a.re * b.re + a.im * b.im) / (b.re * b.re + b.im * b.im), \
                                  (a.im * b.re - a.re * b.im) / (b.re * b.re + b.im * b.im)})
#define cconjf(a) ((f_complex) {a.re, -a.im})
#define cexpf(a) ((f_complex) {expf(a.re) * cosf(a.im),  expf(a.re)*sinf(a.im)})
#define csetf(a, b) ((f_complex) {a, b})

typedef struct fft_2_header{
    int N;
    int P[15];
    int P_n;

    f_complex **W;

    f_complex **time_domain;
    f_complex **freq_domain;
} FFT_2_HEADER;

FFT_2_HEADER *fft_2_new(f_complex **time_domain, f_complex **freq_domain, fft_type N);

void fft_2_free(FFT_2_HEADER *fft);

void fft_2(FFT_2_HEADER *fft, fft_type type);

#endif //ALGORITHM_FFT2_H
