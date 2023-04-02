//
// Created by wcy on 3/31/23.
//

#ifndef ALGORITHM_FFT_H
#define ALGORITHM_FFT_H

#include<immintrin.h>

typedef struct {
    __m256 re;
    __m256 im;
} complex_f8;

typedef struct {
    float re;
    float im;
} complex_f1;

void *avx_malloc(size_t size);

void avx_free(void *ptr);


#ifndef FFT_TYPE
typedef enum {
    Forward = -1,
    Inverse = 1
} fft_type;
#endif


typedef struct avx_fft {
    int N;
    int P[15];
    int P_n;

    complex_f8 *W8;
    complex_f8 *temp[8];

    complex_f1 **W;

    complex_f8 **time_domain;
    complex_f8 **freq_domain;
} AVX_FFT;

AVX_FFT *avx_fft_new(complex_f8 **time_domain, complex_f8 **freq_domain, int N);

void avx_fft_free(AVX_FFT *fft);

void avx_fft(AVX_FFT *fft, fft_type type);


#endif //ALGORITHM_FFT_H
