//
// Created by wcy on 2023/3/31.
//

#ifndef ALGORITHM_AVX_FFT_H
#define ALGORITHM_AVX_FFT_H

#include "../avx/avx.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef union c_to_f {
    __m256 m256;
    float f[8];
} c_to_f;

#define C8_MUL_C1(a, b) (complex_f8){_mm256_sub_ps(_mm256_mul_ps(a.re, _mm256_set1_ps(b.re)), \
                                                    _mm256_mul_ps(a.im, _mm256_set1_ps(b.im))), \
                                     _mm256_add_ps(_mm256_mul_ps(a.re, _mm256_set1_ps(b.im)), \
                                                   _mm256_mul_ps(a.im, _mm256_set1_ps(b.re)))}

#define C8_MUL_C1_CONJ(a, b) (complex_f8){_mm256_add_ps(_mm256_mul_ps(a.re, _mm256_set1_ps(b.re)), \
                                                        _mm256_mul_ps(a.im, _mm256_set1_ps(b.im))), \
                                          _mm256_sub_ps(_mm256_mul_ps(a.im, _mm256_set1_ps(b.re)), \
                                                        _mm256_mul_ps(a.re, _mm256_set1_ps(b.im)))}
#define C8_SET_C1(a_re, a_im) (complex_f8){_mm256_set1_ps(a_re), _mm256_set1_ps(a_im)}


#define C8_DIV_C1(a, b) (complex_f8){_mm256_div_ps(a.re, _mm256_set1_ps(b.re)), \
                                     _mm256_div_ps(a.im, _mm256_set1_ps(b.re))}

typedef struct complex_f1_t {
    float re;
    float im;
} complex_f1;

#define C1_SET_C1(a_re, a_im) (complex_f1){a_re, a_im}

#define C1_CONJ(a) (complex_f1){a.re, -a.im}

#define C1_ADD_C1(a, b) (complex_f1){a.re + b.re, a.im + b.im}
#define C1_SUB_C1(a, b) (complex_f1){a.re - b.re, a.im - b.im}
#define C1_MUL_C1(a, b) (complex_f1){a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}
#define C1_DIV_C1(a, b) (complex_f1){(a.re * b.re + a.im * b.im) / (b.re * b.re + b.im * b.im), \
                                     (a.im * b.re - a.re * b.im) / (b.re * b.re + b.im * b.im)}

#define C1_MUL_I(a) (complex_f1){-a.im, a.re}
#define C1_MUL_F1(a, b) (complex_f1){a.re * b, a.im * b}

#define C8_SET_8C1(a, b, c, d, e, f, g, h) \
    (complex_f8){_mm256_set_ps(a.re, b.re, c.re, d.re, e.re, f.re, g.re, h.re), \
                 _mm256_set_ps(a.im, b.im, c.im, d.im, e.im, f.im, g.im, h.im)}

typedef struct avx_fft_header {
    int N;
    int P[15];
    int P_n;

    complex_f8 *temp;

    complex_f8 *W8;
    complex_f1 **W;

    complex_f8 **time_domain;
    complex_f8 **freq_domain;
} AVX_FFT;

typedef enum fft_type {
    Forward = -1,
    Inverse = 1
} fft_type;

AVX_FFT *avx_fft_new(complex_f8 **time_domain, complex_f8 **freq_domain, fft_type N);

void avx_fft_free(AVX_FFT *fft);

void avx_fft(AVX_FFT *fft, fft_type type);


#endif //ALGORITHM_AVX_FFT_H
