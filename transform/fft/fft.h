//
// Created by wcy on 2023/2/20.
//

#ifndef ALGORITHM_FFT_H
#define ALGORITHM_FFT_H

#include <immintrin.h>

#include <complex.h>

#define PI 3.14159265358979323846
#define VEC_SIZE 8

typedef union {
    __m256 v;
    float f[VEC_SIZE];
} vec8f;
#define V8f_ADD(a, b) _mm256_add_ps(a, b)
#define V8f_SUB(a, b) _mm256_sub_ps(a, b)
#define V8f_MUL(a, b) _mm256_mul_ps(a, b)
#define V8f_DIV(a, b) _mm256_div_ps(a, b)
#define V8f_SET(a0, a1, a2, a3, a4, a5, a6, a7) _mm256_set_ps(a0, a1, a2, a3, a4, a5, a6, a7)

#define V8f_ADD(a, b) _mm256_add_ps(a, b)
#define V8f_SUB(a, b) _mm256_sub_ps(a, b)
#define V8f_MUL(a, b) _mm256_mul_ps(a, b)


typedef union {
    __m256 v;
    float _Complex c[VEC_SIZE / 2];
} vec4c;

#define V4c_ADD(a, b) _mm256_add_ps(a, b)
#define V4c_SUB(a, b) _mm256_sub_ps(a, b)
#define V4c_MUL(out, a, b)      vec4c a1, a2, c1, c2;\
a1.v = _mm256_permutevar8x32_ps(a.v, _mm256_setr_epi32(0, 0, 2, 2, 4, 4, 6, 6));\
a2.v = _mm256_permutevar8x32_ps(a.v, _mm256_setr_epi32(1, 1, 3, 3, 5, 5, 7, 7));\
c1.v = _mm256_mul_ps(a1.v, b.v);\
c2.v = _mm256_mul_ps(a2.v, b.v);\
out.v = _mm256_addsub_ps(c1.v, _mm256_permutevar8x32_ps(c2.v, _mm256_setr_epi32(1, 0, 3, 2, 5, 4, 7, 6)));
#define V4c_DIV(a, b) _mm256_div_ps(a, b)

enum FFT_TYPE {
    Forward = -1,
    Inverse = 1
};

typedef struct {
    int N_float;
    int m_cev;
    vec8f *W;

    int P_n;
    int P[15]
} FFT;

int decompose(int N, int *p);


void fft_test();

typedef struct {
    int N;
    int P[15];
    int P_n;

    float _Complex **W;
} FFT_1_HEADER;

FFT_1_HEADER *fft_1_new(int N);

void fft_1_free(FFT_1_HEADER *fft);

void fft_1(FFT_1_HEADER *fft, float _Complex *__restrict__ in, float _Complex *__restrict__ out, int type);

 void fft_1_paas_4(float _Complex *__restrict__ in,
                  float _Complex *__restrict__ out,
                  float _Complex *__restrict__ W,
                  int P1, int P2, int P3, int type);

 void fft_1_paas_2(float _Complex *__restrict__ in,
                  float _Complex *__restrict__ out,
                  float _Complex *__restrict__ W,
                  int P1, int P2, int P3, int type);

 void fft_1_paas_n(float _Complex *__restrict__ in,
                  float _Complex *__restrict__ out,
                  float _Complex *__restrict__ W,
                  int P1, int P2, int P3, int type);

void fft_1_test();

#endif //ALGORITHM_FFT_H
