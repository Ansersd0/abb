//
// Created by wcy on 2023/3/15.
//

#ifndef ALGORITHM_AVX_H
#define ALGORITHM_AVX_H

#include <immintrin.h>


#ifndef __AVX__
#error "AVX is not supported by this compiler"
#endif

typedef struct complex_c8_t {
    __m256 re;
    __m256 im;
} complex_f8;


#define C8_SET_16F1(a_re, a_im, b_re, b_im, c_re, c_im, d_re, d_im, e_re, e_im, f_re, f_im, g_re, g_im, h_re, h_im) \
    (complex_f8){_mm256_set_ps(a_re, b_re, c_re, d_re, e_re, f_re, g_re, h_re), \
                 _mm256_set_ps(a_im, b_im, c_im, d_im, e_im, f_im, g_im, h_im)}

#define C8_MUL_I(a) (complex_f8){ -a.im, a.re}

#define C8_CONJ(a) (complex_f8){a.re, -a.im}

#define C8_ADD_C8(a, b) (complex_f8){_mm256_add_ps(a.re, b.re), _mm256_add_ps(a.im, b.im)}

#define C8_SUB_C8(a, b) (complex_f8){_mm256_sub_ps(a.re, b.re), _mm256_sub_ps(a.im, b.im)}

#define C8_MUL_C8(a, b) (complex_f8){_mm256_sub_ps(_mm256_mul_ps(a.re, b.re), _mm256_mul_ps(a.im, b.im)), \
                                     _mm256_add_ps(_mm256_mul_ps(a.re, b.im), _mm256_mul_ps(a.im, b.re))}

#define C8_MUL_C8_CONJ(a, b) (complex_f8){_mm256_add_ps(_mm256_mul_ps(a.re, b.re), _mm256_mul_ps(a.im, b.im)), \
                                          _mm256_sub_ps(_mm256_mul_ps(a.im, b.re), _mm256_mul_ps(a.re, b.im))}
#define C8_MUL_F8(a, b) (complex_f8){_mm256_mul_ps(a.re, b), \
                                     _mm256_mul_ps(a.im, b)}


#define C8_MUL_F1(a, b) (complex_f8){_mm256_mul_ps(a.re, _mm256_set1_ps(b)), \
                                     _mm256_mul_ps(a.im, _mm256_set1_ps(b))}

#define C8_DIV_C8(a, b) (complex_f8){_mm256_div_ps(_mm256_sub_ps(_mm256_mul_ps(a.re, b.re), _mm256_mul_ps(a.im, b.im)), \
                                                   _mm256_add_ps(_mm256_mul_ps(b.re, b.re), _mm256_mul_ps(b.im, b.im))), \
                                     _mm256_div_ps(_mm256_add_ps(_mm256_mul_ps(a.re, b.im), _mm256_mul_ps(a.im, b.re)), \
                                                   _mm256_add_ps(_mm256_mul_ps(b.re, b.re), _mm256_mul_ps(b.im, b.im)))}

#define C8_DIV_F8(a, b) (complex_f8){_mm256_div_ps(a.re, b), \
                                     _mm256_div_ps(a.im, b)}

#define C8_DIV_F1(a, b) (complex_f8){_mm256_div_ps(a.re, _mm256_set1_ps(b)), \
                                     _mm256_div_ps(a.im, _mm256_set1_ps(b))}

#define C1_EXP_I_F1(a) (complex_f1){cosf(a), sinf(a)}

void *avx_malloc(size_t size);

void avx_free(void *ptr);

#endif //ALGORITHM_AVX_H
