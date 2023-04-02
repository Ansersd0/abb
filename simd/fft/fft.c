//
// Created by wcy on 3/31/23.
//
#include <math.h>
#include <string.h>

#include "fft.h"

#define C8_F16(a0_re, a0_im, a1_re, a1_im, a2_re, a2_im, a3_re, a3_im, a4_re, a4_im, a5_re, a5_im, a6_re, a6_im, a7_re, a7_im) \
    (complex_f8) { \
        .re = _mm256_setr_ps(a0_re, a1_re, a2_re, a3_re, a4_re, a5_re, a6_re, a7_re), \
        .im = _mm256_setr_ps(a0_im, a1_im, a2_im, a3_im, a4_im, a5_im, a6_im, a7_im) \
    }                                                                                                                         \

#define C8_8C1(c1_0, c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7) \
    (complex_f8) { \
        .re = _mm256_setr_ps(c1_0.re, c1_1.re, c1_2.re, c1_3.re, c1_4.re, c1_5.re, c1_6.re, c1_7.re), \
        .im = _mm256_setr_ps(c1_0.im, c1_1.im, c1_2.im, c1_3.im, c1_4.im, c1_5.im, c1_6.im, c1_7.im) \
    }

#define C8_ADD_C8(a, b) \
    (complex_f8) { \
        .re = _mm256_add_ps(a.re, b.re), \
        .im = _mm256_add_ps(a.im, b.im) \
    }

#define C8_SUB_C8(a, b) \
    (complex_f8) { \
        .re = _mm256_sub_ps(a.re, b.re), \
        .im = _mm256_sub_ps(a.im, b.im) \
    }

#define C8_MUL_C8(a, b) \
    (complex_f8) { \
        .re = _mm256_sub_ps(_mm256_mul_ps(a.re, b.re), _mm256_mul_ps(a.im, b.im)), \
        .im = _mm256_add_ps(_mm256_mul_ps(a.re, b.im), _mm256_mul_ps(a.im, b.re)) \
    }

#define C8_MUL_C8_CONJ(a, b) \
    (complex_f8) { \
        .re = _mm256_add_ps(_mm256_mul_ps(a.re, b.re), _mm256_mul_ps(a.im, b.im)), \
        .im = _mm256_sub_ps(_mm256_mul_ps(a.im, b.re), _mm256_mul_ps(a.re, b.im)) \
    }

#define C8_MUL_I(a) \
    (complex_f8) { \
        .re = _mm256_sub_ps(_mm256_setzero_ps(), a.im), \
        .im = a.re \
    }

#define C8_MUL_C1(a, b) \
    (complex_f8) { \
        .re = _mm256_sub_ps(_mm256_mul_ps(a.re, _mm256_set1_ps(b.re)), _mm256_mul_ps(a.im, _mm256_set1_ps(b.im))), \
        .im = _mm256_add_ps(_mm256_mul_ps(a.re, _mm256_set1_ps(b.im)), _mm256_mul_ps(a.im, _mm256_set1_ps(b.re))) \
    }

#define C8_MUL_C1_CONJ(a, b) \
    (complex_f8) { \
        .re = _mm256_add_ps(_mm256_mul_ps(a.re, _mm256_set1_ps(b.re)), _mm256_mul_ps(a.im, _mm256_set1_ps(b.im))), \
        .im = _mm256_sub_ps(_mm256_mul_ps(a.im, _mm256_set1_ps(b.re)), _mm256_mul_ps(a.re, _mm256_set1_ps(b.im))) \
    }

#define C1_F1(a_re, a_im) \
    (complex_f1) { \
        .re = a_re, \
        .im = a_im \
    }

#define C1_EXP_I_F1(a) \
    (complex_f1) { \
        .re = cosf(a), \
        .im = sinf(a) \
    }

#define C1_CONJ(a) \
    (complex_f1) { \
        .re = a.re, \
        .im = -a.im \
    }


void *avx_malloc(size_t size) {
    return aligned_alloc(32, size);
}

void avx_free(void *ptr) {
    free(ptr);
}


static int decompose(int N, int *p) {
    int n = 0;

    while (N % 4 == 0) {
        N /= 4;
        p[n++] = 4;
    }

    while (N % 2 == 0) {
        N /= 2;
        p[n++] = 2;
    }
    for (int i = 3; i * i <= N; i += 2) {
        while (N % i == 0) {
            N /= i;
            p[n++] = i;
        }
    }
    if (N > 1) {
        p[n++] = N;
    }
    return n;
}

static inline void avx_fft_avx_8(complex_f8 *__restrict in,
                                 complex_f8 *__restrict out,
                                 complex_f8 *__restrict W,
                                 complex_f8 *__restrict temp[8],
                                 int N, float type) {

    float *re, *im, *t_re, *t_im;

    int n = N / 8;
    int m = N % 8;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < n; ++j) {
            t_re = (float *) &temp[i][j].re;
            t_im = (float *) &temp[i][j].im;

            for (int k = 0; k < 8; ++k) {
                t_re[k] = in[(i * N + j * 8 + k) / 8].re[(i * N + j * 8 + k) % 8];
                t_im[k] = in[(i * N + j * 8 + k) / 8].im[(i * N + j * 8 + k) % 8];
            }
        }

        t_re = (float *) &temp[i][n].re;
        t_im = (float *) &temp[i][n].im;
        for (int j = 0; j < m; ++j) {
            t_re[j] = in[(i * N + n * 8 + j) / 8].re[(i * N + n * 8 + j) % 8];
            t_im[j] = in[(i * N + n * 8 + j) / 8].im[(i * N + n * 8 + j) % 8];
        }
    }

    complex_f8 aa[8], bb[8];
    complex_f1 type_1_8 = C1_F1(0.70710678118654752440084436210484903928483593768847f,
                                type * 0.70710678118654752440084436210484903928483593768847f);
    complex_f1 type_I = C1_F1(0, type);
    complex_f1 type_3_8 = C1_F1(-0.70710678118654752440084436210484903928483593768847f,
                                type * 0.70710678118654752440084436210484903928483593768847f);
    for (int i = 0; i < n; ++i) {
        aa[0] = C8_ADD_C8(temp[0][i], temp[4][i]);
        aa[1] = C8_SUB_C8(temp[0][i], temp[4][i]);
        aa[2] = C8_ADD_C8(temp[1][i], temp[5][i]);
        aa[3] = C8_MUL_C1(C8_SUB_C8(temp[1][i], temp[5][i]), type_1_8);
        aa[4] = C8_ADD_C8(temp[2][i], temp[6][i]);
        aa[5] = C8_MUL_C1(C8_SUB_C8(temp[2][i], temp[6][i]), type_I);
        aa[6] = C8_ADD_C8(temp[3][i], temp[7][i]);
        aa[7] = C8_MUL_C1(C8_SUB_C8(temp[3][i], temp[7][i]), type_3_8);

        bb[0] = C8_ADD_C8(aa[0], aa[4]);
        bb[1] = C8_ADD_C8(aa[1], aa[5]);
        bb[2] = C8_SUB_C8(aa[0], aa[4]);
        bb[3] = C8_SUB_C8(aa[1], aa[5]);
        bb[4] = C8_ADD_C8(aa[2], aa[6]);
        bb[5] = C8_ADD_C8(aa[3], aa[7]);
        bb[6] = C8_MUL_C1(C8_SUB_C8(aa[2], aa[6]), type_I);
        bb[7] = C8_MUL_C1(C8_SUB_C8(aa[3], aa[7]), type_I);

        if (N == 1) {
            aa[0] = C8_ADD_C8(bb[0], bb[4]);
            aa[1] = C8_ADD_C8(bb[1], bb[5]);
            aa[2] = C8_ADD_C8(bb[2], bb[6]);
            aa[3] = C8_ADD_C8(bb[3], bb[7]);
            aa[4] = C8_SUB_C8(bb[0], bb[4]);
            aa[5] = C8_SUB_C8(bb[1], bb[5]);
            aa[6] = C8_SUB_C8(bb[2], bb[6]);
            aa[7] = C8_SUB_C8(bb[3], bb[7]);
        } else if (type == -1.0f) {
            aa[0] = C8_MUL_C8(C8_ADD_C8(bb[0], bb[4]), W[8 * i]);
            aa[1] = C8_MUL_C8(C8_ADD_C8(bb[1], bb[5]), W[8 * i + 1]);
            aa[2] = C8_MUL_C8(C8_ADD_C8(bb[2], bb[6]), W[8 * i + 2]);
            aa[3] = C8_MUL_C8(C8_ADD_C8(bb[3], bb[7]), W[8 * i + 3]);
            aa[4] = C8_MUL_C8(C8_SUB_C8(bb[0], bb[4]), W[8 * i + 4]);
            aa[5] = C8_MUL_C8(C8_SUB_C8(bb[1], bb[5]), W[8 * i + 5]);
            aa[6] = C8_MUL_C8(C8_SUB_C8(bb[2], bb[6]), W[8 * i + 6]);
            aa[7] = C8_MUL_C8(C8_SUB_C8(bb[3], bb[7]), W[8 * i + 7]);
        } else {
            aa[0] = C8_MUL_C8_CONJ(C8_ADD_C8(bb[0], bb[4]), W[8 * i]);
            aa[1] = C8_MUL_C8_CONJ(C8_ADD_C8(bb[1], bb[5]), W[8 * i + 1]);
            aa[2] = C8_MUL_C8_CONJ(C8_ADD_C8(bb[2], bb[6]), W[8 * i + 2]);
            aa[3] = C8_MUL_C8_CONJ(C8_ADD_C8(bb[3], bb[7]), W[8 * i + 3]);
            aa[4] = C8_MUL_C8_CONJ(C8_SUB_C8(bb[0], bb[4]), W[8 * i + 4]);
            aa[5] = C8_MUL_C8_CONJ(C8_SUB_C8(bb[1], bb[5]), W[8 * i + 5]);
            aa[6] = C8_MUL_C8_CONJ(C8_SUB_C8(bb[2], bb[6]), W[8 * i + 6]);
            aa[7] = C8_MUL_C8_CONJ(C8_SUB_C8(bb[3], bb[7]), W[8 * i + 7]);
        }

        if (m == 0) {
            for (int j = 0; j < 8; ++j) {
                re = (float *) &out[i * 8 + j].re;
                im = (float *) &out[i * 8 + j].im;
                for (int k = 0; k < 8; ++k) {
                    re[k] = aa[k].re[j];
                    im[k] = aa[k].im[j];
                }
            }
        } else if (i < n - 1) {
            for (int j = 0; j < 8; ++j) {
                re = (float *) &out[i * 8 + j].re;
                im = (float *) &out[i * 8 + j].im;
                for (int k = 0; k < 8; ++k) {
                    re[k] = aa[k].re[j];
                    im[k] = aa[k].im[j];
                }
            }
        } else {
            for (int j = 0; j < m; ++j) {
                re = (float *) &out[i * 8 + j].re;
                im = (float *) &out[i * 8 + j].im;
                for (int k = 0; k < 8; ++k) {
                    re[k] = aa[k].re[j];
                    im[k] = aa[k].im[j];
                }
            }
        }
    }

}


static inline void avx_fft_paas_4(complex_f8 *__restrict in,
                                  complex_f8 *__restrict out,
                                  complex_f1 *__restrict W,
                                  int P1, int P2, int P3, float type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    complex_f8 alpha_0, alpha_1, beta_0, beta_1;

    if (P3 == 1) {
        for (int a = 0; a < P1; ++a) {
            alpha_0 = C8_ADD_C8(in[a], in[2 * P13 + a]);
            alpha_1 = C8_SUB_C8(in[a], in[2 * P13 + a]);
            beta_0 = C8_ADD_C8(in[1 * P13 + a], in[3 * P13 + a]);
            beta_1 = C8_MUL_C1(C8_SUB_C8(in[1 * P13 + a], in[3 * P13 + a]), C1_F1(0, type));
            out[0 * P1 + a] = C8_ADD_C8(alpha_0, beta_0);
            out[1 * P1 + a] = C8_ADD_C8(alpha_1, beta_1);
            out[2 * P1 + a] = C8_SUB_C8(alpha_0, beta_0);
            out[3 * P1 + a] = C8_SUB_C8(alpha_1, beta_1);
        }
    } else {
        if (type == -1) {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    alpha_0 = C8_ADD_C8(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    alpha_1 = C8_SUB_C8(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    beta_0 = C8_ADD_C8(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]);
                    beta_1 = C8_MUL_C1(C8_SUB_C8(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]), C1_F1(0, type));

                    out[c * P12 + a] = C8_ADD_C8(alpha_0, beta_0);
                    out[c * P12 + 1 * P1 + a] = C8_MUL_C1(C8_ADD_C8(alpha_1, beta_1), W[c * (P2 - 1)]);
                    out[c * P12 + 2 * P1 + a] = C8_MUL_C1(C8_SUB_C8(alpha_0, beta_0), W[c * (P2 - 1) + 1]);
                    out[c * P12 + 3 * P1 + a] = C8_MUL_C1(C8_SUB_C8(alpha_1, beta_1), W[c * (P2 - 1) + 2]);
                }
            }
        } else {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    alpha_0 = C8_ADD_C8(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    alpha_1 = C8_SUB_C8(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    beta_0 = C8_ADD_C8(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]);
                    beta_1 = C8_MUL_C1(C8_SUB_C8(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]), C1_F1(0, type));

                    out[c * P12 + +a] = C8_ADD_C8(alpha_0, beta_0);
                    out[c * P12 + 1 * P1 + a] = C8_MUL_C1_CONJ(C8_ADD_C8(alpha_1, beta_1), W[c * (P2 - 1)]);
                    out[c * P12 + 2 * P1 + a] = C8_MUL_C1_CONJ(C8_SUB_C8(alpha_0, beta_0), W[c * (P2 - 1) + 1]);
                    out[c * P12 + 3 * P1 + a] = C8_MUL_C1_CONJ(C8_SUB_C8(alpha_1, beta_1), W[c * (P2 - 1) + 2]);
                }
            }
        }

    }
}

static inline void avx_fft_paas_2(complex_f8 *__restrict in,
                                  complex_f8 *__restrict out,
                                  complex_f1 *__restrict W,
                                  int P1, int P2, int P3, float type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    if (P3 == 1) {
        for (int a = 0; a < P1; ++a) {
            out[0 * P1 + a] = C8_ADD_C8(in[a], in[P13 + a]);
            out[1 * P1 + a] = C8_SUB_C8(in[a], in[P13 + a]);
        }
    } else {
        if (type == -1) {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + 0 * P1 + a] = C8_ADD_C8(in[c * P1 + a], in[P13 + c * P1 + a]);
                    out[c * P12 + 1 * P1 + a] = C8_MUL_C1(C8_SUB_C8(in[c * P1 + a], in[P13 + c * P1 + a]),
                                                          W[c * (P2 - 1)]);
                }
            }
        } else {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + 0 * P1 + a] = C8_ADD_C8(in[c * P1 + a], in[P13 + c * P1 + a]);
                    out[c * P12 + 1 * P1 + a] = C8_MUL_C1_CONJ(C8_SUB_C8(in[c * P1 + a], in[P13 + c * P1 + a]),
                                                               W[c * (P2 - 1)]);
                }
            }
        }
    }
}

static inline void avx_fft_paas_n(complex_f8 *__restrict in,
                                  complex_f8 *__restrict out,
                                  complex_f1 *__restrict W,
                                  int P1, int P2, int P3, int type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    if (P3 == 1) {
        for (int b = 0; b < P2; ++b) {
            for (int a = 0; a < P1; ++a) {
                out[b * P1 + a] = in[a];
                for (int i = 1; i < P2; ++i) {
                    out[b * P1 + a] = C8_ADD_C8(out[b * P1 + a],
                                                C8_MUL_C1(in[i * P13 + a],
                                                          C1_EXP_I_F1(type * 2 * M_PI * (b * i) / P2)));
                }
            }
        }
    } else {
        for (int c = 0; c < P3; ++c) {
            for (int b = 0; b < P2; ++b) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + b * P1 + a] = in[c * P1 + a];
                    for (int i = 1; i < P2; ++i) {
                        out[c * P12 + b * P1 + a] = C8_ADD_C8(out[c * P12 + b * P1 + a],
                                                              C8_MUL_C1(in[i * P13 + c * P1 + a],
                                                                        C1_EXP_I_F1(type * 2 * M_PI * (b * i) / P2)));
                    }
                    out[c * P12 + b * P1 + a] = C8_MUL_C1(out[c * P12 + b * P1 + a],
                                                          C1_EXP_I_F1(type * 2 * M_PI * (b * c) / (P2 * P3)));
                }
            }
        }
    }
}

AVX_FFT *avx_fft_new(complex_f8 **time_domain, complex_f8 **freq_domain, fft_type N) {
    AVX_FFT *fft = malloc(sizeof(AVX_FFT));
    fft->N = N;
    fft->P_n = decompose(N, fft->P);
    fft->W = malloc(sizeof(complex_f1 *) * fft->P_n);
    fft->time_domain = time_domain;
    fft->freq_domain = freq_domain;
    int P3 = N;
    complex_f1 *W_ip;
    for (int i = 0; i < fft->P_n; ++i) {
        P3 /= fft->P[i];
        fft->W[i] = malloc(sizeof(complex_f1) * (fft->P[i] - 1) * P3);
        W_ip = fft->W[i];
        for (int c = 0; c < P3; ++c) {
            for (int b = 1; b < fft->P[i]; ++b) {
                *W_ip++ = C1_EXP_I_F1(-2 * M_PI * b * c / (fft->P[i] * P3));
            }
        }
    }

    int n = N / 8;
    int m = N / 8;
    float *re, *im;
    if (N % 8 == 0) {
        for (int i = 0; i < 8; ++i) {
            fft->temp[i] = avx_malloc(sizeof(complex_f8) * n * 8);
            memset(fft->temp, 0, sizeof(complex_f8) * n * 8);
        }
        fft->W8 = avx_malloc(sizeof(complex_f8) * n * 8);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 8; ++j) {
                re = (float *) &fft->W8[i * 8 + j].re;
                im = (float *) &fft->W8[i * 8 + j].im;
                for (int k = 0; k < 8; ++k) {
                    re[k] = cosf(-2 * M_PI * (j + 8 * i) * k / (8 * N));
                    im[k] = sinf(-2 * M_PI * (j + 8 * i) * k / (8 * N));
                }
            }
        }
    } else {
        for (int i = 0; i < 8; ++i) {
            fft->temp[i] = avx_malloc(sizeof(complex_f8) * (n + 1) * 8);
            memset(fft->temp, 0, sizeof(complex_f8) * (n + 1) * 8);
        }
        fft->W8 = avx_malloc(sizeof(complex_f8) * (n + 1) * 8);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 8; ++j) {
                re = (float *) &fft->W8[i * 8 + j].re;
                im = (float *) &fft->W8[i * 8 + j].im;
                for (int k = 0; k < 8; ++k) {
                    re[k] = cosf(-2 * M_PI * (j + 8 * i) * k / (8 * N));
                    im[k] = sinf(-2 * M_PI * (j + 8 * i) * k / (8 * N));
                }
            }
        }
        for (int j = 0; j < m; ++j) {
            re = (float *) &fft->W8[n * 8 + j].re;
            im = (float *) &fft->W8[n * 8 + j].im;
            for (int k = 0; k < 8; ++k) {
                re[k] = cosf(-2 * M_PI * (j + 8 * n) * k / (8 * N));
                im[k] = sinf(-2 * M_PI * (j + 8 * n) * k / (8 * N));
            }
        }
    }


    return fft;
}

void avx_fft(AVX_FFT *fft, fft_type type) {
    bool flag = 0;
    complex_f8 *p;

    complex_f8 *in = type == Forward ? *fft->time_domain : *fft->freq_domain;
    complex_f8 *out = type == Forward ? *fft->freq_domain : *fft->time_domain;

//    avx_fft_avx_8(in, out, fft->W8, fft->temp, fft->N, type);

    int P1 = 1, P2 = 1, P3 = fft->N;
    for (int ii = 0; ii < fft->P_n; ++ii) {
        P1 *= P2;
        P2 = fft->P[ii];
        P3 /= P2;

        switch (P2) {
            case 4:
                avx_fft_paas_4(in, out, fft->W[ii], P1, P2, P3, type);
                break;
            case 2:
                avx_fft_paas_2(in, out, fft->W[ii], P1, P2, P3, type);
                break;
            default:
                avx_fft_paas_n(in, out, fft->W[ii], P1, P2, P3, type);
                break;
        }

        flag = !flag;
        p = in;
        in = out;
        out = p;
    }

    if (!flag) {
        *fft->time_domain = type == Forward ? out : in;
        *fft->freq_domain = type == Forward ? in : out;
    }
}


void avx_fft_free(AVX_FFT *fft) {
    for (int i = 0; i < fft->P_n; ++i) {
        free(fft->W[i]);
    }
    free(fft->W);
    free(fft);
}