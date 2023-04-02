//
// Created by wcy on 2023/2/20.
//

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "fft.h"


FFT *fft_new(int N) {
    if (N % VEC_SIZE != 0) {
        printf("N must be a multiple of %d", VEC_SIZE);
        return NULL;
    }
    FFT *fft = malloc(sizeof(FFT));
    fft->N_float = N;
    fft->m_cev = N / VEC_SIZE;

    fft->P_n = decompose(fft->m_cev, fft->P);

    fft->W = malloc(sizeof(vec8f) * N);
    float *W_p = fft->W->f;
    for (int i = 0; i < fft->P_n; ++i) {
        *W_p = 1;
        *W_p = cosf(2 * M_PI / fft->P[i]);
    }
    return fft;
}


void fft_f_4_1(float _Complex *in, float _Complex *out, int N) {
    if (N / 4 != 0) {
        return;
    }
    vec4c *p_in = (vec4c *) in;
    vec4c *p_out = (vec4c *) out;
    vec4c temp;

    int block_len = N / 4;
    int block_vec4_n = block_len / VEC_SIZE;
    int block_vec4_remain = block_len % VEC_SIZE;
    for (int i = 0; i < VEC_SIZE; ++i) {
        temp.v = p_in[i].v;
        for (int j = 1; j < block_vec4_n; ++j) {

        }
        p_out[i] = temp;
    }

    // N > 64


    int N_vec = N / 8;
    int P[10];
    int P_n = decompose(N_vec, P);
//    vec8f temp;
    for (int i = 0; i < P_n; ++i) {
        int m_vec = N_vec / P[i];
        for (int j = 0; j < m_vec; ++j) {
            temp.v = p_in[j].v;
            vec8f W[P[i]];
            for (int k = 1; k < P[i]; ++k) {
                W[k].v = _mm256_setr_ps(cosf(2 * M_PI * k / P[i]),
                                        sinf(2 * M_PI * k / P[i]),
                                        cosf(2 * M_PI * k / P[i]),
                                        sinf(2 * M_PI * k / P[i]),
                                        cosf(2 * M_PI * k / P[i]),
                                        sinf(2 * M_PI * k / P[i]),
                                        cosf(2 * M_PI * k / P[i]),
                                        sinf(2 * M_PI * k / P[i]));
            }
            for (int k = 1; k < P[i]; ++k) {
                temp.v += p_in[j + k * m_vec].v * W[k].v;
            }

            p_out[j] = temp;
        }
        p_in = p_out;
    }
}

int decompose(int N, int *p) {
    int n = 0;

//    while (N % 4 == 0) {
//        N /= 4;
//        p[n++] = 4;
//    }

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

typedef enum {
    backward = 1,
    forward = -1
} fft_type;

void my_fft(double _Complex *in, double _Complex *out, int N, fft_type type) {
    int P[20];
    int n = decompose(N, P);

    int flag = 0;
    double _Complex *p;

    int P1 = 1, P2 = 1, P3 = N;
    for (int ii = 0; ii < n; ++ii) {
        P1 *= P2;
        P2 = P[ii];
        P3 /= P2;

        int P13 = P1 * P3;
        int P12 = P1 * P2;

        for (int c = 0; c < P3; ++c) {
            for (int b = 0; b < P2; ++b) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + b * P1 + a] = in[c * P1 + a];
                    for (int i = 1; i < P2; ++i) {
                        out[c * P12 + b * P1 + a]
                                += in[i * P13 + c * P1 + a]
                                   * cexp(type * 2 * M_PI * _Complex_I * (b * i) / P2);
                    }
                    out[c * P12 + b * P1 + a] *= cexp(type * 2 * M_PI * _Complex_I * (b * c) / (P2 * P3));
                }
            }
        }

        flag = !flag;
        p = in;
        in = out;
        out = p;
    }

    if (!flag) {
        memcpy(out, in, sizeof(double _Complex) * N);
    }
}

void fft_test() {
    int N = 7158;
    double _Complex *in, *out;
    in = malloc(sizeof(double _Complex) * N);
    out = malloc(sizeof(double _Complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i] = 1 + i + (2 + i) * _Complex_I;
    }

    my_fft(in, out, N, forward);
    my_fft(out, in, N, backward);

    for (int i = 0; i < N % 51; ++i) {
        printf("in[%d] = %lf + %lf i\n", i, (__real__ in[i]) / N, cimag(in[i]) / N);
    }

    free(in);
    free(out);
}

FFT_1_HEADER *fft_1_new(int N) {
    FFT_1_HEADER *fft = malloc(sizeof(FFT_1_HEADER));
    fft->N = N;
    fft->P_n = decompose(N, fft->P);
    fft->W = malloc(sizeof(float _Complex *) * fft->P_n);
    int P3 = N;
    float _Complex *W_ip;
    for (int i = 0; i < fft->P_n; ++i) {
        P3 /= fft->P[i];
        fft->W[i] = malloc(sizeof(float _Complex) * (fft->P[i] - 1) * P3);
        W_ip = fft->W[i];
        for (int c = 0; c < P3; ++c) {
            for (int b = 1; b < fft->P[i]; ++b) {
                *W_ip++ = cexpf(-2 * M_PI * _Complex_I * (b * c) / (fft->P[i] * P3));
            }
        }
    }
    return fft;
}

void fft_1_free(FFT_1_HEADER *fft) {
    for (int i = 0; i < fft->P_n; ++i) {
        free(fft->W[i]);
    }
    free(fft->W);
    free(fft);
}

void fft_1(FFT_1_HEADER *fft, float _Complex *__restrict__ in, float _Complex *__restrict__ out, int type) {
    int flag = 0;
    float _Complex *p;

    int P1 = 1, P2 = 1, P3 = fft->N;
    for (int ii = 0; ii < fft->P_n; ++ii) {
        P1 *= P2;
        P2 = fft->P[ii];
        P3 /= P2;

        switch (P2) {
            case 4:
                fft_1_paas_4(in, out, fft->W[ii], P1, P2, P3, type);
                break;
            case 2:
                fft_1_paas_2(in, out, fft->W[ii], P1, P2, P3, type);
                break;
            default:
                fft_1_paas_n(in, out, fft->W[ii], P1, P2, P3, type);
                break;
        }

        flag = !flag;
        p = in;
        in = out;
        out = p;
    }

    if (!flag) {
        memcpy(out, in, sizeof(float _Complex) * fft->N);
    }
}

inline void fft_1_paas_4(float _Complex *__restrict__ in,
                         float _Complex *__restrict__ out,
                         float _Complex *__restrict__ W,
                         int P1, int P2, int P3, int type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    float _Complex alpha_0, alpha_1, beta_0, beta_1;

    if (P3 == 1) {
        for (int a = 0; a < P1; ++a) {
            alpha_0 = in[a] + in[2 * P13 + a], alpha_1 = in[a] - in[2 * P13 + a];
            beta_0 = in[1 * P13 + a] + in[3 * P13 + a], beta_1 = in[1 * P13 + a] - in[3 * P13 + a];
            out[0 * P1 + a] = alpha_0 + beta_0;
            out[1 * P1 + a] = alpha_1 + type * beta_1 * _Complex_I;
            out[2 * P1 + a] = alpha_0 - beta_0;
            out[3 * P1 + a] = alpha_1 - type * beta_1 * _Complex_I;
        }
    } else {
        if (type == -1) {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    alpha_0 = in[c * P1 + a] + in[2 * P13 + c * P1 + a];
                    alpha_1 = in[c * P1 + a] - in[2 * P13 + c * P1 + a];
                    beta_0 = in[1 * P13 + c * P1 + a] + in[3 * P13 + c * P1 + a];
                    beta_1 = in[1 * P13 + c * P1 + a] - in[3 * P13 + c * P1 + a];

                    out[c * P12 + +a] = (alpha_0 + beta_0);
                    out[c * P12 + 1 * P1 + a] = (alpha_1 + type * beta_1 * _Complex_I) * W[c * (P2 - 1)];
                    out[c * P12 + 2 * P1 + a] = (alpha_0 - beta_0) * W[c * (P2 - 1) + 1];
                    out[c * P12 + 3 * P1 + a] = (alpha_1 - type * beta_1 * _Complex_I) * W[c * (P2 - 1) + 2];
                }
            }
        } else {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    alpha_0 = in[c * P1 + a] + in[2 * P13 + c * P1 + a];
                    alpha_1 = in[c * P1 + a] - in[2 * P13 + c * P1 + a];
                    beta_0 = in[1 * P13 + c * P1 + a] + in[3 * P13 + c * P1 + a];
                    beta_1 = in[1 * P13 + c * P1 + a] - in[3 * P13 + c * P1 + a];

                    out[c * P12 + +a] = (alpha_0 + beta_0);
                    out[c * P12 + 1 * P1 + a] = (alpha_1 + type * beta_1 * _Complex_I) * conjf(W[c * (P2 - 1)]);
                    out[c * P12 + 2 * P1 + a] = (alpha_0 - beta_0) * conjf(W[c * (P2 - 1) + 1]);
                    out[c * P12 + 3 * P1 + a] = (alpha_1 - type * beta_1 * _Complex_I) * conjf(W[c * (P2 - 1) + 2]);
                }
            }
        }
    }
}


inline void fft_1_paas_2(float _Complex *__restrict__ in,
                         float _Complex *__restrict__ out,
                         float _Complex *__restrict__ W,
                         int P1, int P2, int P3, int type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    if (P3 == 1) {
        for (int a = 0; a < P1; ++a) {
            out[0 * P1 + a] = in[a] + in[P13 + a];
            out[1 * P1 + a] = in[a] - in[P13 + a];
        }
    } else {
        if (type == -1) {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + 0 * P1 + a] = (in[c * P1 + a] + in[P13 + c * P1 + a]);
                    out[c * P12 + 1 * P1 + a] = (in[c * P1 + a] - in[P13 + c * P1 + a]) * W[c * (P2 - 1)];
                }
            }
        } else {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + 0 * P1 + a] = (in[c * P1 + a] + in[P13 + c * P1 + a]);
                    out[c * P12 + 1 * P1 + a] = (in[c * P1 + a] - in[P13 + c * P1 + a]) * conjf(W[c * (P2 - 1)]);
                }
            }
        }
    }
}

inline void fft_1_paas_n(float _Complex *__restrict__ in,
                         float _Complex *__restrict__ out,
                         float _Complex *__restrict__ W,
                         int P1, int P2, int P3, int type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    if (P3 == 1) {
        for (int b = 0; b < P2; ++b) {
            for (int a = 0; a < P1; ++a) {
                out[b * P1 + a] = in[a];
                for (int i = 1; i < P2; ++i) {
                    out[b * P1 + a] += in[i * P13 + a] * cexp(type * 2 * M_PI * _Complex_I * (b * i) / P2);
                }
            }
        }
    } else {
        for (int c = 0; c < P3; ++c) {
            for (int b = 0; b < P2; ++b) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + b * P1 + a] = in[c * P1 + a];
                    for (int i = 1; i < P2; ++i) {
                        out[c * P12 + b * P1 + a] +=
                                in[i * P13 + c * P1 + a] * cexp(type * 2 * M_PI * _Complex_I * (b * i) / P2);
                    }
                    out[c * P12 + b * P1 + a] *= cexp(type * 2 * M_PI * _Complex_I * (b * c) / (P2 * P3));
                }
            }
        }
    }
}

void fft_1_test() {
    int N = 10;
    float _Complex *in, *out;
    in = malloc(sizeof(float _Complex) * N);
    out = malloc(sizeof(float _Complex) * N);

    for (int i = 0; i < N; ++i) {
//        in[i] = 1 + i + (2 + i) * _Complex_I;
        in[i] = 1 + i;
    }

    FFT_1_HEADER *fft = fft_1_new(N);

    // get start time
    struct timeval start, end;
    gettimeofday(&start, NULL);

    fft_1(fft, in, out, forward);

    for (int i = 0; i < N; ++i) {
        printf("%d: %f + %f * i\n", i, creal(out[i]), cimag(out[i]));
    }

//    fft_1(fft, out, in, backward);

    // get end time
    gettimeofday(&end, NULL);
//
//    for (int i = 0; i < N % 51; ++i) {
//        printf("%d: %f + %f * i\n", i, creal(in[i]) / N, cimag(in[i]) / N);
//    }

    // calculate time
    double timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

    printf("time: %f ms\n", timeuse / 1000);
}