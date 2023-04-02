//
// Created by wcy on 2023/3/15.
//

#include <string.h>
#include <malloc.h>
#include <stdio.h>
#include "fft1.h"

#if defined(__GNUC__) || defined(__clang__)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

    my_fft(in, out, N, Forward);
    my_fft(out, in, N, Inverse);

    for (int i = 0; i < N % 51; ++i) {
        printf("in[%d] = %lf + %lf i\n", i, (__real__ in[i]) / N, cimag(in[i]) / N);
    }

    free(in);
    free(out);
}

FFT_2_HEADER *fft_1_new(float _Complex **in, float _Complex **out, int N) {
    FFT_2_HEADER *fft = malloc(sizeof(FFT_2_HEADER));
    fft->N = N;
    fft->P_n = decompose(N, fft->P);
    fft->W = malloc(sizeof(float _Complex *) * fft->P_n);
    fft->time_domain = in;
    fft->freq_domain = out;
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

void fft_1_free(FFT_2_HEADER *fft) {
    for (int i = 0; i < fft->P_n; ++i) {
        free(fft->W[i]);
    }
    free(fft->W);
    free(fft);
}

static inline void fft_1_paas_4(float _Complex *__restrict__ in,
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

static inline void fft_1_paas_2(float _Complex *__restrict__ in,
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


static inline void fft_1_paas_n(float _Complex *__restrict__ in,
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

void fft_1(FFT_2_HEADER *fft, fft_type type) {
    int flag = 0;
    float _Complex *p;
    float _Complex *in = type == Forward ? *fft->time_domain : *fft->freq_domain;
    float _Complex *out = type == Forward ? *fft->freq_domain : *fft->time_domain;

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
        *fft->freq_domain = type == Forward ? in : out;
        *fft->time_domain = type == Forward ? out : in;
    }
}

#endif