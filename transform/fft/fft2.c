//
// Created by wcy on 2023/2/20.
//

#include <math.h>
#include <malloc.h>
#include <stdbool.h>
#include "fft2.h"

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


static inline void fft_2_paas_4(f_complex *__restrict in,
                                f_complex *__restrict out,
                                f_complex *__restrict W,
                                int P1, int P2, int P3, int type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    f_complex alpha_0, alpha_1, beta_0, beta_1;

    if (P3 == 1) {
        for (int a = 0; a < P1; ++a) {
            alpha_0 = caddf(in[a], in[2 * P13 + a]), alpha_1 = csubf(in[a], in[2 * P13 + a]);
            beta_0 = caddf(in[1 * P13 + a], in[3 * P13 + a]), beta_1 = csubf(in[1 * P13 + a], in[3 * P13 + a]);
            out[0 * P1 + a] = caddf(alpha_0, beta_0);
            out[1 * P1 + a] = caddf(alpha_1, cmulf_I(type * beta_1));
            out[2 * P1 + a] = csubf(alpha_0, beta_0);
            out[3 * P1 + a] = csubf(alpha_1, cmulf_I(type * beta_1));
        }
    } else {
        if (type == -1) {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    alpha_0 = caddf(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    alpha_1 = csubf(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    beta_0 = caddf(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]);
                    beta_1 = csubf(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]);

                    out[c * P12 + a] = caddf(alpha_0, beta_0);
                    out[c * P12 + 1 * P1 + a] = cmulf(caddf(alpha_1, cmulf_I(type * beta_1)), W[c * (P2 - 1)]);
                    out[c * P12 + 2 * P1 + a] = cmulf(csubf(alpha_0, beta_0), W[c * (P2 - 1) + 1]);
                    out[c * P12 + 3 * P1 + a] = cmulf(csubf(alpha_1, cmulf_I(type * beta_1)),
                                                      W[c * (P2 - 1) + 2]);
                }
            }
        } else {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    alpha_0 = caddf(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    alpha_1 = csubf(in[c * P1 + a], in[2 * P13 + c * P1 + a]);
                    beta_0 = caddf(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]);
                    beta_1 = csubf(in[1 * P13 + c * P1 + a], in[3 * P13 + c * P1 + a]);

                    out[c * P12 + +a] = caddf(alpha_0, beta_0);
                    out[c * P12 + 1 * P1 + a] = cmulf(caddf(alpha_1, cmulf_I(type * beta_1)),
                                                      cconjf(W[c * (P2 - 1)]));
                    out[c * P12 + 2 * P1 + a] = cmulf(csubf(alpha_0, beta_0), cconjf(W[c * (P2 - 1) + 1]));
                    out[c * P12 + 3 * P1 + a] = cmulf(csubf(alpha_1, cmulf_I(type * beta_1)),
                                                      cconjf(W[c * (P2 - 1) + 2]));
                }
            }
        }

    }
}

static inline void fft_2_paas_2(f_complex *__restrict in,
                                f_complex *__restrict out,
                                f_complex *__restrict W,
                                int P1, int P2, int P3, int type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    if (P3 == 1) {
        for (int a = 0; a < P1; ++a) {
            out[0 * P1 + a] = caddf(in[a], in[P13 + a]);
            out[1 * P1 + a] = csubf(in[a], in[P13 + a]);
        }
    } else {
        if (type == -1) {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + 0 * P1 + a] = caddf(in[c * P1 + a], in[P13 + c * P1 + a]);
                    out[c * P12 + 1 * P1 + a] = cmulf(csubf(in[c * P1 + a], in[P13 + c * P1 + a]), W[c * (P2 - 1)]);
                }
            }
        } else {
            for (int c = 0; c < P3; ++c) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + 0 * P1 + a] = caddf(in[c * P1 + a], in[P13 + c * P1 + a]);
                    out[c * P12 + 1 * P1 + a] = cmulf(csubf(in[c * P1 + a], in[P13 + c * P1 + a]),
                                                      cconjf(W[c * (P2 - 1)]));
                }
            }
        }
    }
}

static inline void fft_2_paas_n(f_complex *__restrict in,
                                f_complex *__restrict out,
                                f_complex *__restrict W,
                                int P1, int P2, int P3, int type) {
    int P13 = P1 * P3;
    int P12 = P1 * P2;

    if (P3 == 1) {
        for (int b = 0; b < P2; ++b) {
            for (int a = 0; a < P1; ++a) {
                out[b * P1 + a] = in[a];
                for (int i = 1; i < P2; ++i) {
                    out[b * P1 + a] = caddf(out[b * P1 + a],
                                            cmulf(in[i * P13 + a], cexpf(smulf_I(type * 2 * M_PI * (b * i) / P2))));
                }
            }
        }
    } else {
        for (int c = 0; c < P3; ++c) {
            for (int b = 0; b < P2; ++b) {
                for (int a = 0; a < P1; ++a) {
                    out[c * P12 + b * P1 + a] = in[c * P1 + a];
                    for (int i = 1; i < P2; ++i) {
                        out[c * P12 + b * P1 + a] = caddf(out[c * P12 + b * P1 + a],
                                                          cmulf(in[i * P13 + c * P1 + a],
                                                                cexpf(smulf_I(type * 2 * M_PI * (b * i) / P2))));
                    }
                    out[c * P12 + b * P1 + a] = cmulf(out[c * P12 + b * P1 + a],
                                                      cexpf(smulf_I(type * 2 * M_PI * (b * c) / (P2 * P3))));
                }
            }
        }
    }
}

FFT_2_HEADER *fft_2_new(f_complex **time_domain, f_complex **freq_domain, fft_type N) {
    FFT_2_HEADER *fft = malloc(sizeof(FFT_2_HEADER));
    fft->N = N;
    fft->P_n = decompose(N, fft->P);
    fft->W = malloc(sizeof(f_complex *) * fft->P_n);
    fft->time_domain = time_domain;
    fft->freq_domain = freq_domain;
    int P3 = N;
    f_complex *W_ip;
    for (int i = 0; i < fft->P_n; ++i) {
        P3 /= fft->P[i];
        fft->W[i] = malloc(sizeof(f_complex) * (fft->P[i] - 1) * P3);
        W_ip = fft->W[i];
        for (int c = 0; c < P3; ++c) {
            for (int b = 1; b < fft->P[i]; ++b) {
                *W_ip++ = cexpf(smulf_I(-2 * M_PI * b * c / (fft->P[i] * P3)));
            }
        }
    }
    return fft;
}

void fft_2(FFT_2_HEADER *fft, fft_type type) {
    bool flag = 0;
    f_complex *p;

    f_complex *in = type == Forward ? *fft->time_domain : *fft->freq_domain;
    f_complex *out = type == Forward ? *fft->freq_domain : *fft->time_domain;
    int P1 = 1, P2 = 1, P3 = fft->N;
    for (int ii = 0; ii < fft->P_n; ++ii) {
        P1 *= P2;
        P2 = fft->P[ii];
        P3 /= P2;

        switch (P2) {
            case 4:
                fft_2_paas_4(in, out, fft->W[ii], P1, P2, P3, type);
                break;
            case 2:
                fft_2_paas_2(in, out, fft->W[ii], P1, P2, P3, type);
                break;
            default:
                fft_2_paas_n(in, out, fft->W[ii], P1, P2, P3, type);
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


void fft_2_free(FFT_2_HEADER *fft) {
    for (int i = 0; i < fft->P_n; ++i) {
        free(fft->W[i]);
    }
    free(fft->W);
    free(fft);
}

