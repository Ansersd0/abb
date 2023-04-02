//
// Created by wcy on 2023/2/27.
//

#ifndef ALGORITHM_REAL_FFT_H
#define ALGORITHM_REAL_FFT_H


#include "common.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define  real float

//#if real == float
//#define  real_sin sinf
//#define  real_cos cosf
//#elif real == double
//#define  real_sin sin
//#define  real_cos cos
//#endif

#define  real_sin sinf
#define  real_cos cosf


typedef struct pack_header {
    int N;
    int P[15];
    int P_n;
    real *W;
    real **time_domain;
    real **freq_domain;
} PACK_HEADER;



PACK_HEADER *pack_new(real **in, real **out, int N);

void pack_fft(PACK_HEADER *header, fft_type type);

void pack_free(PACK_HEADER *header);

void fft_pack_real_test();

#endif //ALGORITHM_REAL_FFT_H
