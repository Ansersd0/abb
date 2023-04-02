//
// Created by wcy on 2023/3/17.
//

#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pft.h"
#include  "../fft/fft1.h"

int main() {
    int N = 50;
    int change = 1;

    double _Complex *in = (double _Complex *) malloc(sizeof(double _Complex) * (N + abs(change)));
    double _Complex *out = (double _Complex *) malloc(sizeof(double _Complex) * (N + abs(change)));
    for (int i = 0; i < N; i++) {
        in[i] = i + 1;
    }
//    PFT_HEADER *p = pft_new(&in, &out, N, change);
//    pft(p);
//    for (int i = 0; i < N + change; i++) {
//        printf("%f\n", in[i]);
//    }

    my_fft(in, out, N, Forward);
    if (change > 0) {
        if (N % 2 == 0) {
            memcpy(out + N / 2 + change,
                   out + N / 2,
                   sizeof(double _Complex) * (N / 2));
            memset(out + N / 2, 0, sizeof(double _Complex) * (change - 1));
        } else {
            memcpy(out + (N + 1) / 2 + 1 + change,
                   out + (N + 1) / 2 + 1,
                   sizeof(double _Complex) * ((N - 1) / 2));
            memset(out + N / 2, 0, sizeof(double _Complex) * change);
        }
    }
    my_fft(out, in, N + change, Inverse);


    for (int i = 0; i < N + change; i++) {
//        printf("%f\n", sqrt(creal(in[i]) * creal(in[i]) + cimag(in[i]) * cimag(in[i])) / N);
        printf("%f + %fi\n", creal(in[i]) / N, cimag(in[i]) / N);
    }

    return 0;
}