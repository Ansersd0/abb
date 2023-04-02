//
// Created by wcy on 2023/3/7.
//

#include <malloc.h>
#include <stdio.h>
#include <time.h>
#include "ntt.h"

int main() {
    int N = 1 << 22;
    ll *in = (ll *) malloc(sizeof(ll) * N);
    ll *out = (ll *) malloc(sizeof(ll) * N);

    for (int i = 0; i < N; i++) {
        in[i] = (i + 1) % 100;
    }

    NTT_HEADER *ntt_header = ntt_new(N);

    // get now time
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);

    ntt_1(ntt_header, in, out, NTT_FORWARD);

    ntt_1(ntt_header, out, in, NTT_INVERSE);

    // get now time
    clock_gettime(CLOCK_MONOTONIC, &finish);

    ntt_sub_inv(in, out, N);

    for (int i = 0; i < 10; i++) {
        printf("%lld ", out[i]);
    }
    printf("\n");

    printf("time: %lf ms\n",
           (double) (finish.tv_sec - start.tv_sec) * 1000.0 + (finish.tv_nsec - start.tv_nsec) / 1000000.0);

    free(in);
    free(out);
    ntt_free(ntt_header);

    return 0;
}

