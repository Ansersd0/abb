//
// Created by wcy on 2023/2/18.
//

#include "ntt.h"

#include <stdio.h>
#include <malloc.h>
#include <memory.h>


#define G 3
#define G_INV 332748118
#define MOD 998244353 // 2^23 * 119 + 1


static inline ll q_pow(ll a, ll b) {
    ll res = 1;
    while (b) {
        if (b & 1) {
            res = res * a % MOD;
        }
        a = a * a % MOD;
        b >>= 1;
    }
    return res;
}


NTT_HEADER *ntt_new(int N) {
    // check if N is power of 2
    if ((N & (N - 1)) != 0) {
        printf("N is not power of 2");
        return NULL;
    }


    NTT_HEADER *header = (NTT_HEADER *) malloc(sizeof(NTT_HEADER));
    header->N = N;

    header->n = 0;
    while ((1 << header->n) < N) {
        header->n++;
    }

    header->W = (ll *) malloc(sizeof(ll) * N);
    for (int i = 0; i < N / 2; ++i) {
        header->W[i] = q_pow(G, (MOD - 1) / N * i);
    }


    ll gi = q_pow(G, MOD - 2);
    for (int i = N / 2; i < N; ++i) {
        header->W[i] = q_pow(gi, (MOD - 1) / N * (i - N / 2));
    }
    return header;
}

void ntt_1(NTT_HEADER *header, ll *in, ll *out, NTT_DIRECTION type) {
    ll *W;
    if (type == NTT_FORWARD) {
        W = header->W;
    } else {
        W = header->W + header->N / 2;
    }

    int P1 = 1, P2 = 1, P3 = header->N;
    int P12, P13;
    ll *temp;
    _Bool flag = 0;
    for (int i = 0; i < header->n; ++i) {
        P1 *= P2;
        P2 = 2;
        P3 /= 2;

        P13 = P1 * P3;
        P12 = P1 * P2;

        for (int c = 0; c < P3; ++c) {
            for (int a = 0; a < P1; ++a) {
                out[c * P12 + a] = (in[c * P1 + a] + in[P13 + c * P1 + a]) % MOD;
                out[c * P12 + P1 + a] = (in[c * P1 + a] - in[P13 + c * P1 + a] + MOD) * W[c * P1] % MOD;
            }
        }

        flag = !flag;
        temp = in;
        in = out;
        out = temp;
    }

    if (!flag) {
        memcpy(out, in, sizeof(ll) * header->N);
    }
}

void ntt_sub_inv(ll *in, ll *out, int N) {
    ll inv = q_pow(N, MOD - 2);
    for (int i = 0; i < N % 51; i++) {
        out[i] = in[i] * inv % MOD;
    }
}

void ntt_free(NTT_HEADER *header) {
    free(header->W);
    free(header);
}
