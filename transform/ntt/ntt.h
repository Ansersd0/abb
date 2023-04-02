//
// Created by wcy on 2023/2/18.
//

#ifndef ALGORITHM_NTT_H
#define ALGORITHM_NTT_H


typedef enum{
    NTT_FORWARD = -1,
    NTT_INVERSE = 1,
} NTT_DIRECTION;


#ifndef ll
typedef long long ll;
#endif

typedef struct {
    int N;
    int n;
    ll *W;
} NTT_HEADER;

NTT_HEADER *ntt_new(int N);

void ntt_1(NTT_HEADER *header, ll *in, ll *out, NTT_DIRECTION type);

void ntt_sub_inv(ll *in, ll *out, int N);

void ntt_free(NTT_HEADER *header);

#endif //ALGORITHM_NTT_H
