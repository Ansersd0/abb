//
// Created by wcy on 2023/3/17.
//
#include <stdlib.h>
#include <string.h>
#include "pft.h"


PFT_HEADER *pft_new(float **in, float **out, int N, int change) {
    PFT_HEADER *pft = malloc(sizeof(PFT_HEADER));
    pft->N = N;
    pft->change = change;
    pft->fft = pack_new(in, out, N);
    pft->i_fft = pack_new(in, out, N + change);
    return pft;
}

void pft(PFT_HEADER *pft) {
    pack_fft(pft->fft, Forward);

    if (pft->change > 0) {
        memset(*(pft->fft->freq_domain) + pft->N, 0, sizeof(float) * pft->change);
    }
    pack_fft(pft->i_fft, Inverse);

    for (int i = 0; i < pft->N + pft->change; i++) {
        *(*(pft->i_fft->time_domain) + i) /= pft->N;
    }

}
