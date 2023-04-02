//
// Created by wcy on 2023/3/17.
//

#ifndef ALGORITHM_PFT_H
#define ALGORITHM_PFT_H

#include "../fft/real_fft.h"

typedef struct {
    int N;
    int change;

    struct pack_header *fft;
    struct pack_header *i_fft;
} PFT_HEADER;

PFT_HEADER *pft_new(float **in, float **out, int N, int change);

void pft_free(PFT_HEADER *pft);

void pft(PFT_HEADER *pft);

void pft_change_fs(PFT_HEADER *pft, int fs);

void pft_change_timbres(PFT_HEADER *pft, int semitone, int cents);

#endif //ALGORITHM_PFT_H
