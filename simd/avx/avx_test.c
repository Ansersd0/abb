//
// Created by wcy on 2023/3/20.
//

#include "avx.h"


int main() {
    complex_f8 *a = avx_malloc(sizeof(complex_f8) * 8);
    complex_f8 *b = avx_malloc(sizeof(complex_f8) * 8);
    complex_f8 *c = avx_malloc(sizeof(complex_f8) * 8);
    for (int i = 0; i < 8; ++i) {
        a[i] = C8_SET_16F1(1, 0,
                           2, 0,
                           3, 0,
                           4, 0,
                           5, 0,
                           6, 0,
                           7, 0,
                           8, 0);
        b[i] = C8_SET_16F1(1, 0,
                           2, 0,
                           3, 0,
                           4, 0,
                           5, 0,
                           6, 0,
                           7, 0,
                           8, 0);
    }

    for (int i = 0; i < 8; ++i) {
        c[i] = C8_ADD_C8(a[i], b[i]);
        c[i] = C8_MUL_C8(a[i], b[i]);
        c[i] = C8_MUL_C8(a[i], b[i]);
        c[i] = C8_MUL_I(a[i]);
        c[i] = C8_ADD_C8(C8_ADD_C8(C8_SUB_C8(a[i], b[i]), C8_MUL_C8(a[i], b[i])), C8_MUL_I(a[i]));
    }
}