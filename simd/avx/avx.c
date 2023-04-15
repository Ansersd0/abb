//
// Created by wcy on 2023/3/15.
//

#include <math.h>
#include "avx.h"

void *avx_malloc(size_t size) {
//    void *ptr = malloc(size + 32);
//    return (void *) (((size_t) ptr + 32) & ~0x1f);
//    return _mm_malloc(size, 32);
    return aligned_alloc(32, size);
}

void avx_free(void *ptr) {
//    free((void *) (((size_t) ptr) & ~0x1f));
//    _mm_free(ptr);
    free(ptr);
}
