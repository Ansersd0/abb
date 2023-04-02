//
// Created by wcy on 3/28/23.
//
#include <stdio.h>
#include "hash.h"

int main() {
    hash_table *hash = hashCreate(10);
    hashInsert(hash, 1);
    hashInsert(hash, 2);
    hashInsert(hash, 3);
    hashInsert(hash, 4);
    hashInsert(hash, 5);

    printf("%d\n", hashSearch(hash, 1));
    printf("%d\n", hashSearch(hash, 2));
    printf("%d\n", hashSearch(hash, 3));

    hashDelete(hash, 1);
    hashDelete(hash, 2);

    printf("%d\n", hashSearch(hash, 1));
    printf("%d\n", hashSearch(hash, 2));

    hashDestroy(hash);
}