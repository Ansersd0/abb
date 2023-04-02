//
// Created by wcy on 3/28/23.
//

#include <malloc.h>
#include "hash.h"

hash_table *hashCreate(int size) {
    hash_table *hash = (hash_table *) malloc(sizeof(hash_table));
    hash->size = size;
    hash->count = 0;
    hash->table = (int *) malloc(sizeof(int) * size);
    for (int i = 0; i < size; i++) {
        hash->table[i] = 0;
    }
    return hash;
}

bool hashInsert(hash_table *hash, int key) {
    int index = key % hash->size;
    while (hash->table[index] != 0) {
        index = (index + 1) % hash->size;
    }
    hash->table[index] = key;
    hash->count++;
    return true;
}

bool hashSearch(hash_table *hash, int key) {
    int index = key % hash->size;
    while (hash->table[index] != key) {
        index = (index + 1) % hash->size;
        if (hash->table[index] == 0 || index == key % hash->size) {
            return false;
        }
    }
    return true;
}

bool hashDelete(hash_table *hash, int key) {
    int index = key % hash->size;
    while (hash->table[index] != key) {
        index = (index + 1) % hash->size;
        if (hash->table[index] == 0 || index == key % hash->size) {
            return false;
        }
    }
    hash->table[index] = 0;
    hash->count--;
    return true;
}

bool hashDestroy(hash_table *hash) {
    free(hash->table);
    free(hash);
    return true;
}
