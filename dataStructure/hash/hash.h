//
// Created by wcy on 3/28/23.
//

#ifndef ALGORITHM_HASH_H
#define ALGORITHM_HASH_H

typedef struct hash_table {
    int size;
    int count;
    int *table;
} hash_table;

hash_table *hashCreate(int size);

bool hashInsert(hash_table *hash, int key);

bool hashSearch(hash_table *hash, int key);

bool hashDelete(hash_table *hash, int key);

bool hashDestroy(hash_table *hash);

#endif //ALGORITHM_HASH_H
