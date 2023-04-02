//
// Created by wcy on 3/28/23.
//

#ifndef ALGORITHM_ARRAY_H
#define ALGORITHM_ARRAY_H

typedef struct {
    int *data;
    int size;
} Array;

Array *ArrayCreate(int size);

int ArrayPush(Array *array, int value);

int ArrayPop(Array *array);

#endif //ALGORITHM_ARRAY_H

