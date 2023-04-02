//
// Created by wcy on 3/28/23.
//

#include <stdio.h>
#include "stack.h"

int main() {
    int N = 10;
    Stack *stack = StackCreate(N);

    for (int i = 0; i < N; ++i) {
        StackPush(stack, i);
    }

    while (!StackEmpty(stack)) {
        printf("%d ", StackPop(stack));
    }
    printf("\n");
}
