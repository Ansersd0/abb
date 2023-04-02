//
// Created by wcy on 3/28/23.
//

#ifndef ALGORITHM_STACK_H
#define ALGORITHM_STACK_H

typedef struct {
    int *data;
    int top;
    int size;
} Stack;

Stack *StackCreate(int size);

int StackPush(Stack *stack, int value);

int StackPop(Stack *stack);

int StackTop(Stack *stack);

int StackSize(Stack *stack);

int StackEmpty(Stack *stack);

void StackDestroy(Stack *stack);


#endif //ALGORITHM_STACK_H
