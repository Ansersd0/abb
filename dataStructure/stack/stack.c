//
// Created by wcy on 3/28/23.
//

#include <malloc.h>
#include "stack.h"

Stack *StackCreate(int size) {
    Stack *stack = malloc(sizeof(Stack));
    stack->data = malloc(sizeof(int) * size);
    stack->top = -1;
    stack->size = size;
    return stack;
}

int StackPush(Stack *stack, int value) {
    if (stack->top == stack->size - 1) {
        return -1;
    }
    stack->data[++stack->top] = value;
}

int StackPop(Stack *stack) {
    if (stack->top == -1) {
        return -1;
    }
    return stack->data[stack->top--];
}

int StackTop(Stack *stack) {
    if (stack->top == -1) {
        return -1;
    }
    return stack->data[stack->top];
}

int StackSize(Stack *stack) {
    return stack->top + 1;
}

int StackEmpty(Stack *stack) {
    return stack->top == -1;
}

void StackDestroy(Stack *stack) {
    free(stack->data);
    free(stack);
}
