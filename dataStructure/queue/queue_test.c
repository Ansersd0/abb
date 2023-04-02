//
// Created by wcy on 3/28/23.
//
#include <stdio.h>
#include "queue.h"

int main() {
    int N = 10;
    Queue *queue = QueueCreate(N);

    for (int i = 0; i < N; ++i) {
        QueuePush(queue, i);
    }

    while (!QueueEmpty(queue)) {
        printf("%d ", QueuePop(queue));
    }
    printf("\n");
}