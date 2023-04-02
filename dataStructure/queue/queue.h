//
// Created by wcy on 3/28/23.
//

#ifndef ALGORITHM_QUEUE_H
#define ALGORITHM_QUEUE_H

typedef struct {
    int *data;
    int head;
    int tail;
    int size;
} Queue;

Queue *QueueCreate(int size);

int QueuePush(Queue *queue, int value);

int QueuePop(Queue *queue);

int QueueFront(Queue *queue);

int QueueSize(Queue *queue);

int QueueEmpty(Queue *queue);

void QueueDestroy(Queue *queue);

#endif //ALGORITHM_QUEUE_H
