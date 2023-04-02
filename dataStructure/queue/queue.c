//
// Created by wcy on 3/28/23.
//

#include <malloc.h>
#include "queue.h"

Queue *QueueCreate(int size) {
    Queue *queue = malloc(sizeof(Queue));
    queue->data = malloc(sizeof(int) * size);
    queue->head = 0;
    queue->tail = 0;
    queue->size = size;
    return queue;
}

int QueuePush(Queue *queue, int value) {
    if (queue->tail == queue->size) {
        if (queue->head == 0) {
            return 0;
        }
        for (int i = queue->head; i < queue->tail; ++i) {
            queue->data[i - queue->head] = queue->data[i];
        }
        queue->tail -= queue->head;
        queue->head = 0;
    }
    queue->data[queue->tail++] = value;
    return 1;
}

int QueuePop(Queue *queue) {
    if (queue->head == queue->tail) {
        return 0;
    }
    return queue->data[queue->head++];
}

int QueueFront(Queue *queue) {
    if (queue->head == queue->tail) {
        return 0;
    }
    return queue->data[queue->head];
}

int QueueSize(Queue *queue) {
    return queue->tail - queue->head;
}

int QueueEmpty(Queue *queue) {
    return queue->head == queue->tail;
}

void QueueDestroy(Queue *queue) {
    free(queue->data);
    free(queue);
}
