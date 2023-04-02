//
// Created by wcy on 3/28/23.
//

#include <malloc.h>
#include "linkList.h"

LinkList *LinkListCreate() {
    LinkList *linkList = malloc(sizeof(LinkList));
    linkList->head = NULL;
    linkList->tail = NULL;
    linkList->size = 0;
    return linkList;
}

int LinkListAppend(LinkList *linkList, int value) {
    LinkListNode *node = malloc(sizeof(LinkListNode));
    node->value = value;
    node->next = NULL;
    if (linkList->size == 0) {
        linkList->head = node;
        linkList->tail = node;
    } else {
        linkList->tail->next = node;
        linkList->tail = node;
    }
    linkList->size++;
    return 1;
}

int LinkListInsert(LinkList *linkList, int index, int value) {
    if (index < 0 || index > linkList->size) {
        return 0;
    }
    LinkListNode *node = malloc(sizeof(LinkListNode));
    node->value = value;
    if (index == 0) {
        node->next = linkList->head;
        linkList->head = node;
    } else if (index == linkList->size) {
        node->next = NULL;
        linkList->tail->next = node;
        linkList->tail = node;
    } else {
        LinkListNode *prev = linkList->head;
        for (int i = 0; i < index - 1; ++i) {
            prev = prev->next;
        }
        node->next = prev->next;
        prev->next = node;
    }
    linkList->size++;
    return 1;
}

int LinkListDelete(LinkList *linkList, int index) {
    if (index < 0 || index >= linkList->size) {
        return 0;
    }
    LinkListNode *node;
    if (index == 0) {
        node = linkList->head;
        linkList->head = node->next;
        if (linkList->size == 1) {
            linkList->tail = NULL;
        }
    } else {
        LinkListNode *prev = linkList->head;
        for (int i = 0; i < index - 1; ++i) {
            prev = prev->next;
        }
        node = prev->next;
        prev->next = node->next;
        if (index == linkList->size - 1) {
            linkList->tail = prev;
        }
    }

    linkList->size--;
    free(node);
    return 1;
}

int LinkListGet(LinkList *linkList, int index) {
    if (index < 0 || index >= linkList->size) {
        return 0;
    }
    LinkListNode *node = linkList->head;
    for (int i = 0; i < index; ++i) {
        node = node->next;
    }
    return node->value;
}

int LinkListSize(LinkList *linkList) {
    return linkList->size;
}

int LinkListEmpty(LinkList *linkList) {
    LinkListNode *node = linkList->head;
    while (node != NULL) {
        LinkListNode *next = node->next;
        free(node);
        node = next;
    }
    return linkList->size == 0;
}

void LinkListDestroy(LinkList *linkList) {
    LinkListNode *node = linkList->head;
    while (node != NULL) {
        LinkListNode *next = node->next;
        free(node);
        node = next;
    }
    free(linkList);
}
