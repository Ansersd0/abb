//
// Created by wcy on 3/28/23.
//

#ifndef ALGORITHM_LINKLIST_H
#define ALGORITHM_LINKLIST_H

typedef struct LinkListNode {
    int value;
    struct LinkListNode *next;
} LinkListNode;

typedef struct {
    LinkListNode *head;
    LinkListNode *tail;
    int size;
} LinkList;

LinkList *LinkListCreate();

int LinkListAppend(LinkList *linkList, int value);

int LinkListInsert(LinkList *linkList, int index, int value);

int LinkListDelete(LinkList *linkList, int index);

int LinkListGet(LinkList *linkList, int index);

int LinkListSize(LinkList *linkList);

int LinkListEmpty(LinkList *linkList);

void LinkListDestroy(LinkList *linkList);

#endif //ALGORITHM_LINKLIST_H
