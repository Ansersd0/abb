//
// Created by wcy on 3/28/23.
//
#include <stdio.h>
#include "linkList.h"

int main() {
    LinkList *linkList = LinkListCreate();
    LinkListAppend(linkList, 1);
    LinkListAppend(linkList, 2);
    LinkListAppend(linkList, 3);
    LinkListAppend(linkList, 4);

    for (int i = 0; i < LinkListSize(linkList); ++i) {
        printf("%d ", LinkListGet(linkList, i));
    }
    printf("\n");

    LinkListInsert(linkList, 0, 0);
    LinkListInsert(linkList, 2, 5);

    for (int i = 0; i < LinkListSize(linkList); ++i) {
        printf("%d ", LinkListGet(linkList, i));
    }
    printf("\n");

    LinkListDelete(linkList, 0);
    LinkListDelete(linkList, 2);

    for (int i = 0; i < LinkListSize(linkList); ++i) {
        printf("%d ", LinkListGet(linkList, i));
    }
    printf("\n");

    LinkListDestroy(linkList);
}