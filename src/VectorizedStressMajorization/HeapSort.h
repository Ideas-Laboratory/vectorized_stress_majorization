
#ifndef HeapSort_H
#define HeapSort_H

#include <iostream>
#include <algorithm>
#include "BasicHeader.h"

void HeapAdjust(std::vector<float> &a, std::vector<int> &ids, int i, int size);
void BuildHeap(std::vector<float> &a, std::vector<float> &ids, int size);
void HeapSort(std::vector<float> &a, int size, std::vector<int> &ids);

#endif