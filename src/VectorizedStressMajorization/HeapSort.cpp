
#include "HeapSort.h"
using namespace std;

void HeapAdjust(std::vector<float> &a, std::vector<int> &ids, int i, int size)  //调整堆 
{
	int lchild = 2 * i;       //i的左孩子节点序号 
	int rchild = 2 * i + 1;     //i的右孩子节点序号 
	int max = i;            //临时变量 
	//std::cout << "heap adjust " << a.size() << "," << max << "," << lchild << "," << rchild << std::endl;
	if (i <= size / 2)          //如果i是叶节点就不用进行调整 
	{
		if (lchild <= size&&a[lchild - 1]>a[max - 1])
		{
			max = lchild;
		}
		if (rchild <= size&&a[rchild - 1]>a[max - 1])
		{
			max = rchild;
		}
		if (max != i)
		{
			swap(a[i - 1], a[max - 1]);
			swap(ids[i - 1], ids[max - 1]);
			HeapAdjust(a, ids, max, size);    //避免调整之后以max为父节点的子树不是堆 
		}
	}
}

void BuildHeap(std::vector<float> &a, std::vector<int> &ids, int size)    //建立堆 
{
	int i;
	for (i = size / 2; i >= 1; i--)    //非叶节点最大序号值为size/2 
	{
		HeapAdjust(a, ids, i, size);
	}
}


void HeapSort(std::vector<float> &a, int size, std::vector<int> &child_list)    //堆排序 
{
	int i;
	BuildHeap(a, child_list, size);
	for (i = size; i >= 1; i--)
	{
		swap(a[0], a[i - 1]);           //交换堆顶和最后一个元素，即每次将剩余元素中的最大者放到最后面 
		swap(child_list[0], child_list[i - 1]);
		//BuildHeap(a,i-1);        //将余下元素重新建立为大顶堆 
		HeapAdjust(a, child_list, 1, i - 1);      //重新调整堆顶节点成为大顶堆 
	}
}
