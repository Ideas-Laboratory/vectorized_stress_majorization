#ifndef ShortestPath_H
#define ShortestPath_H

#include "BasicHeader.h"
#include <string> 
using namespace std;


const int maxnum = 100000;
//const int maxint = 99;
class ShortestPath
{
public:
	ShortestPath(int max);
	~ShortestPath();
	void Dijkstra(int n, int v, int *dist, int *prev, MatrixXf c);
	MatrixXf edgelistToShortestPath(int n, Edges& edgelist);
	void initToMaxint(int maxint, int n, int*& dist);
	int maxint;
	//int n, line;             // 图的结点数和路径数
};
#endif