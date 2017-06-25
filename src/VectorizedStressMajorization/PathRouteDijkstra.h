#ifndef PathRouteDijkstra_H
#define PathRouteDijkstra_H

#include "BasicHeader.h"
#include <string> 
using namespace std;


class PathRouteDijkstra
{
public:
	PathRouteDijkstra();
	~PathRouteDijkstra();
	//int n, line;             // 图的结点数和路径数 
	void init(int n, int source, int target, MatrixXf adjMatrix);
	vector<int>  PathRoute(int n, int source, int target, MatrixXf adjMatrix);
	int* getAllBCValues(int n, MatrixXf adjMatrix);
	void getPathitoj(vector<int> &pathitoj, int i, int j);
	int getClosestUnmarkedNode();
	void printPath(int node);

	int* predecessor;
	int* distance;
	bool* mark;
	int N;
	int source;
	int target;
	MatrixXf adjMatrix;
	vector<int> path;
};
#endif