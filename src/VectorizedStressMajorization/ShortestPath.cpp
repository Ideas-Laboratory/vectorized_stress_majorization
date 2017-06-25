#include "stdafx.h"
#include "ShortestPath.h"
using namespace std;

ShortestPath::ShortestPath(int max)
{
	this->maxint = max;
}

ShortestPath::~ShortestPath()
{
}



// n -- n nodes
// v -- the source node
// dist[] -- the distance from the ith node to the source node
// prev[] -- the previous node of the ith node
// c[][] -- every two nodes' distance
void ShortestPath::Dijkstra(int n, int v, int *dist, int *prev, MatrixXf c)
{
	int *s = new int[n];    // 判断是否已存入该点到S集合中
	for (int i = 0; i < n; i++)
	{
		dist[i] = c(v, i);
		s[i] = 0;     // 初始都未用过该点
		if (dist[i] == maxint)
			prev[i] = 0;
		else
			prev[i] = v;
	}
	dist[v] = 0;
	s[v] = 1;

	// 依次将未放入S集合的结点中，取dist[]最小值的结点，放入结合S中
	// 一旦S包含了所有V中顶点，dist就记录了从源点到所有其他顶点之间的最短路径长度
	// 注意是从第二个节点开始，第一个为源点
	for (int i = 1; i < n; i++)
	{
		int tmp = maxint;
		int u = v;
		// 找出当前未使用的点j的dist[j]最小值
		for (int j = 0; j < n; j++)
		if ((!s[j]) && dist[j]<tmp)
		{
			u = j;              // u保存当前邻接点中距离最小的点的号码
			tmp = dist[j];
		}
		s[u] = 1;    // 表示u点已存入S集合中

		// 更新dist
		for (int j = 0; j < n; j++)
		if ((!s[j]) && c(u, j)<maxint)
		{
			int newdist = dist[u] + c(u, j);
			if (newdist < dist[j])
			{
				dist[j] = newdist;
				prev[j] = u;
			}
		}
	}
}


//undirected, weight = 1 if the two vertexes is neighbour
MatrixXf ShortestPath::edgelistToShortestPath(int n, Edges& edgelist) {
	 
	int *dist = new int[n];     // 表示当前点到源点的最短路径长度 
	MatrixXf c(n, n);			// 记录图的两点间路径长度 邻接矩阵  
	int *prev = new int[n];     // 记录当前点的前一个结点
	 

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			c(i, j) = maxint; 
		}
	}
	//building the adjacency matrix c
	int sourceid = -1;
	int targetid = -1;
	for (int i = 0; i < edgelist.size(); i++)
	{
		sourceid = edgelist[i].first;
		targetid = edgelist[i].second;
		c(sourceid, targetid) = 1;
		c(targetid, sourceid) = 1; 
	}

	initToMaxint(maxint, n, dist);

	for (int i = 0; i < n; i++) {
		Dijkstra(n, i, dist, prev, c);
		for (int j = 0; j < n; j++) {
			c(i, j) = dist[j];
		}
		initToMaxint(maxint, n, dist);
	}

	return c;
}

void ShortestPath::initToMaxint(int maxint, int n, int*& dist) {
	for (int i = 0; i < n; i++)
	{
		dist[i] = maxint;

	}
}