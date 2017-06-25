#include "stdafx.h"
#include "PathRouteDijkstra.h"
#define INFINITY 999

using namespace std;

PathRouteDijkstra::PathRouteDijkstra()
{
}

PathRouteDijkstra::~PathRouteDijkstra()
{
}

void PathRouteDijkstra::init(int n, int source, int target, MatrixXf adjMatrix){

	this->N = n;
	this->source = source;
	this->target = target;
	this->adjMatrix = adjMatrix;

	predecessor = new int[N];
	distance = new int[N];
	mark = new bool[N]; //keep track of visited node

	for (int i = 0; i<N; i++) {
		mark[i] = false;
		predecessor[i] = -1;
		distance[i] = INFINITY;
	}
	distance[source] = 0;

}

int PathRouteDijkstra::getClosestUnmarkedNode(){
	int minDistance = INFINITY;
	int closestUnmarkedNode;
	for (int i = 0; i<N; i++) {
		if ((!mark[i]) && (minDistance >= distance[i])) {
			minDistance = distance[i];
			closestUnmarkedNode = i;
		}
	}
	return closestUnmarkedNode;
}


vector<int> PathRouteDijkstra::PathRoute(int n, int source, int target, MatrixXf adjMatrix){
	init(n, source, target, adjMatrix);

	int closestUnmarkedNode;
	int count = 0;
	while (count < N) {
		closestUnmarkedNode = getClosestUnmarkedNode();
		mark[closestUnmarkedNode] = true;
		for (int i = 0; i < N; i++) {
			if ((!mark[i]) && (adjMatrix(closestUnmarkedNode, i)>0)) {
				if (distance[i] > distance[closestUnmarkedNode] + adjMatrix(closestUnmarkedNode, i)) {
					distance[i] = distance[closestUnmarkedNode] + adjMatrix(closestUnmarkedNode, i);
					predecessor[i] = closestUnmarkedNode;
				}
			}
		}
		count++;
	}
	printPath(target);
	return path;
}


void PathRouteDijkstra::printPath(int node){
	if (node == source)
		//cout << target << "..";
		path.push_back(source);
	else if (predecessor[node] == -1)
		//cout << "No path from "<< source <<" to " << node << endl;
		path.push_back(-2);
	else {
		printPath(predecessor[node]);
		//cout << node << "..";
		path.push_back(node);
	}
}

int* PathRouteDijkstra::getAllBCValues(int n, MatrixXf adjMatrix){

	this->N = n;
	this->adjMatrix = adjMatrix;

	vector<vector<int>> paths;
	int* BC_all = new int[N];
	for (int i = 0; i < N; i++){
		BC_all[i] = 0;
	}

	predecessor = new int[N];
	distance = new int[N];
	mark = new bool[N]; //keep track of visited node
	vector<int> pathitoj;
	int count = 0;

	for (int i = 0; i < this->N; i++){

		// for every source, the predecessor and distance are the same, calculate tnem once
		for (int k = 0; k<N; k++) {
			mark[k] = false;
			predecessor[k] = -1;
			distance[k] = INFINITY;
		}
		distance[i] = 0;

		int closestUnmarkedNode;
		count = 0;
		while (count < N) {
			closestUnmarkedNode = getClosestUnmarkedNode();
			mark[closestUnmarkedNode] = true;
			for (int t = 0; t < N; t++) {
				if ((!mark[t]) && (adjMatrix(closestUnmarkedNode, t)>0)) {
					if (distance[t] > distance[closestUnmarkedNode] + adjMatrix(closestUnmarkedNode, t)) {
						distance[t] = distance[closestUnmarkedNode] + adjMatrix(closestUnmarkedNode, t);
						predecessor[t] = closestUnmarkedNode;
					}
				}
			}
			count++;
		}

		//for every i to j, corresponding to a path
		for (int j = i + 1; j < this->N; j++){
			pathitoj.clear();
			getPathitoj(pathitoj, i, j);
			if (pathitoj.size() != 0)
				paths.push_back(pathitoj);
			for (int k = 1; k < pathitoj.size() - 1; k++){
				BC_all[pathitoj[k]]++;
			}
		}
	}

	return BC_all;
}

void PathRouteDijkstra::getPathitoj(vector<int> &pathitoj, int i, int j){
	if (j == i)
		pathitoj.push_back(i);
	else if (predecessor[j] == -1)
		pathitoj.push_back(-2);
	else {
		getPathitoj(pathitoj, i, predecessor[j]);
		pathitoj.push_back(j);
	}
}