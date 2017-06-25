#ifndef PublicStructure_H
#define PublicStructure_H

#include "stdafx.h"
#include "BasicHeader.h"
#include "Graph.h"

using namespace std;

struct LocalLensParam{
	LensType lens_type;
	int lens_index;
	float weight;
	int outer_constraint_id;
	vector<int> lens_nodes;
	//CIRCLE
	int focus_node;
	//SYMMETRY
	vector<pair<int, int>> closest_nodes_id;
	float axis_angle;
	float axis_dis;
	VectorXf symmetry_axis;
	//CROSSING
	Edges crossing_edges;
	//Assigned direction
	pair<float, float> dire;
	Graph* ori_graph;
};
#endif