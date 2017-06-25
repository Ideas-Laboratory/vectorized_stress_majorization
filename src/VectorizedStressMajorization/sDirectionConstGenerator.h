
#ifndef sDirectionConstGenerator_H
#define sDirectionConstGenerator_H

#include "Tool.h"
#include "BasicHeader.h" 
#include "AtomConst.h"

class Solver;

/***************************************************************************************************************************
sDirectionConstGenerator is used for transforming direction constraints into AtomConsts
Each sDirectionConstGenerator actions on an pair of node (u,v)
Given idea direction of (u,v)
***************************************************************************************************************************/
class sDirectionConstGenerator
{
public:
	sDirectionConstGenerator();
	virtual ~sDirectionConstGenerator();
	virtual void init();
	void initShape(Edge edge, float edge_length, pair<float, float> dire_vector, float weightScale, int const_id, int outer_extra_const_id);
	float getWi();
	virtual void setSolver(Solver* solver);

	int length_type_num; //floor(n/2) types in a circle. first = 1   
	float wi;
	VectorXf vertex_list;
	vector<int> vertex_ids;
	float edge_length;
	float weightScale;
	Solver* solver;
	Edge edge;
	pair<float, float> dire_vector;
	AtomConst* generateStressConst(bool isbigger);
	int const_id;
	int outer_extra_const_id;
};

#endif