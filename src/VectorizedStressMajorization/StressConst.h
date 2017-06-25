
#ifndef StressConst_H
#define StressConst_H

#include "Constraint.h"
#include "BasicHeader.h" // to make life easier
#include "cuSolver.cuh"
#include "Tool.h"

class Solver;

/***************************************************************************************************************************
StressConst is basic for the initial graph layout with stress.
One StressConst is action on a pair of nodes (u,v). 
There will be |V|*|V-1|/2 pair of nodes in the graph G.
The idea length is the shortest path between u and v.
The default direction is obtained from last iteration.
In the subsequent proccess will combine these StresConsts with other constriants
***************************************************************************************************************************/

class StressConst : public Constraint
{
public:
	StressConst();
	virtual ~StressConst();

	//const_id is the id of this constraint. extra_constraint begin with 0.
	void initShape(Edge edge, float r_length, pair<float, float> dire_vec, float wi, int const_id);

	virtual void init();
	virtual void getRightHand(VectorXf& right_hand);
	virtual void setSolver(Solver* solver);
	virtual vector<int> getpid();
	int getType() {
		return type;
	}
	int Copy(void *linear_constraint);
	float getWi();
	virtual int getdnum();
	virtual int  getConstId();
	virtual void  setConstId(int const_id);
	virtual void setOuterExtraConstId(int outer_extra_const_id);
	virtual int getOuterExtraConstId();
	virtual void setWi(float wi);
	virtual void setLMatrix(MatrixXf &sysm);
	void setLMatrix(MatrixXf &sysm, float **distMatrix);
	virtual void minusLMatrix(MatrixXf &sysm);

	Edge edge;
	Solver* solver;
	vector<int> pid;
	float r_length;
	float wi;
	pair<float, float> dire_vec = make_pair(0, 0);
	virtual Edges getEdges();

	//const_id is the id of this constraint. 
	//extra_constraint or temple_new_constriant begin from 0 
	int const_id;

	//the outer id of extra_constraint. A angle constraint includes several stress constraints, the outer_extra_const_id is corresponding to the the whole angle constraint
	//basic stress constraint id is -1
	int outer_extra_const_id = -1;

	static const int type = Stubs::STRESSCONSTRAINTS;
	int const_type = 0;


private:
	VectorXf d_vector;
	StressConst(const StressConst&);
	void operator = (const StressConst&);
};

#endif