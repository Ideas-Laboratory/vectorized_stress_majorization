#ifndef NonOverlapConst_H
#define NonOverlapConst_H

#include "Constraint.h"
#include "BasicHeader.h" 
#include "cuSolver.cuh"
#include "Tool.h"

class Solver;

/************************************************************************
NonOverlapConst is used for separating the overlapping rectangle nodes
************************************************************************/

class NonOverlapConst : public Constraint
{
public:
	NonOverlapConst();
	virtual ~NonOverlapConst();

	//const_id is the id of this constraint. extra_constraint begin with 0.
	void initShape(Edge edge, float r_length, float width, float height, float wi, int const_id);

	virtual void init();
	virtual void getRightHand(VectorXf& right_hand);
	virtual void setSolver(Solver* solver);
	virtual vector<int> getpid();
	virtual void setpid(vector<int> pid);
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
	virtual Edges getEdges();
	void setLMatrix(MatrixXf &sysm);
	virtual void minusLMatrix(MatrixXf &sysm);

	Edge edge;
	Solver* solver;
	vector<int> pid;
	float r_length;
	float wi;
	//pair<float, float> dire_vec = make_pair(0, 0);
	float node_width, node_height;

	//const_id is the id of this constraint. 
	//extra_constraint or temple_new_constriant begin from 0 
	int const_id;

	//the outer id of extra_constraint. A angle constraint includes several stress constraints, the outer_extra_const_id is corresponding to the the whole angle constraint
	//basic stress constraint id is -1
	int outer_extra_const_id = -1;

	static const int type = Stubs::NOVERLAP;

private:

	VectorXf d_vector;
	NonOverlapConst(const NonOverlapConst&);
	void operator = (const NonOverlapConst&);
};

#endif