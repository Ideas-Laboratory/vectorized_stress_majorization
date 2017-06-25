
#ifndef AtomConst_H
#define AtomConst_H

#include "Constraint.h"
#include "BasicHeader.h" // to make life easier
#include "cuSolver.cuh"
#include "Tool.h"

/************************************************************
** Atomconstraint is used for establishing other complex constraints.
** Atomconstraint is similar to stress constraint
** Except that it may have particular direction and length
**************************************************************/

class Solver;

class AtomConst : public Constraint
{
public:
	AtomConst();
	virtual ~AtomConst();

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

	static const int type = Stubs::ATOMCONSTRAINTS;
	int const_type = 0;


private:
	VectorXf d_vector;
	AtomConst(const AtomConst&);
	void operator = (const AtomConst&);
};

#endif