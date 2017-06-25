
#ifndef CrossingRemoveConst_H
#define CrossingRemoveConst_H

#include "Constraint.h"
#include "BasicHeader.h"
#include "Tool.h"


/************************************************************************
CrossingRemoveConst is used for removing the constraints on each pair of edges
*************************************************************************/
class Solver;

class CrossingRemoveConst : public Constraint
{
public:
	CrossingRemoveConst();
	virtual ~CrossingRemoveConst();
	void initShape(Edges edges, vector<float> r_lengths, float wi, vector<int> degrees, int const_id);

	virtual void init();
	virtual void getRightHand(VectorXf& right_hand);
	virtual void setSolver(Solver* solver);
	virtual vector<int> getpid();
	virtual void setpid(vector<int> pid);
	float getWi();
	virtual int getdnum();
	int Copy(void *constraints);
	int getType() {
		return type;
	}
	virtual int getConstId();
	virtual void  setConstId(int const_id);
	virtual void setOuterExtraConstId(int outer_extra_const_id);
	virtual int getOuterExtraConstId();
	virtual void setWi(float wi);
	virtual Edges getEdges();
	void setLMatrix(MatrixXf &sysm);
	virtual void minusLMatrix(MatrixXf &sysm);

	Edges edges;
	Solver* solver;
	MatrixXf N;
	vector<int> pid;
	vector<float> r_lengths;
	vector<float> edge_weights;
	float wi;

	//const_id is the id of this constraint. 
	//extra_constraint or temple_new_constriant begin from 0 
	int const_id;

	//the outer id of extra_constraint. A angle constraint includes several stress constraints, the outer_extra_const_id is corresponding to the the whole angle constraint
	//basic stress constraint id is -1
	int outer_extra_const_id = -1;

	static const int type = Stubs::CROSSINGREMOVAL;
private:

	VectorXf d_vector;
	CrossingRemoveConst(const CrossingRemoveConst&);
	void operator = (const CrossingRemoveConst&);
};

#endif