
#ifndef CircleConst_H
#define CircleConst_H

#include "Constraint.h"
#include "BasicHeader.h" // to make life easier
#include "Tool.h"

class Solver;
/************************************************************************************
CircleConst is used for making a ring(or defined node sequence) more circular
**************************************************************************************/
class CircleConst : public Constraint
{
public:
	CircleConst();
	virtual ~CircleConst();
	void initShape(vector<int> pid, vector<float> r_lengths, float wi, int const_id);
	virtual void init();
	virtual void getRightHand(VectorXf& right_hand);
	virtual void setSolver(Solver* solver);
	virtual vector<int> getpid();
	virtual void setpid(vector<int> pid);
	float getWi();
	virtual int getdnum();
	virtual int getConstId();
	virtual void setConstId(int const_id);
	virtual void setOuterExtraConstId(int outer_extra_const_id);
	virtual int getOuterExtraConstId();
	virtual void setWi(float wi);
	virtual Edges getEdges();
	vector<float> getLocalRelativeShape();
	void setLMatrix(MatrixXf &sysm);
	virtual void minusLMatrix(MatrixXf &sysm);

	int getType() {
		return type;
	}
	int Copy(void *linear_constraint);

	Edges edges;
	Solver* solver;
	vector<int> pid;//ordered anticlockwise
	vector<float> r_lengths;
	float wi;
	float differ_angle;
	Edge anchor_edge;
	pair<float, float> dire_vec = make_pair(0, 0);
	vector<float> relative_shape;

	//const_id is the id of this constraint. 
	//extra_constraint or temple_new_constriant begin from 0 
	int const_id;

	//the outer id of extra_constraint. A angle constraint includes several stress constraints, the outer_extra_const_id is corresponding to the the whole angle constraint
	//basic stress constraint id is -1
	int outer_extra_const_id = -1;

	static const int type = Stubs::CIRCLE;
private:
	VectorXf d_vectors;

private:
	CircleConst(const CircleConst&);
	void operator = (const CircleConst&);
};

#endif