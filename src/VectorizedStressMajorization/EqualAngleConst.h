
#ifndef EqualAngleConst_H
#define EqualAngleConst_H

#include "Constraint.h"
#include "BasicHeader.h" // to make life easier
#include "Tool.h"
#include "HeapSort.h"

class Solver;

/*************************************************************************************
EqualAngleConst is used for maximizing the minimum angle in a star-structure
In the other word, it makes the nodes leaving a center node with even angles
*******************************************************************************/

class EqualAngleConst : public Constraint
{
public:
	EqualAngleConst();
	virtual ~EqualAngleConst();
	void initShape(int parentid, vector<int> child_list, float plength, float wi, int const_id);

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
	float plength;
	float wi_coe;
	float wi;
	float differ_angle;
	Edge anchor_edge;
	pair<float, float> dire_vec = make_pair(0, 0);
	vector<float> relative_shape;
	int begin_id;
	int parentid;
	vector<int> child_list;

	//const_id is the id of this constraint. 
	//extra_constraint or temple_new_constriant begin from 0 
	int const_id;

	//the outer id of extra_constraint. A angle constraint includes several stress constraints, the outer_extra_const_id is corresponding to the the whole angle constraint
	//basic stress constraint id is -1
	int outer_extra_const_id = -1;

	static const int type = Stubs::EUQALANGLE;
private:
	VectorXf d_vectors;

private:
	EqualAngleConst(const EqualAngleConst&);
	void operator = (const EqualAngleConst&);
};

#endif