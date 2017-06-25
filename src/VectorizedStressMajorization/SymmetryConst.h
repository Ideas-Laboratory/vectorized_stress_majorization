
#ifndef SymmetryConst_H
#define SymmetryConst_H

#include "Constraint.h"
#include "BasicHeader.h" // to make life easier
#include "cuSolver.cuh"
#include "Tool.h"

class Solver;

/***************************************************************************************************************************
SymmetryConst is used for making a selected area symmetry
SymmetryConst actions on a set of nodes, given the symmetry-axis
***************************************************************************************************************************/

class SymmetryConst : public Constraint
{
public:
	SymmetryConst();
	virtual ~SymmetryConst();

	//const_id is the id of this constraint. extra_constraint begin with 0.
	void initShape(vector<pair<int, int>> closest_node_id_pairs, float* shortest_path, int ori_graph_n, float rotate_angle, float move_dis, float wi, int const_id);

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
	virtual Edges getEdges();
	virtual void setWi(float wi);
	void setLMatrix(MatrixXf &sysm);
	virtual void minusLMatrix(MatrixXf &sysm);

	Edges edges, edges_model;
	Solver* solver;
	vector<int> pid;
	vector<float> r_lengths;
	float wi;
	vector<pair<int, int>> closest_node_id_pairs;
	float rotate_angle;
	float move_dis;

	//const_id is the id of this constraint. 
	//extra_constraint or temple_new_constriant begin from 0 
	int const_id;

	//the outer id of extra_constraint. A angle constraint includes several stress constraints, the outer_extra_const_id is corresponding to the the whole angle constraint
	//basic stress constraint id is -1
	int outer_extra_const_id = -1;

	static const int type = Stubs::SYMMETRYCONSTRAINTS;


private:

	VectorXf d_vector;
	SymmetryConst(const SymmetryConst&);
	void operator = (const SymmetryConst&);
};

#endif