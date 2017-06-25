#ifndef Solver_H
#define Solver_H


#include "BasicHeader.h"
#include "Tool.h"
#include "cuSolver.cuh"


class Constraint;

class Solver
{
public:
	Solver();
	~Solver();
	cuDWORD *linear_constraints = 0;
	void addConstraint(Constraint* constraint);
	void addExtraConstraint(Constraint* constraint);
	void clearConstraints();
	/*Init the laplacian matrix(system_matrix)*/
	void initMatrix();
	void solve(int iters, Edges edgelist, float *shortest_path);
	/*Set the laplacian matrix(system_matrix). */
	void setSystemMatrix();
	/*Re-init the laplacian matrix(system matrix)
	* when new constraints are added to the system. */
	void reInitMatrix(vector<Constraint*> newExtraConsts);
	void reInitMatrix(Constraint* newExtraConst);
	void reset();
	void initCuda();
	/*Edit the weight of constraints. */
	void EditConstraintWeight(int unit_const_begin,
		int unit_const_num, float new_weight);


	/*If it is the first time runing*/
	bool firstSetSysm = true, firstInitCuda = true;
	/*The result position of nodes*/
	VectorXf P_Opt_x, P_Opt_y;
	/*Laplacian matrix*/
	MatrixXf system_matrix;
	float **dists = 0;
	VectorXf right_hand_x, right_hand_y;
	std::vector<Constraint* > constraints;
	std::vector<Constraint* > extra_constraints;
	CudaSolver* cuda_solver;


private:

	Solver(const Solver&);
	void operator = (const Solver&);
};

#endif


