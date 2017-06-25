
#ifndef Constraint_H
#define Constraint_H

#include "BasicHeader.h"
#include "cuSolver.cuh"
//#include "CrossingRemovalConst_cuProxy.cuh"
class Solver;

class Constraint
{
public:
	Constraint();
	virtual ~Constraint();

	// the interface for solver
	virtual void init() = 0;
	virtual int Copy(void *) = 0;
	virtual void getRightHand(VectorXf& right_hand) = 0;
	virtual int getType() { return type; }
	virtual void setSolver(Solver* solver) = 0;
	virtual vector<int> getpid() = 0;
	virtual float getWi() = 0;
	virtual int getdnum() = 0;
	virtual int  getConstId() = 0;
	virtual void  setConstId(int const_id) = 0;
	virtual void setOuterExtraConstId(int outer_extra_const_id) = 0;
	virtual int getOuterExtraConstId() = 0;
	virtual void setWi(float wi) = 0;
	virtual void setLMatrix(MatrixXf &sysm) = 0;
	virtual void minusLMatrix(MatrixXf &sysm) = 0;
	virtual Edges getEdges() = 0;

	int const_type = -1;
	static const int type = -1;// undefined

private:
	Constraint(const Constraint&); // not implemented
	void operator = (const Constraint&); // not implemented

};

#endif
