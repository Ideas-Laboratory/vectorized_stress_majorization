#include "Solver.cuh"
#include "Constraint.h" 
#include "stdafx.h"
#include <set>
vector<Constraint*>back;
int cnt = 0;
Solver::Solver()
{
	cuda_solver = new CudaSolver();
}

Solver::~Solver()
{
	delete cuda_solver;
	if (linear_constraints != 0){
		delete linear_constraints;
		linear_constraints = 0;
	}
}

void Solver::addConstraint(Constraint* constraint)
{
	this->constraints.push_back(constraint);
}
void Solver::addExtraConstraint(Constraint* constraint)
{
	this->extra_constraints.push_back(constraint);
}
void Solver::clearConstraints(){
	this->constraints.clear();
}

void Solver::initMatrix()
{
	if (this->constraints.empty())
	{
		std::cout << "empty constraints" << endl;
		exit(4);
	}
	std::cout << "prepare system matrix.\n";
	this->setSystemMatrix();
}

void Solver::solve(int iters, Edges edgelist, float *shortest_path)
{
	VectorXi es = VectorXi::Zero(edgelist.size() * 2);
	for (int i = 0; i < es.size() / 2; i++){
		es[2 * i + 0] = edgelist[i].first;
		es[2 * i + 1] = edgelist[i].second;
	}
	cuda_solver->Solve(this->right_hand_x.data(), this->right_hand_y.data(),
		this->P_Opt_x.data(), this->P_Opt_y.data(), iters, es.data(), shortest_path);


	return;
}

void Solver::setSystemMatrix()
{
	if (this->constraints.empty())
	{
		cout << "NO CONSTRAINTS" << endl;
		exit(2);
	}
	if (firstSetSysm)
	{
		linear_constraints = new cuDWORD[2009715200];
		this->system_matrix
			= MatrixXf::Zero(this->P_Opt_x.size(), this->P_Opt_x.size());

		if (dists) {
			delete[] dists[0];
			delete[] dists;
		}
		int n_data = this->P_Opt_x.size();
		dists = new float*[n_data];
		dists[0] = new float[this->P_Opt_x.size() * this->P_Opt_x.size()];
		std::fill(dists[0], dists[0] + n_data * n_data, 0);
		system_matrix(0, 0) += 1;//the anchor
		firstSetSysm = false;
	}
	else {
		std::fill(this->system_matrix.data(),
			this->system_matrix.data() + this->P_Opt_x.size()*this->P_Opt_x.size(), 0);
	}
	for (int i = 1; i < P_Opt_x.size(); i++)
	{
		dists[i] = dists[i - 1] + P_Opt_x.size();
	}

	for (decltype(this->constraints.size()) i = 0; i < this->constraints.size(); ++i)
	{
		if (this->constraints[i]->getType() == Stubs::STRESSCONSTRAINTS ||
			firstSetSysm)
			((StressConst*)this->constraints[i])->setLMatrix(system_matrix, dists);
		else
			this->constraints[i]->setLMatrix(system_matrix);
	}

	//Extra constraints 
	if (!extra_constraints.empty()){
		//L+L'
		for (decltype(this->extra_constraints.size()) i = 0;
			i < this->extra_constraints.size(); ++i)
		{
			this->extra_constraints[i]->setLMatrix(system_matrix);
		}
	}

	this->initCuda();
}
void Solver::reInitMatrix(vector<Constraint*> newExtraConsts) {
	//L+L'
	for (decltype(newExtraConsts.size()) i = 0; i < newExtraConsts.size(); ++i)
	{
		//L+L'
		newExtraConsts[i]->setLMatrix(system_matrix);
	}
	for (int i = 0; i < newExtraConsts.size(); i++) {
		this->extra_constraints.push_back(newExtraConsts[i]);
	}
}
void Solver::reInitMatrix(Constraint* newExtraConst) {
	//L+L'
	newExtraConst->setLMatrix(system_matrix);
	this->extra_constraints.push_back(newExtraConst);
	cout << "constraint" << endl;
}
void Solver::EditConstraintWeight(int unit_const_begin,
	int unit_const_num, float new_weight) {
	for (int i = unit_const_begin; i <unit_const_begin + unit_const_num; i++) {
		//L- these constraints L
		this->extra_constraints[i]->minusLMatrix(system_matrix);

		extra_constraints[i]->setWi(new_weight);
		this->extra_constraints[i]->setLMatrix(system_matrix);

	}

}

void Solver::reset() {
	for (decltype(this->extra_constraints.size()) i = 0;
		i < this->extra_constraints.size(); ++i)
	{
		if (this->constraints[i]->getType() == Stubs::STRESSCONSTRAINTS)
		{

		}
		else
		{
			this->extra_constraints[i]->minusLMatrix(system_matrix);
		}
	}
	this->extra_constraints.clear();
	this->initCuda();
}

void Solver::initCuda() {
	int *idx_lc = new int[extra_constraints.size() + 1];
	int n_ilc = 0;

	cuDWORD *lc_ptr = (cuDWORD*)linear_constraints;
	int i = 0;
	for (auto var : extra_constraints)
	{
		i++;
		idx_lc[n_ilc++] = lc_ptr - linear_constraints;
		lc_ptr += var->Copy(lc_ptr);

	}
	idx_lc[extra_constraints.size()] = lc_ptr - linear_constraints;

	cuda_solver->init(
		system_matrix.data(), this->system_matrix.rows(),
		linear_constraints, idx_lc,
		extra_constraints.size(),
		dists
		);

}
