#include "stdafx.h"
#include "StressConst.h"
#include "Solver.cuh" 
#include <iterator>

#include <iostream>

StressConst::StressConst(){
	this->init();
}

StressConst::~StressConst(){}
void StressConst::init(){
	// init necessary variable here
	this->r_length = 1.0;
	this->d_vector = VectorXf::Zero(2);
}
void StressConst::initShape(Edge edge, float r_length, pair<float, float> dire_vec, float wi, int const_id){
	this->edge = edge;
	this->r_length = r_length;
	pid.push_back(this->edge.first);
	pid.push_back(this->edge.second);
	pid.shrink_to_fit();
	this->wi = wi;
	this->dire_vec = dire_vec;
	this->const_id = const_id;
}
int StressConst::Copy(void *linear_constraint) {
	cuDWORD *constraints = static_cast<cuDWORD *>(linear_constraint);
	constraints[0] = getType();
	constraints[1] = getdnum() * 2;
	memcpy(&constraints[2], &dire_vec.first, 4);
	memcpy(&constraints[3], &dire_vec.second, 4);
	constraints[4] = edge.first;
	constraints[5] = edge.second;

	float para = wi;
	memcpy(&constraints[6], &para, 4);
	return 7;
}// 8*4 bytes per constraints;

void StressConst::getRightHand(VectorXf& right_hand){
	right_hand = this->d_vector; 
}
int StressConst::getdnum(){
	return 1;
}
void StressConst::setSolver(Solver* solver){
	this->solver = solver;
}
vector<int> StressConst::getpid(){
	return pid;
}
float StressConst::getWi(){
	return this->wi;
}
void StressConst::setWi(float wi) {
	this->wi = wi;
}
int StressConst::getConstId() {
	return this->const_id;
}
void StressConst::setConstId(int const_id) {
	this->const_id = const_id;
}
void StressConst::setOuterExtraConstId(int outer_extra_const_id) {
	this->outer_extra_const_id = outer_extra_const_id;
}
int StressConst::getOuterExtraConstId() {
	return outer_extra_const_id;
}
void StressConst::setLMatrix(MatrixXf &sysm, float **distMatrix = NULL) {
	if (distMatrix) {
		distMatrix[this->edge.first][this->edge.second] = r_length;//  1/w
		distMatrix[this->edge.second][this->edge.first] = r_length;//  1/w
	}
	float k = this->wi / (this->r_length*this->r_length);
	sysm(this->edge.first, this->edge.first) += k;
	sysm(this->edge.first, this->edge.second) += -k;
	sysm(this->edge.second, this->edge.first) += -k;
	sysm(this->edge.second, this->edge.second) += k;
}
void StressConst::setLMatrix(MatrixXf &sysm){
	setLMatrix(sysm, 0);
}
void StressConst::minusLMatrix(MatrixXf &sysm) {
	float k = this->wi / (this->r_length*this->r_length);
	sysm(this->edge.first, this->edge.first) += -k;
	sysm(this->edge.first, this->edge.second) += k;
	sysm(this->edge.second, this->edge.first) += k;
	sysm(this->edge.second, this->edge.second) += -k;
}
Edges StressConst::getEdges(){
	Edges es;
	es.push_back(edge);
	return es;
}