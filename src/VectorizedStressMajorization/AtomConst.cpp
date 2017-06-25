#include "stdafx.h"
#include "AtomConst.h"
#include "Solver.cuh" 
#include <iterator>

#include <iostream>


AtomConst::AtomConst(){
	this->init();
}
AtomConst::~AtomConst(){}
void AtomConst::init(){
	// init necessary variable here
	this->r_length = 1.0;
	this->d_vector = VectorXf::Zero(2);
}
void AtomConst::initShape(Edge edge, float r_length, pair<float, float> dire_vec, float wi, int const_id){
	this->edge = edge;
	this->r_length = r_length;
	pid.push_back(this->edge.first);
	pid.push_back(this->edge.second);
	pid.shrink_to_fit();
	this->wi = wi;
	this->dire_vec = dire_vec;
	this->const_id = const_id;
}
int AtomConst::Copy(void *linear_constraint) {
	cuDWORD *constraints = static_cast<cuDWORD *>(linear_constraint);
	constraints[0] = getType();
	constraints[1] = getdnum() * 2;
	memcpy(&constraints[2], &dire_vec.first, 4);
	memcpy(&constraints[3], &dire_vec.second, 4);
	constraints[4] = edge.first;
	constraints[5] = edge.second;

	float para = wi / r_length;
	memcpy(&constraints[6], &para, 4);
	return 7;
}

void AtomConst::getRightHand(VectorXf& right_hand){
	right_hand = this->d_vector; //(1 / (this->wi * this->wi * 2))* 
}
int AtomConst::getdnum(){
	return 1;
}
void AtomConst::setSolver(Solver* solver){
	this->solver = solver;

}
vector<int> AtomConst::getpid(){
	return pid;
}
float AtomConst::getWi(){
	return this->wi;
}
void AtomConst::setWi(float wi) {
	this->wi = wi;
}
int AtomConst::getConstId() {
	return this->const_id;
}
void AtomConst::setConstId(int const_id) {
	this->const_id = const_id;
}
void AtomConst::setOuterExtraConstId(int outer_extra_const_id) {
	this->outer_extra_const_id = outer_extra_const_id;
}
int AtomConst::getOuterExtraConstId() {
	return outer_extra_const_id;
}
void AtomConst::setLMatrix(MatrixXf &sysm){
	float k = this->wi / (this->r_length);
	sysm(this->edge.first, this->edge.first) += k;
	sysm(this->edge.first, this->edge.second) += -k;
	sysm(this->edge.second, this->edge.first) += -k;
	sysm(this->edge.second, this->edge.second) += k;
}
void AtomConst::minusLMatrix(MatrixXf &sysm) {
	float k = this->wi / (this->r_length);
	sysm(this->edge.first, this->edge.first) += -k;
	sysm(this->edge.first, this->edge.second) += k;
	sysm(this->edge.second, this->edge.first) += k;
	sysm(this->edge.second, this->edge.second) += -k;
}
Edges AtomConst::getEdges(){
	Edges es;
	es.push_back(edge);
	return es;
}