
#include "stdafx.h"
#include "NonOverlapConst.h"
#include "Solver.cuh" 
#include <iterator>

#include <iostream>

NonOverlapConst::NonOverlapConst(){
	this->init();
}

NonOverlapConst::~NonOverlapConst(){
}
void NonOverlapConst::init(){
	// init necessary variable here
	this->r_length = 1.0;
	this->d_vector = VectorXf::Zero(2);
}
void NonOverlapConst::initShape(Edge edge, float r_length, float width, float height, float wi,  int const_id){
	this->edge = edge;
	this->r_length = r_length;
	pid.push_back(this->edge.first);
	pid.push_back(this->edge.second);
	pid.shrink_to_fit();
	this->wi = wi;
	this->node_width = width;
	this->node_height = height;
	this->const_id = const_id;

 
}
int NonOverlapConst::Copy(void *linear_constraint) {
	cuDWORD *constraints = static_cast<cuDWORD *>(linear_constraint);
	constraints[0] = NonOverlapConst::type;
	constraints[1] = getdnum() * 2;
	memcpy(&constraints[2], &node_width, 4);
	memcpy(&constraints[3], &node_height, 4);
	constraints[4] = edge.first;
	constraints[5] = edge.second;
	float para = wi;
	memcpy(&constraints[6], &para, 4);
	memcpy(&constraints[7], &r_length, 4);
	return 8;
}// 7*4 bytes per constraints;
void NonOverlapConst::getRightHand(VectorXf& right_hand)
{
	right_hand = this->d_vector;
}
void NonOverlapConst::setLMatrix(MatrixXf &sysm)
{
	float k = this->wi / (this->r_length*this->r_length);
	//float k = this->wi / (this->r_length);
	sysm(this->edge.first, this->edge.first) += k;
	sysm(this->edge.first, this->edge.second) += -k;
	sysm(this->edge.second, this->edge.first) += -k;
	sysm(this->edge.second, this->edge.second) += k;
}
void NonOverlapConst::minusLMatrix(MatrixXf &sysm) {
	float k = this->wi / (this->r_length*this->r_length);
	//float k = this->wi / (this->r_length);
	sysm(this->edge.first, this->edge.first) += -k;
	sysm(this->edge.first, this->edge.second) += k;
	sysm(this->edge.second, this->edge.first) += k;
	sysm(this->edge.second, this->edge.second) += -k;
}
int NonOverlapConst::getdnum(){
	return 1;
}
void NonOverlapConst::setSolver(Solver* solver){
	this->solver = solver;
}
vector<int> NonOverlapConst::getpid(){
	return pid;
}
float NonOverlapConst::getWi(){
	return this->wi;
}
void NonOverlapConst::setWi(float wi){
	this->wi = wi;
}
int NonOverlapConst::getConstId() {
	return this->const_id;
}
void NonOverlapConst::setConstId(int const_id) {
	this->const_id = const_id;
}
void NonOverlapConst::setOuterExtraConstId(int outer_extra_const_id) {
	this->outer_extra_const_id = outer_extra_const_id;
}
int NonOverlapConst::getOuterExtraConstId() {
	return outer_extra_const_id;
}
Edges NonOverlapConst::getEdges(){
	Edges es;
	es.push_back(edge);
	return es;
}
void NonOverlapConst::setpid(vector<int> pid){
	this->pid = pid;
	this->edge.first = pid[0];
	this->edge.second = pid[1];
}
