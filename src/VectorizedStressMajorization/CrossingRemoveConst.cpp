
#include "stdafx.h"
#include "CrossingRemoveConst.h"
#include "Solver.cuh" 
#include <iterator>
#include <iostream>

CrossingRemoveConst::CrossingRemoveConst(){
	this->init();
}
CrossingRemoveConst::~CrossingRemoveConst(){}
void CrossingRemoveConst::init(){}
void CrossingRemoveConst::initShape(Edges edges, vector<float> r_lengths, float wi, vector<int> degrees, int const_id){
	this->const_id = const_id;
	this->edges = edges;
	this->r_lengths = r_lengths;
	this->d_vector = VectorXf::Zero(2 * edges.size());
	/*the last element in degrees is the total degree of all nodes in edges*/
	for (int i = 0; i < degrees.size()-1; i++){
		this->edge_weights.push_back((float)degrees[i] / (float)degrees[degrees.size() - 1]);
	}

	for (int i = 0; i < edges.size(); i++){
		this->pid.push_back(edges[i].first);
		this->pid.push_back(edges[i].second);
	}
	this->pid.shrink_to_fit();
	this->wi = wi;
}
int CrossingRemoveConst::Copy(void * linear_constraint){
	cuDWORD *constraints = static_cast<cuDWORD *>(linear_constraint);
	constraints[0] = CrossingRemoveConst::type;
	constraints[1] = getdnum() * 2;
	constraints[2] = edges.size();
	memcpy(&constraints[3], &wi, 4);
	int idx = 4;
	/*edge*/
	for (int i = 0; i < this->edges.size(); i++) {
		constraints[idx++] = this->edges[i].first;
		constraints[idx++] = this->edges[i].second;
	}
	/*rest length*/
	float tmp;
	for (int i = 0; i < this->edges.size(); i++) {
		tmp = this->r_lengths[i]; 
		memcpy(&constraints[idx++], &tmp, 4);
	}
	/*edge weights*/
	for (int i = 0; i < this->edges.size(); i++){
		memcpy(&constraints[idx++], &this->edge_weights[i], 4);
	}
	return idx;
}
void CrossingRemoveConst::getRightHand(VectorXf& right_hand){
	right_hand = this->d_vector; 
}
void CrossingRemoveConst::setSolver(Solver* solver){
	this->solver = solver;
}
vector<int> CrossingRemoveConst::getpid(){
	return pid;
}
float CrossingRemoveConst::getWi(){
	return this->wi;
}
void CrossingRemoveConst::setWi(float wi){
	this->wi = wi;
}
void CrossingRemoveConst::setLMatrix(MatrixXf &sysm){
	float k = 0;
	for (int i = 0; i < edges.size(); i++){
		float k = edge_weights[i] * edges.size() * wi / (r_lengths[i] * r_lengths[i]);
		sysm(edges[i].first, edges[i].first) += k;
		sysm(edges[i].first, edges[i].second) += -k;
		sysm(edges[i].second, edges[i].first) += -k;
		sysm(edges[i].second, edges[i].second) += k;
	}
}
void CrossingRemoveConst::minusLMatrix(MatrixXf &sysm) {
	float k = 0;
	for (int i = 0; i < edges.size(); i++){
		float k = edge_weights[i] * edges.size() * wi / (r_lengths[i] * r_lengths[i]);
		sysm(edges[i].first, edges[i].first) += -k;
		sysm(edges[i].first, edges[i].second) += k;
		sysm(edges[i].second, edges[i].first) += k;
		sysm(edges[i].second, edges[i].second) += -k;
	}
}
int CrossingRemoveConst::getdnum(){
	return this->edges.size();
}
int CrossingRemoveConst::getConstId() {
	return const_id;
}
void CrossingRemoveConst::setConstId(int const_id) {
	this->const_id = const_id;
}
void CrossingRemoveConst::setOuterExtraConstId(int outer_extra_const_id) {
	this->outer_extra_const_id = outer_extra_const_id;
}
int CrossingRemoveConst::getOuterExtraConstId() {
	return outer_extra_const_id;
}
Edges CrossingRemoveConst::getEdges(){
	return edges;
}
void CrossingRemoveConst::setpid(vector<int> pid){
	this->pid = pid;
	this->edges.clear();
	for (int i = 0; i < this->pid.size() / 2; i++){
		this->edges.push_back(make_pair(this->pid[2 * i + 0], this->pid[2 * i + 1]));
	}
}