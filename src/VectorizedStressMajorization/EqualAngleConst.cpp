
#include "stdafx.h"
#include "EqualAngleConst.h"
#include "Solver.cuh" 
#include <iterator>

#include <iostream>

EqualAngleConst::EqualAngleConst(){
	this->init();
}
EqualAngleConst::~EqualAngleConst(){}
void EqualAngleConst::init(){ }
void EqualAngleConst::initShape(int parentid, vector<int> child_list, float plength, float wi,  int const_id){
	this->parentid = parentid;
	this->child_list = child_list;
	this->pid = child_list;
	this->pid.push_back(parentid);

	for (int i = 0; i < this->child_list.size(); i++){
		this->edges.push_back(make_pair(this->child_list[i], parentid));
	}

	this->plength = plength;
	this->wi = wi;
	this->const_id = const_id;

	this->differ_angle = 2 * M_PI / this->child_list.size();
	this->anchor_edge = this->edges[0]; // initial set
	this->d_vectors = VectorXf::Zero(2 * this->child_list.size());
	relative_shape = this->getLocalRelativeShape();
}
void EqualAngleConst::getRightHand(VectorXf& right_hand){
	right_hand = this->d_vectors;
}
void EqualAngleConst::setLMatrix(MatrixXf &sysm){
	float k = 0;
	k = this->wi / (plength*plength);
	for (int i = 0; i < this->edges.size(); i++){
		sysm(this->edges[i].first, this->edges[i].first) += k;
		sysm(this->edges[i].first, this->edges[i].second) += -k;
		sysm(this->edges[i].second, this->edges[i].first) += -k;
		sysm(this->edges[i].second, this->edges[i].second) += k;
	}
}
void EqualAngleConst::minusLMatrix(MatrixXf &sysm) {
	float k = 0;
	k = this->wi / (plength*plength);
	for (int i = 0; i < this->edges.size(); i++){
		sysm(this->edges[i].first, this->edges[i].first) += -k;
		sysm(this->edges[i].first, this->edges[i].second) += k;
		sysm(this->edges[i].second, this->edges[i].first) += k;
		sysm(this->edges[i].second, this->edges[i].second) += -k;
	}
}
void EqualAngleConst::setSolver(Solver* solver){
	this->solver = solver;
}
vector<int> EqualAngleConst::getpid(){
	return pid;
}
int EqualAngleConst::getdnum(){
	return this->edges.size();
}
float EqualAngleConst::getWi(){
	return this->wi;
}
void EqualAngleConst::setWi(float wi){
	this->wi = wi;
}
int EqualAngleConst::Copy(void * linear_constraint){
	cuDWORD *constraints = static_cast<cuDWORD *>(linear_constraint);
	constraints[0] = EqualAngleConst::type;
	constraints[1] = getdnum() * 2;
	memcpy(&constraints[2], &differ_angle, 4);
	constraints[3] = pid.size();
	constraints[4] = edges.size(); //=angles size
	constraints[5] = plength;
	float param = wi / (float)plength;
	memcpy(&constraints[6], &param, 4);
	int idx = 7;
	for (int i = 0; i < this->pid.size(); i++) {
		constraints[idx++] = this->pid[i];
	}
	//relative shape
	float tmp;
	for (int i = 0; i < this->relative_shape.size(); i++){
		tmp = relative_shape[i];
		memcpy(&constraints[idx++], &tmp, 4);
	}

	return pid.size() + relative_shape.size() + 7;
}
int EqualAngleConst::getConstId() {
	return const_id;
}
void EqualAngleConst::setConstId(int const_id) {
	this->const_id = const_id;
}
void EqualAngleConst::setOuterExtraConstId(int outer_extra_const_id) {
	this->outer_extra_const_id = outer_extra_const_id;
}
int EqualAngleConst::getOuterExtraConstId() {
	return outer_extra_const_id;
}
Edges EqualAngleConst::getEdges(){
	return edges;
}
void EqualAngleConst::setpid(vector<int> pid){
	this->pid = pid;

	this->parentid = pid[pid.size()-1];
	this->child_list.clear();
	for (int i = 0; i < pid.size()-1; i++){
		this->child_list.push_back(pid[i]);
	}
	this->edges.clear();
	for (int i = 0; i < this->child_list.size(); i++){
		this->edges.push_back(make_pair(this->child_list[i], this->parentid));
	}
}
vector<float> EqualAngleConst::getLocalRelativeShape(){
	std::vector<float> relativePositions;

	for (int i = 0; i < child_list.size(); i++){
		float rotate_angle = (i)*differ_angle;
		float edge_length = 1.0f;
		float re_x = edge_length * cos(rotate_angle);// +relativePositions[relativePositions.size() - 2];// this->solver->P_Opt[2 * edges[i].first + 0];
		float re_y = edge_length * sin(rotate_angle);// +relativePositions[relativePositions.size() - 1];

		if (abs(re_x) < 0.0000001){
			re_x = 0;
		}
		if (abs(re_y) < 0.0000001){
			re_y = 0;
		}
		relativePositions.push_back(re_x);
		relativePositions.push_back(re_y);
	}

	//the center node:
	relativePositions.push_back(0.0f);
	relativePositions.push_back(0.0f);
	relativePositions.shrink_to_fit();

	return relativePositions;
}