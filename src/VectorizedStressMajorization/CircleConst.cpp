
#include "stdafx.h"
#include "CircleConst.h"
#include "Solver.cuh" 
#include <iterator>

#include <iostream>

CircleConst::CircleConst(){
	this->init();
}
CircleConst::~CircleConst(){}

void CircleConst::init(){
	// init necessary variable here  
}
void CircleConst::initShape(vector<int> pid, vector<float> r_lengths, float wi, int const_id)
{
	this->pid = pid;

	for (int i = 1; i < this->pid.size(); i++){
		this->edges.push_back(make_pair(this->pid[i - 1], this->pid[i]));
	}
	this->edges.push_back(make_pair(this->pid[this->pid.size() - 1], this->pid[0]));

	this->r_lengths = r_lengths;
	this->wi = wi;
	this->const_id = const_id;

	this->differ_angle = M_PI - M_PI*(this->pid.size() - 2) / this->pid.size();
	this->anchor_edge = this->edges[0]; // initial set
	this->d_vectors = VectorXf::Zero(2 * this->pid.size());

	relative_shape = this->getLocalRelativeShape();
	
}


void CircleConst::getRightHand(VectorXf& right_hand){
	right_hand = this->d_vectors; 
}

void CircleConst::setLMatrix(MatrixXf &sysm){
	float k = 0;
	for (int i = 0; i < this->edges.size(); i++){
		float k = this->wi / (this->r_lengths[i] * this->r_lengths[i]);
		sysm(this->edges[i].first, this->edges[i].first) += k;
		sysm(this->edges[i].first, this->edges[i].second) += -k;
		sysm(this->edges[i].second, this->edges[i].first) += -k;
		sysm(this->edges[i].second, this->edges[i].second) += k;
	}
}
void CircleConst::minusLMatrix(MatrixXf &sysm) {
	float k = 0;
	for (int i = 0; i < this->edges.size(); i++){
		float k = this->wi / (this->r_lengths[i] * this->r_lengths[i]);
		sysm(this->edges[i].first, this->edges[i].first) += -k;
		sysm(this->edges[i].first, this->edges[i].second) += k;
		sysm(this->edges[i].second, this->edges[i].first) += k;
		sysm(this->edges[i].second, this->edges[i].second) += -k;
	}
}
void CircleConst::setSolver(Solver* solver){
	this->solver = solver;

}
vector<int> CircleConst::getpid(){
	return pid;
}
int CircleConst::getdnum(){
	return this->edges.size();
}
float CircleConst::getWi(){
	return this->wi;
}
void CircleConst::setWi(float wi){
	this->wi = wi;
}
int CircleConst::Copy(void * linear_constraint)
{
	cuDWORD *constraints = static_cast<cuDWORD *>(linear_constraint);
	constraints[0] = CircleConst::type;
	constraints[1] = getdnum() * 2;
	memcpy(&constraints[2], &differ_angle, 4);
	constraints[3] = pid.size();
	constraints[4] = r_lengths.size();
	memcpy(&constraints[5], &wi, 4);
	int idx = 6;
	//pids
	for (int i = 0; i < this->pid.size(); i++) {
		constraints[idx++] = this->pid[i];
	}

	float tmp;
	//r_lengths
	for (int i = 0; i < this->edges.size(); i++) {
		tmp = this->r_lengths[i];
		memcpy(&constraints[idx++], &tmp, 4);
	}
	//relative shapes
	for (int i = 0; i < this->relative_shape.size(); i++){
		tmp = relative_shape[i];
		memcpy(&constraints[idx++], &tmp, 4);
	}
	/*the length of this constriant*/
	return idx;
}
int CircleConst::getConstId() {
	return const_id;
}
void CircleConst::setConstId(int const_id) {
	this->const_id = const_id;
}
void CircleConst::setOuterExtraConstId(int outer_extra_const_id) {
	this->outer_extra_const_id = outer_extra_const_id;
}
int CircleConst::getOuterExtraConstId() {
	return outer_extra_const_id;
}
Edges CircleConst::getEdges(){
	return edges;
}
void CircleConst::setpid(vector<int> pid){
	this->pid = pid;
	this->edges.clear();
	for (int i = 1; i < this->pid.size(); i++){
		this->edges.push_back(make_pair(this->pid[i - 1], this->pid[i]));
	}
	this->edges.push_back(make_pair(this->pid[this->pid.size() - 1], this->pid[0]));
}
vector<float> CircleConst::getLocalRelativeShape(){
	std::vector<float> relativePositions; 
	relativePositions.push_back(1.0f);
	relativePositions.push_back(0.0f);

	for (int i = 0; i < edges.size() - 1; i++){
		float rotate_angle = (i + 1)*differ_angle;
		float edge_length = 1.0f;

		float re_x = edge_length * cos(rotate_angle);
		float re_y = edge_length * sin(rotate_angle);

		if (abs(re_x) < 0.0000001){
			re_x = 0;
		}
		if (abs(re_y) < 0.0000001){
			re_y = 0;
		}

		relativePositions.push_back(re_x);
		relativePositions.push_back(re_y);
	}

	relativePositions.shrink_to_fit();

	return relativePositions;
}