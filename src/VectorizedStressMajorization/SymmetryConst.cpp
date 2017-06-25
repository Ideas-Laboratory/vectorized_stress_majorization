#include "stdafx.h"
#include "SymmetryConst.h"
#include "Solver.cuh" 
#include <iterator>
#include <iostream>

SymmetryConst::SymmetryConst(){
	this->init();
}

SymmetryConst::~SymmetryConst(){}
void SymmetryConst::init(){
}
void SymmetryConst::initShape(vector<pair<int, int>> closest_node_id_pairs, float* shortest_path, int ori_graph_n,
	float rotate_angle, float move_dis, float wi,  int const_id){
	this->closest_node_id_pairs = closest_node_id_pairs;
	for (int i = 0; i < closest_node_id_pairs.size(); i++) {
		this->pid.push_back(closest_node_id_pairs[i].first);
		this->pid.push_back(closest_node_id_pairs[i].second);
		for (int j = i + 1; j < closest_node_id_pairs.size(); j++) {
			this->edges_model.push_back(make_pair(closest_node_id_pairs[i].first, closest_node_id_pairs[j].first));
			this->edges.push_back(make_pair(closest_node_id_pairs[i].second, closest_node_id_pairs[j].second));
			float edge_x = solver->P_Opt_x[closest_node_id_pairs[i].first] - solver->P_Opt_x[closest_node_id_pairs[j].first];
			float edge_y = solver->P_Opt_y[closest_node_id_pairs[i].first] - solver->P_Opt_y[closest_node_id_pairs[j].first];
			float r_length = sqrt(edge_x *edge_x + edge_y*edge_y);
			this->r_lengths.push_back(r_length);
		}
	}
	this->d_vector = VectorXf::Zero(2 * this->edges.size());
	this->const_id = const_id;
	this->wi = wi;
	this->rotate_angle = rotate_angle;
	this->move_dis = move_dis;
}
int SymmetryConst::Copy(void *linear_constraint) {
	cuDWORD *constraints = static_cast<cuDWORD *>(linear_constraint);
	constraints[0] = SymmetryConst::type;
	constraints[1] = getdnum() * 2;
	memcpy(&constraints[2], &rotate_angle, 4);
	memcpy(&constraints[3], &move_dis, 4);
	constraints[4] = this->closest_node_id_pairs.size();//nodes number
	constraints[5] = this->closest_node_id_pairs.size()*(this->closest_node_id_pairs.size() - 1) / 2;//r_lengths number
	memcpy(&constraints[6], &wi, 4);
	int index = 7;
	/*model nodes*/
	for (int i = 0; i < closest_node_id_pairs.size(); i++) {
		constraints[index++] = closest_node_id_pairs[i].first;//model nodes
	}
	/*nodes in the other side*/
	for (int i = 0; i < closest_node_id_pairs.size(); i++) {
		constraints[index++] = closest_node_id_pairs[i].second;
	}
	/*rest lengths on every pair of nodes*/
	for (int i = 0; i < r_lengths.size(); i++) {
		memcpy(&constraints[index++], &r_lengths[i], 4);
	}
	return index;
}
void SymmetryConst::getRightHand(VectorXf& right_hand){
	right_hand = this->d_vector; 
}
void SymmetryConst::setLMatrix(MatrixXf &sysm){
	float k = 0;
	for (int i = 0; i < this->edges.size(); i++){
		k = this->wi / this->r_lengths[i];
		sysm(this->edges[i].first, this->edges[i].first) += k;
		sysm(this->edges[i].first, this->edges[i].second) += -k;
		sysm(this->edges[i].second, this->edges[i].first) += -k;
		sysm(this->edges[i].second, this->edges[i].second) += k;
	}
}
void SymmetryConst::minusLMatrix(MatrixXf &sysm) {
	float k = 0;
	for (int i = 0; i < this->edges.size(); i++){
		k = this->wi / this->r_lengths[i];
		sysm(this->edges[i].first, this->edges[i].first) += -k;
		sysm(this->edges[i].first, this->edges[i].second) += k;
		sysm(this->edges[i].second, this->edges[i].first) += k;
		sysm(this->edges[i].second, this->edges[i].second) += -k;
	}

}
int SymmetryConst::getdnum() {
	return r_lengths.size();
}
void SymmetryConst::setSolver(Solver* solver){
	this->solver = solver;
}
vector<int> SymmetryConst::getpid() {
	return pid;
}
float SymmetryConst::getWi() {
	return this->wi;
}
void SymmetryConst::setWi(float wi){
	this->wi = wi;
}
int SymmetryConst::getConstId() {
	return this->const_id;
}
void SymmetryConst::setConstId(int const_id) {
	this->const_id = const_id;
}
void SymmetryConst::setOuterExtraConstId(int outer_extra_const_id) {
	this->outer_extra_const_id = outer_extra_const_id;
}
int SymmetryConst::getOuterExtraConstId() {
	return outer_extra_const_id;
}
Edges SymmetryConst::getEdges(){
	return edges;
}
void SymmetryConst::setpid(vector<int> pid){
	this->pid = pid;
	for (int i = 0; i < pid.size() / 2; i++){
		this->closest_node_id_pairs.push_back(make_pair(pid[2 * i + 0], pid[2 * i + 1]));
	}
	this->edges_model.clear();
	this->edges.clear();
	for (int i = 0; i < closest_node_id_pairs.size(); i++) {
		for (int j = i + 1; j < closest_node_id_pairs.size(); j++) {
			this->edges_model.push_back(make_pair(closest_node_id_pairs[i].first, closest_node_id_pairs[j].first));
			this->edges.push_back(make_pair(closest_node_id_pairs[i].second, closest_node_id_pairs[j].second));
		}
	}
}
