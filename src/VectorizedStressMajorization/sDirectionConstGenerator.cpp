
#include "sDirectionConstGenerator.h"
#include "Solver.cuh" 

sDirectionConstGenerator::sDirectionConstGenerator()
{
	this->init();
	//std::cout << "create a direction constraint.\n";
}

sDirectionConstGenerator::~sDirectionConstGenerator()
{

}
void sDirectionConstGenerator::init()
{
	this->wi = 1.0f;
}


void sDirectionConstGenerator::initShape(Edge edge, float edge_length, pair<float, float> dire_vector, float weightScale, int const_id, int outer_extra_const_id)
{
	this->vertex_ids = vertex_ids;
	this->weightScale = weightScale;
	this->dire_vector = dire_vector;
	this->edge_length = edge_length;
	this->edge = edge;
	this->const_id = const_id;
	this->outer_extra_const_id = outer_extra_const_id;
}

AtomConst* sDirectionConstGenerator::generateStressConst(bool isbigger){

	vector<AtomConst*> SCs;

	AtomConst *SC = new AtomConst();
	SC->setSolver(this->solver);

	if (isbigger){
		SC->initShape(edge, this->edge_length, this->dire_vector, this->weightScale, const_id);
	}
	else{
		SC->initShape(edge, this->edge_length, this->dire_vector, this->weightScale, const_id);
	}
	SC->setOuterExtraConstId(outer_extra_const_id);

	return SC;
}

void sDirectionConstGenerator::setSolver(Solver* solver)
{
	this->solver = solver;

}