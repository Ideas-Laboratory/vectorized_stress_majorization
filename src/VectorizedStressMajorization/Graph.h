#ifndef Graph_H
#define Graph_H 
#include "BasicHeader.h"
#include "Tool.h"  
#include "Solver.cuh"
#include "sDirectionConstGenerator.h"
#include "ShortestPath.h"

/*Generate the stress constraints*/
vector<StressConst*> createStressConsts(Edges edges, float* shortest_path, 
	Solver* solver, int outer_extra_const_id, int node_num);
/*Generate Direction Constraints for the edges in edgelist, 
*The idea direction is dire_vec, with the weight = weightScale*/
vector<Constraint*> createDireConsts(Edges edgelist, float* distances,
	Solver *solver, pair <float, float> dire_vec, float weightScale,
	int outer_extra_const_id, int node_num);
/*Generate the Community(cluster) Non-overlap Constraints
* Which consists of Atom Constraints*/
vector<Constraint*> createClusterNOConsts(vector<vector<int>> communities, 
	Solver *solver, vector<pair <float, float>> mpds, float weightScale,
	int outer_extra_const_id, int c_index);

class Graph
{
public:
	Graph();
	~Graph();
	/*Read the edges, node names, distances, circles, communities from file
	* Where edges and distances file is prerequisite
	* All the files must be named with some rules:
	* They are contained in a folder named filename,
	* edge file: filename.txt
	* distances file: matrix_filename.txt
	* circles file: filename_circles.txt
	* communities file: filename_communities.txt
	* node names file: filename_namelist.txt*/
	bool GraphInitFromFile(string filename);
	/*Init the solver and the coordinates*/
	void SolverInit();
	/*Initial stress constraints*/
	void SetupInitialConstraints();
	/*Solve the graph with stress constraints*/
	void SolveStress(int iters);
	/*Make the graph fit the screen*/
	void FitScreen();
	/*Init adjacent matrix and edgelist_all, 
	* which contains all the pairs of nodes*/
	void InitMatrix();

	/*Number of nodes in the graph*/
	int n;
	/*The shortest paths in graph */
	float *distances;
	/*The actual edges */
	Edges edgelist;
	/*All node pairs, size = n*(n-1)/2 */
	Edges edgelist_all;
	/*(i, j) = 1, node i and j is adjacent
	  (i, j) = 0, otherwise */
	MatrixXf adjMatrix;
	/*Circles read from file*/
	vector<vector<int>> circles_all;
	/*Communities read from file*/
	vector<vector<int>> communities;
	Solver *solver;
	/*Basic stress constraints*/
	vector<StressConst*> StressConsts;
	/*The name of nodes*/
	vector<string> name_list;
	/*The number of unit constraints in a complex constraint */
	vector<int> extra_outer_constraint_include_inner_nums;
};

#endif


