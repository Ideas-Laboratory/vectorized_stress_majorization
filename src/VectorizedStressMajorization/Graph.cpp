#include "Graph.h" 
#include "stdafx.h" 
#include "iostream"
#include "AtomConst.h"

Graph::Graph(){ }
Graph::~Graph(){ }
bool Graph::GraphInitFromFile(string graph_name)
{
	string filename = "../../data/" + graph_name; // Data path if alternative data location is used
	VertexNameList vertexnamelist;
	string filepath = filename + "/" + graph_name + ".txt";
	string str = "";
	ifstream fin(filepath);
	if (!fin){
		str = "No edge file called " + graph_name + ".txt";
		cout << str << endl;
		return 0;
	}
	fin.close();
	FromTxtEdgeFile(filepath, edgelist, n);

	filepath = filename + "/matrix_" + graph_name + ".txt";
	ifstream fin1(filepath);
	if (!fin1){
		str = "No distance matrix file called matrix_" + graph_name + ".txt";
		cout << str << endl;
		return 0;
	}
	fin1.close(); 
	distances = readAdjFromFileToShortestPath(filename + "/matrix_" + graph_name + ".txt", n);

	string communitiespath = filename + "/" + graph_name + "_communities.txt";
	ifstream fin3(communitiespath);
	if (!fin3){
		str = "No community file called " + graph_name + "_communities.txt";
		cout << str << endl;
	}
	else
	{
		communities = communitiesFromFile(communitiespath);
	}
	fin3.close();

	string circlesfilepath = filename + "/" + graph_name + "_circles.txt";
	ifstream fin2(circlesfilepath);
	if (!fin2){
		str = "No circles file called " + graph_name + "_circles.txt";
		cout << str << endl;
	}
	else
	{
		readCirclesFileToCircles(circlesfilepath, circles_all);
	}
	fin2.close();


	filepath = filename + "/" + graph_name + "_namelist.txt";
	ifstream fin4(filepath);
	if (!fin4) {}
	else {
		readNameList(filepath, name_list);
	}
	fin4.close();

	InitMatrix();
	return true;
}
void Graph::SolverInit() {
	solver = new Solver();
	randomCoordinates(n, solver->P_Opt_x, solver->P_Opt_y);
}
void Graph::SetupInitialConstraints() {
	cout << "setting up initial cosntriants......." << endl;
	StressConsts = createStressConsts(edgelist_all, distances, solver, -1, n);//edgelist_notInE edgelist_all
	cout << "create constraint done " << StressConsts.size() << endl;
	for (int i = 0; i < StressConsts.size(); i++) {
		 solver->addConstraint(StressConsts[i]);
	}
	cout << "setting up initial constraints 3" << endl;
	StressConsts.clear();
	edgelist_all.~vector();
	cout << "set up constraints done!" << endl;

}
void Graph::SolveStress(int iters) {	cout << "init matrix" << endl;
	time_t start_cho = clock();
	solver->initMatrix();
	time_t end_cho = clock();
	printf("the initing matrix time is : %f\n",		double(end_cho - start_cho) / CLOCKS_PER_SEC);	cout << "solving" << endl;
	time_t start = clock();
	solver->solve(iters, edgelist, 0);//distances
	time_t end = clock();
	printf("the solving time is : %f\n", double(end - start) / CLOCKS_PER_SEC);

	FitScreen();
}

void Graph::FitScreen() {
	Vector3f move = MovePoints(solver->P_Opt_x, solver->P_Opt_y);

	for (int i = 0; i < solver->P_Opt_x.size(); i++) {
		solver->P_Opt_x[i] = (solver->P_Opt_x[i] / move[2] + move[0])*0.9f;
		solver->P_Opt_y[i] = (solver->P_Opt_y[i] / move[2] + move[1])*0.9f;
	}
}

void Graph::InitMatrix() {
	adjMatrix = MatrixXf::Zero(n, n);
	for (int i = 0; i < n; i++){
		for (int j = i+1; j < n; j++){
			if (findDistance(distances, n, i, j) == 1){
				adjMatrix(i, j) = 1;
				adjMatrix(j, i) = 1;
			}
			edgelist_all.push_back(make_pair(i, j));
		}
	}
}
vector<StressConst*> createStressConsts(Edges edges, float* distances, 
	Solver* solver, int outer_extra_const_id, int node_num) {

	vector<StressConst*> constraints;

	size_t id_edge = 0;
	int count = 0;
	int const_id = 0;
	for (auto& i : edges)
	{
		StressConst *SC = new StressConst();

		SC->setSolver(solver);

		float rest_length = findDistance(distances, node_num, i.first, i.second);
		SC->initShape(i, rest_length, make_pair(0, 0), 1.0f, const_id);
		SC->setOuterExtraConstId(outer_extra_const_id);
		constraints.push_back(SC);
		const_id++;
	}
	return constraints;
}
vector<Constraint*> createDireConsts(Edges edgelist, float* distances, 
	Solver *solver, pair <float, float> dire_vec, float weightScale, 
	int outer_extra_const_id, int node_num) {
	vector<Constraint*> constraints;
	int const_id = 0;
	int d_idx = 0;
	for (auto& i : edgelist)
	{
		sDirectionConstGenerator *DCge = new sDirectionConstGenerator();

		float rest_length = findDistance(distances, node_num, i.first, i.second);

		DCge->setSolver(solver);
		DCge->initShape(i, rest_length, dire_vec, weightScale, const_id, outer_extra_const_id);
		AtomConst* SC = DCge->generateStressConst(false);
		constraints.push_back(SC);
		d_idx++;
		const_id++;
	}
	return constraints;
}
vector<Constraint*> createClusterNOConsts(vector<vector<int>> communities, Solver *solver,
vector<pair <float, float>> mpds, float weightScale, int outer_extra_const_id, int c_index) {
	vector<Constraint*> constraints;
	int const_id = 0;

	for (int i = 0; i < communities.size(); i++) {
		for (int j = i + 1; j < communities.size(); j++) {

			int mpds_index = i*(2 * communities.size() - i - 1) / 2 + (j - i - 1);

			if (mpds[mpds_index].first != 0 && mpds[mpds_index].second != 0){//if it is overlapping
				for (int m = 0; m < communities[i].size(); m++) {
					for (int n = 0; n < communities[j].size(); n++) {
						sDirectionConstGenerator *DCge = new sDirectionConstGenerator();
						DCge->setSolver(solver);
						Edge e = make_pair(communities[i][m], communities[j][n]);
						pair<float, float> ori_v = make_pair(solver->P_Opt_x[communities[j][n]] - solver->P_Opt_x[communities[i][m]], solver->P_Opt_y[communities[j][n]] - solver->P_Opt_y[communities[i][m]]);
						pair<float, float> dire = make_pair(mpds[mpds_index].first + ori_v.first, mpds[mpds_index].second + ori_v.second);

						dire.first = -dire.first;
						dire.second = -dire.second;
						float r_length = sqrtf(dire.first*dire.first + dire.second*dire.second);
						if (r_length == 0.0f){
							
							r_length = sqrtf(mpds[mpds_index].first*mpds[mpds_index].first + mpds[mpds_index].second*mpds[mpds_index].second);
							dire.first = mpds[mpds_index].first / r_length;
							dire.second = mpds[mpds_index].second / r_length;
						}
						else{
							dire.first /= r_length;
							dire.second /= r_length;
						}
						float weight = weightScale;
						r_length *= 5;
						DCge->initShape(e, (float)r_length, dire, weight, c_index + const_id, outer_extra_const_id);
						AtomConst* SC = DCge->generateStressConst(true);
						constraints.push_back(SC);
						const_id++;
					}
				}
			}

		}
	}

	return constraints;
}

