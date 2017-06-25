#include "stdafx.h"
#include "../VectorStressMajorization/VectorStressMajorization.h"

/*Mouse click action*/
void mouseClick(int btn, int state, int x, int y);
vector<int> hightlight_nodes;
vector<pair<int, int>> hightlight_edges;
/*Center nodes for Equal Angle Constraints*/
vector<int> parents;

int main(int argc, char* argv[])
{
	graph = new Graph;
	graph->GraphInitFromFile("fpga_dcop1220");

	graph->SolverInit();
	graph->SetupInitialConstraints();
	graph->SolveStress(25);

	int constraint_id = 0;

	graph->solver->initCuda();
	graph->solver->solve(50, graph->edgelist, graph->distances);

	/*Community Non-overlap Constraints*/
	/*The boundaries of Minkowski Differentce of every pair of A and B
	* Where A and B are different communities(clusters)*/
	vector<vector<float>>  minusBoundaries = getMinusBoundaries(
		graph->communities, graph->solver->P_Opt_x, graph->solver->P_Opt_y);
	/*The minimal penetration depths*/
	vector<pair<float, float>> mpds = nearestToOrigins(minusBoundaries);
	vector<Constraint*> ClusterNOConstraints = createClusterNOConsts(
		graph->communities, graph->solver, mpds, 2.0f,
		graph->extra_outer_constraint_include_inner_nums.size(),
		constraint_id);
	graph->extra_outer_constraint_include_inner_nums.push_back(ClusterNOConstraints.size());
	constraint_id += ClusterNOConstraints.size();

	/*Equal Angle Constraint*/
	parents = { 1152, 214, 215, 216, 217,218, 219, 220, 221 };
	hightlight_nodes = parents;
	vector<Constraint*> equal_angle_constraints;
	for (int i = 0; i < parents.size(); i++) {
		vector<int> children;
		float weight = 30.0f;
		if (parents[i] == 1152) weight = 80.0f;
		for (int j = 0; j < graph->n; j++) {
			if (findDistance(graph->distances, graph->n, parents[i], j) == 1) {
				hightlight_nodes.push_back(j);
				hightlight_edges.push_back(make_pair(parents[i], j));
				children.push_back(j);
			}
		}
		EqualAngleConst *EAC = new EqualAngleConst();
		EAC->setSolver(graph->solver);
		EAC->initShape(parents[i], children, 1.0f, weight, constraint_id);
		constraint_id++;
		equal_angle_constraints.push_back(EAC);
	}
	graph->solver->reInitMatrix(equal_angle_constraints);

	/*Circle Constraints*/
	/*Circles*/
	vector<vector<int>> circles 
		= { {3,233,1,236,1087,235,5,245,1088,243 },
			{41,417,1100,416,46,400,1099,397},
			{577,80,574,74,571,1111} };
	vector<Constraint*> circle_constraints;
	for (int i = 0; i < circles.size(); i++) {
		/*The idea distances of adjacent nodes on the ring*/
		vector<float> r_lengths;
		float weight = 40.0f;
		for (int j = 0; j < circles[i].size(); j++) {
			hightlight_nodes.push_back(circles[i][j]);
			hightlight_edges.push_back(make_pair(circles[i][j],
						circles[i][(j + 1) % circles[i].size()]));
			r_lengths.push_back(1.0f);
		}
		CircleConst *CC = new CircleConst();
		CC->setSolver(graph->solver);
		CC->initShape(circles[i], r_lengths, weight, constraint_id);
		circle_constraints.push_back(CC);
	}
	graph->solver->reInitMatrix(circle_constraints);
	

	if (ClusterNOConstraints.size() != 0) {
		graph->solver->reInitMatrix(ClusterNOConstraints);
	}
	graph->solver->initCuda();
	graph->solver->solve(20, graph->edgelist, graph->distances);
	graph->FitScreen();

	//opengl
	glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("Large Graph");
	glutDisplayFunc(renderScene);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glutMouseFunc(mouseClick);
	glutMainLoop();

	system("pause");
	return 0;
}
void renderScene(void) {

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glMatrixMode(GL_MODELVIEW);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glLoadIdentity();
	glBegin(GL_COLOR_BUFFER_BIT);


	for (int i = 0; i < graph->edgelist.size(); i++)
	{
		glColor4f(0.78f, 0.78f, 0.78f, 0.7f);

		float pointSource[2] = { graph->solver->P_Opt_x[graph->edgelist[i].first], 
			graph->solver->P_Opt_y[graph->edgelist[i].first] };
		float pointTarget[2] = { graph->solver->P_Opt_x[graph->edgelist[i].second], 
			graph->solver->P_Opt_y[graph->edgelist[i].second] };

		glLineWidth(1.5f);
		glBegin(GL_LINES);
		glVertex2f(pointSource[0], pointSource[1]);
		glVertex2f(pointTarget[0], pointTarget[1]);
		glEnd();

	}

	for (int i = 0; i < graph->solver->P_Opt_x.size(); i++) {
		glColor4f(0.5, 0.5, 0.5, 1.0f);

		for (int u = 0; u < graph->communities.size(); u++) {
			for (int y = 0; y < graph->communities[u].size(); y++) {
				if (graph->communities[u][y] == i) {
					glColor4f(colors[3 * u + 0], 
							colors[3 * u + 1],
							colors[3 * u + 2], 1.0f);
					break;
				}
			}
		}

		glPointSize(8);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[i], graph->solver->P_Opt_y[i]);
		glEnd();
	}


	/*highlight edges*/
	for (int i = 0; i < hightlight_edges.size(); i++) {
		glColor4f(1.0f, 0.0f, 0.0f, 0.6f);
		float pointSource[2] = { graph->solver->P_Opt_x[hightlight_edges[i].first], 
			graph->solver->P_Opt_y[hightlight_edges[i].first] };
		float pointTarget[2] = { graph->solver->P_Opt_x[hightlight_edges[i].second], 
			graph->solver->P_Opt_y[hightlight_edges[i].second] };

		glLineWidth(2.0f);
		glBegin(GL_LINES);
		glVertex2f(pointSource[0], pointSource[1]);
		glVertex2f(pointTarget[0], pointTarget[1]);
		glEnd();
	}
	/*hightlight nodes*/
	for (int i = 0; i < hightlight_nodes.size(); i++) {
		glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
		glPointSize(11);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[hightlight_nodes[i]], 
			graph->solver->P_Opt_y[hightlight_nodes[i]]);
		glEnd();

		for (int u = 0; u < graph->communities.size(); u++) {
			for (int y = 0; y < graph->communities[u].size(); y++) {
				if (graph->communities[u][y] == hightlight_nodes[i]) {
					glColor4f(colors[3 * u + 0], colors[3 * u + 1], colors[3 * u + 2], 1.0f);
					break;
				}
			}
		}
		glPointSize(8);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[hightlight_nodes[i]],
			graph->solver->P_Opt_y[hightlight_nodes[i]]);
		glEnd();
	}

	glutSwapBuffers();
}

void mouseClick(int btn, int state, int x, int y) {
	float _x = (float)x / 500.0f - 1.0f;
	float _y = -((float)y / 500.0f - 1.0f);
	if (btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		for (int i = 0; i < graph->solver->P_Opt_x.size(); i++)
		{
			float	xp = graph->solver->P_Opt_x[i];
			float	yp = graph->solver->P_Opt_y[i];
			if (_x < xp + 0.02f && _x > xp - 0.02f)
			{
				if (_y < yp + 0.02f && _y > yp - 0.02f)
				{
					cout << "click on =======  " << i << endl;
				}
			}
		}
	}

}