#include "stdafx.h"  
#include "../VectorizedStressMajorization/VectorizedStressMajorization.h"



/*Mouse click action*/
void mouseClick(int btn, int state, int x, int y);
/*Draw lasso*/
void drawLasso(vector<float> area);
/*Draw symmetry axis*/
void drawSymmetryAxis(VectorXf axis);
/*Generate Symmetry Constraint for the area_nodes with weight*/
vector<Constraint*> GenerateSymmetryConst(
	vector<int> area_nodes,
	double weight);
vector<Constraint*>  removeCrossing(
	vector<int> lasso_area_nodes, double weight);
/*Symmetry lasso area, constrained the nodes in this area to be symmetry*/
vector<float> symmetry_lasso_area;
/*Symemtry axis*/
VectorXf symmetry_axis = VectorXf::Zero(4);
/*Edge crossing lasso area, remove edge crossings in this area*/
vector<float> crossing_lasso_area;
Edges intersect_edges;
/*Closest node pairs,  in each pair:
* the first element is the template node, the second is the constrained node.
* Constrained nodes should be as similar as possible to the template nodes*/
vector<pair<int, int>> closest_node_id_pairs;

int _tmain(int argc, _TCHAR* argv[])
{
	graph = new Graph;
	graph->GraphInitFromFile("bus1138");

	graph->SolverInit();
	graph->SetupInitialConstraints();
	graph->SolveStress(50);
	graph->FitScreen();

	symmetry_lasso_area =
	{ 0.122f, 0.964f, -0.12f, 0.764f, -0.08f, 0.488f, 0.144f, 0.216f,
		0.578f, 0.116f, 0.718f, 0.164f, 0.776f, 0.262f,0.736f, 0.568f,
		0.662f, 0.858f, 0.33f, 0.954f };
	symmetry_axis[0] = 0.51f;
	symmetry_axis[1] = 0.806f;
	symmetry_axis[2] = 0.156f;
	symmetry_axis[3] = 0.272f;

	vector<int> symmetry_lasso_area_nodes;
	AddLassoNode(graph->solver->P_Opt_x, graph->solver->P_Opt_y,
		symmetry_lasso_area, symmetry_lasso_area_nodes);
	vector<Constraint*> symmetry_constraints
		= GenerateSymmetryConst(symmetry_lasso_area_nodes, 4.0f);

	if (symmetry_constraints.size() != 0) {
		graph->solver->reInitMatrix(symmetry_constraints);
	}

	/*Remove Edge Crossing Constriants*/
	vector<int> crossing_lasso_area_nodes;
	crossing_lasso_area =
	{ -0.282f,-0.154f,-0.362f,-0.136f,-0.428f,-0.166f,
		-0.446f,-0.214f,-0.404f,-0.248f,-0.298f,-0.24f,
		-0.258f,-0.168f };
	AddLassoNode(graph->solver->P_Opt_x, graph->solver->P_Opt_y,
		crossing_lasso_area, crossing_lasso_area_nodes);
	vector<Constraint*> remove_crossing_constraints = removeCrossing(crossing_lasso_area_nodes, 60.0f);
	if (remove_crossing_constraints.size() != 0) {
		graph->solver->reInitMatrix(remove_crossing_constraints);
	}


	/*To check the layout before adding Symmetry Constraints and
	* Edge Crossing Removing Constraints to the system,
	* note the following 3 rows code */
	graph->solver->initCuda();
	graph->solver->solve(5, graph->edgelist, graph->distances);
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
		glPointSize(8);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[i], graph->solver->P_Opt_y[i]);
		glEnd();
	}

	/*hightlight symmetry nodes.
	* Node in red: templates
	* Node in green: constrained nodes */
	for (int i = 0; i < closest_node_id_pairs.size(); i++) {
		glColor4f(0.8f, 0.5f, 0.8f, 1.0f);
		glPointSize(11);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[closest_node_id_pairs[i].first],
			graph->solver->P_Opt_y[closest_node_id_pairs[i].first]);
		glEnd();

		glColor4f(0.5, 0.5, 0.5, 1.0f);
		glPointSize(8);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[closest_node_id_pairs[i].first],
			graph->solver->P_Opt_y[closest_node_id_pairs[i].first]);
		glEnd();

		glColor4f(0.5f, 0.8f, 0.8f, 1.0f);
		glPointSize(11);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[closest_node_id_pairs[i].second],
			graph->solver->P_Opt_y[closest_node_id_pairs[i].second]);
		glEnd();

		glColor4f(0.5, 0.5, 0.5, 1.0f);
		glPointSize(8);
		glBegin(GL_POINTS);
		glVertex2f(graph->solver->P_Opt_x[closest_node_id_pairs[i].second],
			graph->solver->P_Opt_y[closest_node_id_pairs[i].second]);
		glEnd();
	}

	/*Highlight intersected edges*/
	for (int i = 0; i < intersect_edges.size(); i++) {
		glColor4f(0.7f, 0.0f, 0.0f, 0.4f);

		float pointSource[2] = { graph->solver->P_Opt_x[intersect_edges[i].first],
			graph->solver->P_Opt_y[intersect_edges[i].first] };
		float pointTarget[2] = { graph->solver->P_Opt_x[intersect_edges[i].second],
			graph->solver->P_Opt_y[intersect_edges[i].second] };
		glLineWidth(2.0f);
		glBegin(GL_LINES);
		glVertex2f(pointSource[0], pointSource[1]);
		glVertex2f(pointTarget[0], pointTarget[1]);
		glEnd();
	}

	drawLasso(symmetry_lasso_area);
	drawSymmetryAxis(symmetry_axis);
	drawLasso(crossing_lasso_area);

	glutSwapBuffers();
}

void mouseClick(int btn, int state, int x, int y) {
	float _x = (float)x / 500.0f - 1.0f;
	float _y = -((float)y / 500.0f - 1.0f);
	cout << _x << "," << _y << endl;
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
void drawLasso(vector<float> area) {
	glColor4f(0.0f, 0.0f, 0.5, 0.5f);
	glPointSize(1.5f);

	for (int i = 0; i < area.size() / 2; i++) {
		glBegin(GL_LINES);
		glVertex2f(area[2 * i + 0], area[2 * i + 1]);
		glVertex2f(area[2 * (i + 1) % area.size() + 0],
			area[2 * (i + 1) % area.size() + 1]);
		glEnd();
	}
}
void drawSymmetryAxis(VectorXf axis) {
	glColor4f(0.0f, 0.5f, 0.0f, 0.5f);
	glPointSize(1.5f);
	glBegin(GL_LINES);
	glVertex2f(axis[0], axis[1]);
	glVertex2f(axis[2], axis[3]);
	glEnd();
}
vector<Constraint*> GenerateSymmetryConst(vector<int> area_nodes, double weight) {
	vector<Constraint*> SyCs;
	float symmetry_dis_x = symmetry_axis[0] - symmetry_axis[2];
	float symmetry_dis_y = symmetry_axis[1] - symmetry_axis[3];
	/*Rotate and Move in order to find the closest node pairs.*/
	float rotate_angle = acosf((symmetry_dis_y)
		/ sqrt(symmetry_dis_x*symmetry_dis_x
			+ symmetry_dis_y*symmetry_dis_y));
	if (symmetry_axis[0] > symmetry_axis[2]) {
		rotate_angle = -rotate_angle;
	}
	_rotate(rotate_angle, graph->solver->P_Opt_x, graph->solver->P_Opt_y);
	_rotate(rotate_angle, symmetry_axis);
	float move_dis = symmetry_axis[0];
	moveCenterToZero(move_dis, graph->solver->P_Opt_x, graph->solver->P_Opt_y);
	moveCenterToZero(move_dis, symmetry_axis);
	/*Find the closest node pairs with greedy strategy.
	* Those who do not has a closest one in the other side will be eliminated.*/
	mapAreaRightToLeft(graph->solver->P_Opt_x, graph->solver->P_Opt_y,
		area_nodes, closest_node_id_pairs);

	if (closest_node_id_pairs.size() <= 1) {
		return SyCs;
	}
	int  constraint_id = graph->solver->extra_constraints.size();
	SymmetryConst *SyC = new SymmetryConst();
	SyC->setSolver(graph->solver);
	SyC->initShape(closest_node_id_pairs, graph->distances, graph->n,
		rotate_angle, move_dis, weight, constraint_id);

	moveCenterToZero(-move_dis, graph->solver->P_Opt_x, graph->solver->P_Opt_y);
	moveCenterToZero(-move_dis, symmetry_axis);
	_rotate(-rotate_angle, graph->solver->P_Opt_x, graph->solver->P_Opt_y);
	_rotate(-rotate_angle, symmetry_axis);

	SyCs.push_back(SyC);

	return SyCs;
}

vector<Constraint*>  removeCrossing(vector<int> lasso_area_nodes, double weight)
{
	vector<Constraint*> CRC0s;

	Edges es;
	for (int i = 0; i < lasso_area_nodes.size(); i++) {
		for (int j = i + 1; j < lasso_area_nodes.size(); j++) {
			if (findDistance(graph->distances, graph->n, lasso_area_nodes[i], lasso_area_nodes[j]) == 1) {
				es.push_back(make_pair(lasso_area_nodes[i], lasso_area_nodes[j]));
			}
		}
	}
	for (int i = 0; i < es.size(); i++) {
		for (int j = i + 1; j < es.size(); j++) {
			pair<float, float> edge1_source = make_pair(graph->solver->P_Opt_x[es[i].first], graph->solver->P_Opt_y[es[i].first]);
			pair<float, float> edge1_target = make_pair(graph->solver->P_Opt_x[es[i].second], graph->solver->P_Opt_y[es[i].second]);
			pair<float, float> edge2_source = make_pair(graph->solver->P_Opt_x[es[j].first], graph->solver->P_Opt_y[es[j].first]);
			pair<float, float> edge2_target = make_pair(graph->solver->P_Opt_x[es[j].second], graph->solver->P_Opt_y[es[j].second]);

			float	line1[4] = { edge1_source.first, edge1_source.second, edge1_target.first, edge1_target.second };
			float	line2[4] = { edge2_source.first, edge2_source.second, edge2_target.first, edge2_target.second };
			if (!(es[i].first == es[j].first || es[i].second == es[j].second || es[i].second == es[j].first || es[i].first == es[j].second)) {
				if (lineSegIntersect(line1, line2))
				{
					CrossingRemoveConst *CRC = new CrossingRemoveConst();
					CRC->setSolver(graph->solver);
					int const_id = graph->solver->extra_constraints.size();
					vector<int> degrees = { 0, 0, 0 };
					for (int m = 0; m < graph->n; m++) {
						if (findDistance(graph->distances, graph->n, m, es[i].first) == 1.0f ||
							findDistance(graph->distances, graph->n, m, es[i].second) == 1.0f)
							degrees[0]++;
						if (findDistance(graph->distances, graph->n, m, es[j].first == 1.0f) ||
							findDistance(graph->distances, graph->n, m, es[j].second) == 1.0f)
							degrees[1]++;
					}
					degrees[2] = degrees[0] + degrees[1];
					CRC->initShape(Edges{ es[i], es[j] }, vector<float> { 1, 1 }, weight, degrees, const_id);
					CRC->setOuterExtraConstId(graph->extra_outer_constraint_include_inner_nums.size());
					CRC0s.push_back(CRC);

					intersect_edges.push_back(es[i]);
					intersect_edges.push_back(es[j]);
				}
			}

		}
	}
	return CRC0s;

}