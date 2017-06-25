#include "BasicHeader.h"
#include "StressConst.h" 
#include <fstream> 
#include <iostream> 
#include "HeapSort.h"

using namespace std;

//random nodes positions
int randomLayout(int num_of_points, VertexList& vertextlist_x, VertexList& vertextlist_y);//Points& points
float RandFloatRange();
void randomCoordinates(int node_num, VectorXf &xs, VectorXf &ys);

//rendering
Vector3f MovePoints(VectorXf &ps_x, VectorXf &ps_y);

//constraints tools
vector<float> generateLengthsOfCircle(int typesnum, float edge_length, float inter_angle, int P_Num);

//file reading tool
vector<vector<int>> communitiesFromFile(string communitiespath);
void readCoords(VectorXf& xs, VectorXf& ys, string filepath, int n);
void readCoordsAndWH(string filepath, Points &points, vector<float> &wh);
vector<string> split(const string &s, const string &seperator);
vector<string> split2(string str, string pattern);
MatrixXf readAdjFromFile(string adjMatrixfilePath, int num);
float* readAdjFromFileToShortestPath(string adjMatrixfilePath, int num);
void FromTxtEdgeFile(string filepath, Edges &edges, int &numOfPoints);
void readCirclesFileToCircles(string filepath, vector<vector<int>> &circles);
void readNameList(string filepath, vector<string>& namelist);

//community noverlap tools
pair<float, float> nearestToOrigin(vector<float> boundary);
vector<pair<float, float>> nearestToOrigins(vector<vector<float>> &boundaries);
vector<float>  convexHull(vector<Point_hull> points, int n);
vector<float> CalcConvexHull(vector<float> vertex_list);
vector<vector<float>> CalcConvexHulls(vector<vector<int>> communities, VectorXf vertex_list);
bool pointInHull(vector<float> A, pair<float, float> p);
float xmulti(pair<float, float> p1, pair<float, float> p2, pair<float, float> p0);
vector<Point_hull> MinkowskiDifferenceA_B(vector<float> setA, vector<float> setB);
vector<vector<float>> getMinusBoundaries(vector<vector<int>> communities, VectorXf vertex_list_x, VectorXf vertex_list_y);

//detect line intersect
bool lineSegIntersect(float line1[], float line2[]);

//math tools
float sqrt2(const float x);
float computeTheta(Vector2f edge);
void keepRelativeDire(float edge_x, float edge_y, pair<float, float>& dire);

//others
bool intInVector(int a, vector<int> A);
bool nodeInPolygon(pair<float, float> p, vector<float> area);
float findDistance(float* shortest_path, int node_num, int i, int j);

//symmetry constraint tools
void _rotate(float alpha, VectorXf& aixs);
void _rotate(float alpha, VectorXf &popt_x, VectorXf &popt_y);
void moveCenterToZero(float center_x, VectorXf& popt);
void moveCenterToZero(float center_x, VectorXf &popt_x, VectorXf &popt_y);
void mapAreaRightToLeft(VectorXf popt_x, VectorXf popt_y, vector<int> pids, vector<pair<int, int>> &closest_node_id_pairs);
void drawSymmetryAxis(VectorXf axis);
int orientation(Point_hull p, Point_hull q, Point_hull r);


/*Find the nodes with x-coordinates in xs and y-coordinates in ys
* which is in lasso_area
* the result indexes are in lasso_area_nodes
*return: true -- if there is node in the lasso area
*        false -- if there is no node in the lasso area*/
bool AddLassoNode(VectorXf xs, VectorXf ys, vector<float> lasso_area, vector<int>& lasso_area_nodes);