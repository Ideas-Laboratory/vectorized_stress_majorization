#include "Tool.h" 

int randomLayout(int num_of_points, VertexList& vertextlist_x, VertexList& vertextlist_y) {//Points& points
	float x, y;
	for (int i = 0; i < num_of_points; i++) {
		x = RandFloatRange();
		y = RandFloatRange();
		vertextlist_x.push_back(x);
		vertextlist_y.push_back(y);
	}

	return 1;
}
float RandFloatRange() {
	float r = rand() / (float)(RAND_MAX); //[0,1]
	r = (r * 2 - 1.0f) / 4.0f;
	return r;
}

//rendering
Vector3f MovePoints(VectorXf &ps_x, VectorXf &ps_y) {
	int num = ps_x.size();
	float max_x = -9999999.0f, max_y = -9999999.0f, min_x = 9999999.0f, min_y = 9999999.0f;
	for (int i = 0; i < num; ++i) {
		if (max_x < ps_x[i]) {
			max_x = ps_x[i];
		}
		if (min_x - ps_x[i] > 0) {
			min_x = ps_x[i];
		}
		if (max_y < ps_y[i]) {
			max_y = ps_y[i];

		}
		if (min_y > ps_y[i]) {
			min_y = ps_y[i];
		}
	}
	float scale1 = 1;
	if (max_x - min_x > max_y - min_y) {
		scale1 = (max_x - min_x) / 2.0f;  //2 is the width and height of the dispalyer
	}
	else {
		scale1 = (max_y - min_y) / 2.0f;
	}
	max_x = max_x / scale1;
	min_x = min_x / scale1;
	max_y = max_y / scale1;
	min_y = min_y / scale1;

	float move1 = 0.0f, move2 = 0.0f;
	if (max_x - 1.0f > 0.0f) {
		move1 = 1.0f - max_x;
	}if (min_x + 1.0f < 0.0f) {
		move1 = -min_x - 1.0f;
	}if (max_y - 1.0f > 0.0f) {
		move2 = 1.0f - max_y;
	}if (min_y + 1.0f < 0.0f) {
		move2 = -min_y - 1.0f;
	}
	Vector3f movePoints = { move1, move2, scale1 };
	return movePoints;
}
//Constraints tools
float computeTheta(Vector2f edge){
	float angle;

	float productValue = edge[0];
	float edge_l = sqrt(edge[0] * edge[0] + edge[1] * edge[1]);
	float cosValue = productValue / (edge_l);
	angle = acos(cosValue);
	if (edge[1] <= 0){
		angle = 2 * M_PI - angle;
	}
	return angle;
}
vector<float> generateLengthsOfCircle(int typesnum, float edge_length, float inter_angle, int P_Num){
	vector<float> lengths;
	//first length is the edge length
	lengths.push_back(edge_length);

	//second length
	float l2 = 2 * edge_length * edge_length - 2 * edge_length*edge_length* cosf(inter_angle);//(M_PI - inter_angle) / 2
	lengths.push_back(sqrtf(l2));

	for (int i = 2; i < typesnum; i++){
		float edge_length2 = powf(edge_length, 2);


		float angle_in_cos = inter_angle - acosf((edge_length2 + lengths[i - 1] * lengths[i - 1] - lengths[i - 2] * lengths[i - 2]) / (2 * lengths[i - 1] * edge_length));
		float li = lengths[i - 1] * lengths[i - 1] + edge_length2 - 2 * lengths[i - 1] * edge_length * cosf(angle_in_cos);
		lengths.push_back(sqrtf(li));

	}

	int t = 0;
	if (P_Num % 2 == 0){
		t = typesnum - 2;
	}
	else{
		t = typesnum - 1;
	}
	for (int i = t; i >= 0; i--){
		lengths.push_back(lengths[i]);
	}

	return lengths;
}

//file reading tool
void FromTxtEdgeFile(string filepath, Edges &edges, int &numOfPoints){
	numOfPoints = 0;
	ifstream fin(filepath);
	string s = "";

	std::vector<std::string> ss;
	numOfPoints = 0;

	while (getline(fin, s))
	{
		ss = split2(s, " ");
		if (ss.size() == 2){
			int source = std::stoi(ss[0]);
			int target = std::stoi(ss[1]);
			if (source != target){
				edges.push_back(make_pair(target, source));
			}

			if (source > numOfPoints || target > numOfPoints){
				if (source > target){
					numOfPoints = source;
				}
				else{
					numOfPoints = target;
				}
			}
		}
	}
	numOfPoints = numOfPoints + 1;
}

vector<string> split2(string str, string pattern)
{
	string::size_type pos;
	vector<string> result;
	str += pattern;
	int size = str.size();

	for (int i = 0; i<size; i++)
	{
		pos = str.find(pattern, i);
		if (pos<size)
		{
			string s = str.substr(i, pos - i);
			result.push_back(s);
			i = pos + pattern.size() - 1;
		}
	}
	return result;
}
vector<string> split(const string &s, const string &seperator) {
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while (i != s.size()) {
		int flag = 0;
		while (i != s.size() && flag == 0) {
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x]) {
					++i;
					flag = 0;
					break;
				}
		}
		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0) {
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x]) {
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j) {
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}

MatrixXf readAdjFromFile(string adjMatrixfilePath, int num) {
	ifstream fin(adjMatrixfilePath);
	MatrixXf c(num, num);

	string s = "";
	int i = 0;
	while (getline(fin, s))
	{
		std::vector<std::string> ss = split(s, ",");
		for (int j = 0; j < num; j++) {
			c(i, j) = std::stoi(ss[j]);
		}
		i++;
	}
	fin.close();
	return c;
}
float* readAdjFromFileToShortestPath(string adjMatrixfilePath, int num) {
	float* sp = new float[num*(num - 1) / 2];
	ifstream fin(adjMatrixfilePath);

	string s = "";
	int i = 0;
	while (getline(fin, s))
	{
		std::vector<std::string> ss = split(s, ",");
		for (int j = i + 1; j < num; j++) {
			sp[(2 * (num - 1) - i + 1)*i / 2 + (j - i - 1)] = std::stof(ss[j]);
			if (std::stof(ss[j]) >= num) {
				sp[(2 * (num - 1) - i + 1)*i / 2 + (j - i - 1)] = 15.0f;
			}
		}
		i++;
	}
	fin.close();
	return sp;
}

void readCirclesFileToCircles(string filepath, vector<vector<int>> &circles){
	ifstream fin(filepath);
	string s = "";
	int index = 0;
	while (getline(fin, s))
	{
		vector<string> ss = split(s, ",");
		vector<int> circle;
		for (int i = 0; i < ss.size(); i++){
			circle.push_back(atoi(ss[i].c_str()));
		}

		circles.push_back(circle);
		index++;
	}
}
vector<vector<int>> communitiesFromFile(string communitiespath){
	ifstream fin(communitiespath);
	int com_num = 0;
	int com_index = 0;
	vector<vector<int>> communities;

	string s = "";
	while (getline(fin, s))
	{
		std::vector<std::string> ss = split(s, " ");
		if (ss.size() == 3 && ss[1] == "Communities:"){
			com_num = atoi(ss[2].c_str());
			communities.resize(com_num);
		}
		ss = split(s, "\t");
		if (ss.size() == 2){
			com_index = atoi(ss[1].c_str());
			communities[com_index].push_back(atoi(ss[0].c_str()));
		}

	}
	fin.close();
	return communities;
}
void readNameList(string filepath, vector<string>& namelist) {
	ifstream fin(filepath);
	string s = "";
	while (getline(fin, s))
	{
		if (s != "") namelist.push_back(s);
	}
	fin.close();
}


pair<float, float> nearestToOrigin(vector<float> boundary){
	if (boundary.size() == 0){
		return make_pair(1, 1);
	}
	pair<float, float> result = make_pair(0, 0);
	float length = boundary[2 * 0 + 0] * boundary[2 * 0 + 0] + boundary[2 * 0 + 1] * boundary[2 * 0 + 1];
	result.first = boundary[2 * 0 + 0];
	result.second = boundary[2 * 0 + 1];

	for (int i = 1; i < boundary.size() / 2; i++){
		float length_t = boundary[2 * i + 0] * boundary[2 * i + 0] + boundary[2 * i + 1] * boundary[2 * i + 1];
		if (length > length_t){
			length = length_t;
			result.first = boundary[2 * i + 0];
			result.second = boundary[2 * i + 1];
		}
	}
	return result;
}
vector<pair<float, float>> nearestToOrigins(vector<vector<float>>& boundaries){
	vector<pair<float, float>> results;
	for (int i = 0; i < boundaries.size(); i++){
		vector<float> angles;
		vector<int> boun_node_ids;
		for (int j = 0; j < boundaries[i].size() / 2; j++){
			float angle = computeTheta(Vector2f(boundaries[i][2 * j + 0], boundaries[i][2 * j + 1]));
			angles.push_back(angle);
			boun_node_ids.push_back(j);
		}
		HeapSort(angles, angles.size(), boun_node_ids);
		vector<float> boundary = boundaries[i];
		for (int j = 0; j < boun_node_ids.size(); j++){
			boundary[2 * j + 0] = boundaries[i][2 * boun_node_ids[j] + 0];
			boundary[2 * j + 1] = boundaries[i][2 * boun_node_ids[j] + 1];
		}
		boundaries[i] = boundary;

		if (nodeInPolygon(make_pair(0.0f, 0.0f), boundary))
			results.push_back(nearestToOrigin(boundary));
		else
			results.push_back(make_pair(0, 0));
	}
	return results;

}
struct POINT1 { double x; double y; };

typedef vector<POINT1> PTARRAY;

bool operator==(const POINT1 &pt1, const POINT1 &pt2) {
	return (pt1.x == pt2.x && pt1.y == pt2.y);
}
bool CompareVector(const POINT1 &pt1, const POINT1 &pt2) {

	double m1 = sqrt((double)(pt1.x * pt1.x + pt1.y * pt1.y));
	double m2 = sqrt((double)(pt2.x * pt2.x + pt2.y * pt2.y));

	double v1 = pt1.x / m1, v2 = pt2.x / m2;
	return (v1 > v2 || (v1 == v2 && m1 < m2));
}

vector<float> CalcConvexHull(vector<float> vertex_list) {

	PTARRAY vecSrc;
	for (int i = 0; i < vertex_list.size(); i = i + 2){
		POINT1 p = { vertex_list[i], vertex_list[i + 1] };
		vecSrc.push_back(p);
	}

	if (vecSrc.size() < 3) {
		return vertex_list;
	}

	POINT1 ptBase = vecSrc.front();
	for (PTARRAY::iterator i = vecSrc.begin() + 1; i != vecSrc.end(); ++i) {

		if (i->y < ptBase.y || (i->y == ptBase.y && i->x > ptBase.x)) {

			ptBase = *i;
		}
	}

	for (PTARRAY::iterator i = vecSrc.begin(); i != vecSrc.end();) {

		if (*i == ptBase) {
			i = vecSrc.erase(i);
		}
		else {

			i->x -= ptBase.x, i->y -= ptBase.y;
			++i;
		}
	}

	sort(vecSrc.begin(), vecSrc.end(), &CompareVector);

	vecSrc.erase(unique(vecSrc.begin(), vecSrc.end()), vecSrc.end());

	for (PTARRAY::reverse_iterator ri = vecSrc.rbegin();
		ri != vecSrc.rend() - 1; ++ri) {
		PTARRAY::reverse_iterator riNext = ri + 1;

		ri->x -= riNext->x, ri->y -= riNext->y;
	}

	for (PTARRAY::iterator i = vecSrc.begin() + 1; i != vecSrc.end(); ++i) {

		for (PTARRAY::iterator iLast = i - 1; iLast != vecSrc.begin();) {
			int v1 = i->x * iLast->y, v2 = i->y * iLast->x;

			if (v1 < v2 || (v1 == v2 && i->x * iLast->x > 0 &&
				i->y * iLast->y > 0)) {
				break;
			}

			i->x += iLast->x, i->y += iLast->y;
			iLast = (i = vecSrc.erase(iLast)) - 1;
		}
	}

	vecSrc.front().x += ptBase.x, vecSrc.front().y += ptBase.y;
	for (PTARRAY::iterator i = vecSrc.begin() + 1; i != vecSrc.end(); ++i) {
		i->x += (i - 1)->x, i->y += (i - 1)->y;
	}
	vecSrc.push_back(ptBase);
	vector<float> boundary_vertex;
	for (int i = 0; i < vecSrc.size(); i++){
		boundary_vertex.push_back(vecSrc[i].x);
		boundary_vertex.push_back(vecSrc[i].y);
	}
	return boundary_vertex;
}
vector<vector<float>> CalcConvexHulls(vector<vector<int>> communities, VectorXf vertex_list){
	vector<vector<float>> boundaries;
	for (int i = 0; i < communities.size(); i++){
		vector<Point_hull> point_community;
		for (int j = 0; j < communities[i].size(); j++){
			Point_hull p;
			p.x = vertex_list[2 * communities[i][j] + 0];
			p.y = vertex_list[2 * communities[i][j] + 1];

			point_community.push_back(p);
		}
		vector<float> boundary = convexHull(point_community, point_community.size());
		boundaries.push_back(boundary);
	}

	return boundaries;
}
bool lineSegIntersect(float line1[], float line2[]){
	float a = line1[2] - line1[0];
	float b = line2[0] - line2[2];
	float c = line1[3] - line1[1];
	float d = line2[1] - line2[3];
	float g = line2[0] - line1[0];
	float h = line2[1] - line1[1];
	float f = a * d - b * c;
	if (fabs(f) < 1.0e-06){
		return false;
	}

	float t = (d*g - b*h) / f;
	float s = (-c*g + a*h) / f;
	if ((0 > t) || (t > 1)){
		return false;
	}

	if ((0 > s) || (s > 1)){
		return false;
	}

	return true;
}


#define SQRT_MAGIC_F 0x5f3759df 
float  sqrt2(const float x)
{
	const float xhalf = 0.5f*x;

	union // get bits for floating value
	{
		float x;
		int i;
	} u;
	u.x = x;
	u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
	return x*u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy 
}


float xmulti(pair<float, float> p1, pair<float, float> p2, pair<float, float> p0)  
{
	return (p1.first - p0.first)*(p2.second - p0.second) - (p2.first - p0.first)*(p1.second - p0.second);
}

#define eps 1e-10
bool pointInHull(vector<float> A, pair<float, float> p){

	int pnum = A.size() / 2;
	//judging if the node is in the polygon by dichotomy. Not on the margin
	//for (i = 1; i <= m; i++){    //B[i]
	if (xmulti(p, make_pair(A[2], A[3]), make_pair(A[0], A[1])) <= eps || xmulti(p, make_pair(A[2 * (pnum - 1)], A[2 * (pnum - 1) + 1]), make_pair(A[0], A[1])) >= -eps)    //在第一个点为起点的扇形之外或在边上
		return false;
	int left = 1, right = pnum - 1;
	while (right - left != 1){
		int mid = (left + right) / 2;
		if (xmulti(p, make_pair(A[2 * mid], A[2 * mid + 1]), make_pair(A[0], A[1]))>eps)
			left = mid;
		else
			right = mid;
	}
	if (xmulti(p, make_pair(A[2 * right], A[2 * right + 1]), make_pair(A[2 * left], A[2 * left + 1])) <= eps)    //在边之外或在边上
		return false;
	//}
	return true;
}
vector<Point_hull> MinkowskiDifferenceA_B(vector<float> setA, vector<float> setB){
	//B is a boundary, A is a set of points ,A - B's boundary = (A-B)'s boundary
	vector<Point_hull> C;
	for (int i = 0; i < setA.size() / 2; i++){
		for (int j = 0; j < setB.size() / 2; j++){
			Point_hull p;
			p.x = setA[2 * i] - setB[2 * j];
			p.y = setA[2 * i + 1] - setB[2 * j + 1];
			C.push_back(p);
		}
	}
	return C;
}

vector<vector<float>> getMinusBoundaries(vector<vector<int>> communities, VectorXf vertex_list_x, VectorXf vertex_list_y){
	vector<vector<float>> minusBoundaries;
	vector<vector<float>> vertex_list_communities;
	for (int i = 0; i < communities.size(); i++){
		vector<float> vertex_list_community;
		for (int j = 0; j < communities[i].size(); j++){
			vertex_list_community.push_back(vertex_list_x[communities[i][j]]);
			vertex_list_community.push_back(vertex_list_y[communities[i][j]]);
		}
		vertex_list_communities.push_back(vertex_list_community);
	}

	for (int i = 0; i < communities.size(); i++){
		for (int j = i + 1; j < communities.size(); j++){
			vector<Point_hull> minusSet;
			minusSet = MinkowskiDifferenceA_B(vertex_list_communities[i], vertex_list_communities[j]);
			vector<float> boundary;
			boundary = convexHull(minusSet, minusSet.size());
			minusBoundaries.push_back(boundary);
		}
	}

	return minusBoundaries;
}

bool intInVector(int a, vector<int> A){
	for (int i = 0; i < A.size(); i++){
		if (a == A[i]){
			return true;
		}
	}
	return false;
}

bool nodeInPolygon(pair<float, float> p, vector<float> area) {
	int count = area.size() / 2;
	if (count < 3) {
		return false;
	}
	bool result = false;
	for (int i = 0, j = count - 1; i < count; i++) {
		if (area[2 * i + 0] < p.first && area[2 * j + 0] >= p.first
			|| area[2 * j + 0] < p.first && area[2 * i + 0] >= p.first) {

			if (area[2 * i + 1] + (p.first - area[2 * i + 0]) / (area[2 * j + 0] - area[2 * i + 0])*(area[2 * j + 1] - area[2 * i + 1]) < p.second) {
				result = !result;
			}
		}
		j = i;
	}
	return result;
}
void _rotate(float alpha, VectorXf& popt) {
	for (int i = 0; i < popt.size()/2; i++) {
		float thow = (sqrtf(popt[2 * i + 0] * popt[2 * i + 0] + popt[2 * i + 1] * popt[2 * i + 1]));
		float p_angle = 0; // x==0, y==0
		if (popt[2 * i + 0] > 0) {
			p_angle = atanf(popt[2 * i + 1] / popt[2 * i + 0]);
		}
		else if (popt[2 * i + 0] < 0 && popt[2 * i + 1] >= 0) {
			p_angle = atanf(popt[2 * i + 1] / popt[2 * i + 0]) + M_PI;
		}
		else if (popt[2 * i + 0] < 0 && popt[2 * i + 1] < 0) {
			p_angle = atanf(popt[2 * i + 1] / popt[2 * i + 0]) - M_PI;
		}
		else if (popt[2 * i + 0] == 0 && popt[2 * i + 1] < 0) {
			p_angle = -M_PI / 2;
		}
		else if (popt[2 * i + 0] == 0 && popt[2 * i + 1] > 0) {
			p_angle = M_PI / 2;
		}
		popt[2 * i + 0] = thow * cos(p_angle - alpha);
		popt[2 * i + 1] = thow * sin(p_angle - alpha);
	}
}
void _rotate(float alpha, VectorXf& popt_x, VectorXf &popt_y) {
	for (int i = 0; i < popt_x.size(); i++) {
		float thow = (sqrtf(popt_x[i] * popt_x[i] + popt_y[i] * popt_y[i]));
		float p_angle = 0; // x==0, y==0
		if (popt_x[i] > 0) {
			p_angle = atanf(popt_y[i] / popt_x[i]);
		}
		else if (popt_x[i] < 0 && popt_y[i] >= 0) {
			p_angle = atanf(popt_y[i] / popt_x[i]) + M_PI;
		}
		else if (popt_x[i] < 0 && popt_y[i] < 0) {
			p_angle = atanf(popt_y[i] / popt_x[i]) - M_PI;
		}
		else if (popt_x[i] == 0 && popt_y[i] < 0) {
			p_angle = -M_PI / 2;
		}
		else if (popt_x[i] == 0 && popt_y[i] > 0) {
			p_angle = M_PI / 2;
		}
		popt_x[i] = thow * cos(p_angle - alpha);
		popt_y[i] = thow * sin(p_angle - alpha);
	}
}
void moveCenterToZero(float center_x, VectorXf& popt) {
	for (int i = 0; i < popt.size()/2; i++) {
		popt[i] -= center_x;
	}
}
void moveCenterToZero(float center_x, VectorXf &popt_x, VectorXf &popt_y) {
	for (int i = 0; i < popt_x.size(); i++) {
		popt_x[i] -= center_x;
	}
}
/*greedy mapping*/
void mapAreaRightToLeft(VectorXf popt_x, VectorXf popt_y, vector<int> pids, vector<pair<int, int>> &closest_node_id_pairs) {
	closest_node_id_pairs.clear();
	vector<int> left_node_ids, right_node_ids;
	for (int i = 0; i < pids.size(); i++) {
		if (popt_x[pids[i]] < 0) {
			left_node_ids.push_back(pids[i]);
		}
		else if (popt_x[pids[i]] > 0) {
			right_node_ids.push_back(pids[i]);
		}
	}
	while (!(left_node_ids.size() == 0 || right_node_ids.size() == 0)) {
		float shortest_dis_sq = 9999999999;
		pair<int, int> current_closest_pair;
		int closest_left_node_idx, closest_right_node_idx;
		for (int i = 0; i < left_node_ids.size(); i++) {
			for (int j = 0; j < right_node_ids.size(); j++) {
				float x_dis = popt_x[left_node_ids[i]] - (-popt_x[right_node_ids[j]]);
				float y_dis = popt_y[left_node_ids[i]] - (popt_y[right_node_ids[j]]);
				float c_dis = x_dis*x_dis + y_dis*y_dis;
				if (shortest_dis_sq > c_dis) {
					shortest_dis_sq = c_dis;
					current_closest_pair = make_pair(left_node_ids[i], right_node_ids[j]);
					closest_left_node_idx = i;
					closest_right_node_idx = j;
				}
			}
		}
		closest_node_id_pairs.push_back(current_closest_pair);
		left_node_ids.erase(left_node_ids.begin() + closest_left_node_idx);
		right_node_ids.erase(right_node_ids.begin() + closest_right_node_idx);
	}
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point_hull p, Point_hull q, Point_hull r)
{
	float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

	if (val == 0)
		return 0; // colinear
	return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// Finds convex hull of a set of n points.
vector<float> convexHull(vector<Point_hull> points, int n)
{
	vector<float> output_boundary;
	// There must be at least 3 points
	if (n < 3){
		for (int i = 0; i < points.size(); i++){
			output_boundary.push_back(points[i].x);
			output_boundary.push_back(points[i].y);
		}
		return output_boundary;
	}

	// Initialize Result
	//int next[n];
	VectorXi next = VectorXi::Zero(n);
	for (int i = 0; i < n; i++)
		next[i] = -1;

	// Find the leftmost point
	int l = 0;
	for (int i = 1; i < n; i++)
		if (points[i].x < points[l].x)
			l = i;

	// Start from leftmost point, keep moving counterclockwise
	// until reach the start point again
	int p = l, q;
	do
	{
		// Search for a point 'q' such that orientation(p, i, q) is
		// counterclockwise for all points 'i'
		q = (p + 1) % n;
		for (int i = 0; i < n; i++)
			if (orientation(points[p], points[i], points[q]) == 2)
				q = i;

		next[p] = q; // Add q to result as a next point of p
		p = q; // Set p as q for next iteration
	} while (p != l);

	// Result
	for (int i = 0; i < n; i++)
	{
		if (next[i] != -1){
			output_boundary.push_back(points[i].x);
			output_boundary.push_back(points[i].y);
		}
	}
	return output_boundary;
}
bool AddLassoNode(VectorXf xs, VectorXf ys, vector<float> lasso_area, vector<int>& lasso_area_nodes){
	for (int i = 0; i < xs.size(); i++) {
		pair<float, float> p = make_pair(xs[i], ys[i]);
		if (nodeInPolygon(p, lasso_area)) {
			lasso_area_nodes.push_back(i);
		}
	}
	if (lasso_area_nodes.size() == 0)
		return false;
	return true;
}
void randomCoordinates(int node_num, VectorXf &xs, VectorXf &ys){
	VectorXf result_xs(node_num);
	VectorXf result_ys(node_num);

	float x, y;
	cout << RandFloatRange() << "," << RandFloatRange() << endl;
	for (int i = 0; i < node_num; i++) {
		x = RandFloatRange();
		y = RandFloatRange();
		result_xs[i] = x;
		result_ys[i] = y;
	}
	xs = result_xs;
	ys = result_ys;
}
float findDistance(float* shortest_path, int node_num, int i, int j){
	if (i < j)
		return shortest_path[(2 * (node_num - 1) - i + 1)*i / 2 + (j - i - 1)];
	else
		return shortest_path[(2 * (node_num - 1) - j + 1)*j / 2 + (i - j - 1)];
}
void keepRelativeDire(float edge_x, float edge_y, pair<float, float>& dire){
	if (edge_x * dire.first + edge_y * dire.second < 0){
		dire.first = -dire.first;
		dire.second = -dire.second;
	}
}