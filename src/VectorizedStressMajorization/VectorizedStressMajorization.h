#pragma once

#include <GL/glut.h> 

#include "ShortestPath.h"  
#include "Solver.cuh"
#include <fstream>

#include <iostream> 
#include "BasicHeader.h" 
#include "StressConst.h"
#include "CrossingRemoveConst.h" 
#include "PathRouteDijkstra.h"
#include <ctime> 
#include "CircleConst.h" 
#include "EqualAngleConst.h" 
#include "NonOverlapConst.h" 
#include "SymmetryConst.h" 
#include "Graph.h" 

static Graph *graph; 
static float scalecolors[] = {
	1.0, 0.918, 0.506,
	0.914, 0.780, 0.129,
	0.914, 0.686, 0.129,
	0.914, 0.58, 0.129,
	0.914, 0.455, 0.129,
	0.914, 0.318, 0.129,
	0.914, 0.176, 0.129
}; 
static vector<float> colors{
	241.0f / 255, 89.0f / 255, 89.0f / 255,//watermalon red 0
	86.0f / 255, 153.0f / 255, 102.0f / 255, //green 3
	23.0f / 255, 182.0f / 255, 222.0f / 255,//light blue 4
	60.0f / 255, 126.0f / 255, 200.0f / 255,//dark blue 5
	120.0f / 255, 109.0f / 255, 233.0f / 255,//purple 6
	232.0f / 255, 143.0f / 255, 206.0f / 255,//pink 7
	169.0f / 255, 32.0f / 255, 32.0f / 255,//blood red 8

	0.68627f, 0.57647f, 1.0f,
	0.57647f, 0.57647f, 1.0f,
	0.57647f, 0.68627f, 1.0f,
	0.61569f, 0.8549f, 1.0f,
	0.74902f, 1.0f, 0.94902f,
	0.74902f, 1.0f, 0.74902f,
	0.85098f, 1.0f, 0.74902f,
	1.0f, 0.843137f, 0.576470f,
	1.0f, 0.7647f, 0.59216f,
	1.0f, 0.7098f, 0.61176f,
	1.0f, 0.63137f, 0.58431f,
	1.0f, 0.57647f, 0.63137f,
	1.0f, 0.57647f, 1.0f,
	0.8f, 0.8f, 0.8f,
	0.78f, 0.94f, 0.94f,
	0.87f, 0.875f, 0.937f
};
void renderScene(void);
