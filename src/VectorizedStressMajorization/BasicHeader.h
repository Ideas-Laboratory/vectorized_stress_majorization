#ifndef BasicHeader_H
#define BasicHeader_H

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <string>


using namespace std;


typedef std::pair<int, int> Edge;
typedef std::vector<std::pair<int, int> > Edges;
typedef std::vector<int> FaceList;
typedef std::vector<float> VertexList;
typedef std::vector<float> NormalList;
typedef std::vector<float> STLVectorf;
typedef std::vector<int> STLVectori;
typedef std::vector<std::vector<int> > AdjList;


typedef std::vector<std::string> VertexNameList;
typedef std::vector<int> VertexValueList;
typedef std::vector<std::pair<float, float> > Points;

typedef Eigen::SparseMatrix<float> SparseMatrix;
typedef Eigen::Triplet<float> Triplet;
typedef std::vector<Eigen::Triplet<float> > TripletList;
typedef Eigen::VectorXf VectorXf;
typedef Eigen::VectorXi VectorXi;
typedef Eigen::Vector3f Vector3f;
typedef Eigen::Matrix3f Matrix3f;
typedef Eigen::Vector2f Vector2f;
typedef Eigen::Vector4f Vector4f;
typedef Eigen::Matrix2f Matrix2f;
typedef Eigen::Matrix4f Matrix4f;
typedef Eigen::MatrixXf MatrixXf;
typedef Eigen::MatrixXi MatrixXi;
typedef Eigen::SimplicialCholesky<Eigen::SparseMatrix<float> > SimplicialCholesky;

typedef Eigen::MatrixXd MatrixXd;


struct Point_hull
{
	float x;
	float y;
};

enum DIM{ XDIM = 0, YDIM = 1 };

#define MAX_CHAR 128
#endif
