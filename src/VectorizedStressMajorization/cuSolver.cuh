#ifndef CU_SOLVER
#define CU_SOLVER
#include "stdafx.h"
#include <vector>
//#include "BasicHeader.h"
//using namespace Eigen;
//using namespace viennacl;
class CudaCore;
struct cuSparseMat{
	int *RowPtr, *ColIdx;
	float *Values;
	int const_size;
	cuSparseMat(int* RowPtr, int *ColIdx, float *Values, int const_size)
		: RowPtr(RowPtr), ColIdx(ColIdx), Values(Values), const_size(const_size)
	{}
	~cuSparseMat() {
		delete[]RowPtr;
		delete[]ColIdx;
		delete[]Values;
	}
};
class CudaSolver {//proxy
public:
	void init(float *, int, cuDWORD *, int*, int, float **);
	//void init(float *mat, int n_mat, cuSparseMat *J, cuDWORD* constraints, int* const_idx, std::vector<int> *closest_idx);
	void Solve(float *, float *, float *, float *, int, int*, float*);
	CudaSolver();
	~CudaSolver();

	void setMaxIter(int mi);

private:
	CudaCore *core;
};

#endif