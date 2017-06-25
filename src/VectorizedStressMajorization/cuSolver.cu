#include "cuSolver.cuh"
#include "stdafx.h"

#ifndef VIENNACL_WITH_CUDA
#define VIENNACL_WITH_CUDA
#endif

#ifdef min
#undef min
#endif  
#ifdef __INTELLISENSE__
#define __CUDACC__
#endif

#include <iostream>
#include <vector>
#include "viennacl/linalg/cg.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/matrix_proxy.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"

#ifdef __INTELLISENSE__
#undef __global__
#undef __device__
#undef __constant__
#undef __forceinline__
#undef __shared__
#undef __restrict__

#define __global__
#define __device__
#define __constant__
#define __forceinline__
#define __shared__
#define __restrict__
#endif

#define dist(i ,j) ((i)<(j)?dist[((j) - ((i)+1) + ((i)*(2*n-(i)-1))/2)]:dist[((i) - ((j)+1) + ((j)*(2*n-(j)-1))/2)])


/*Rotate and move for Symmetry Constraints*/
__inline__ __device__ void _rotate(float alpha, float &edge_x, float &edge_y);
__inline__ __device__ void moveCenterToZero(float center_x,
	float &edge_x, float &edge_y);
/*Shape matching
* Rotate and transformate the relative_shape
* to map the corresponding nodes in P_Opt*/
__inline__ __device__ void ShapeMatching(float *relative_shape,
	float *P_Opt_x, float *P_Opt_y, int *pid, int node_num);

/*For Equal Angle Constraints*/
__forceinline__ __device__ int EqualAngleKernel(
	void* __restrict__ constraints,
	float* __restrict__ P_Opt_x,
	float* __restrict__ P_Opt_y,
	float* __restrict__ _right_hand_x,
	float* __restrict__ _right_hand_y) {

	int*   i_constraints = (int  *)constraints;
	float* f_constraints = (float *)constraints;
	int pid_size = i_constraints[1], pid_begin = 5;
	float p_length = f_constraints[3];
	float param = f_constraints[4];
	int relative_size = 2 * pid_size, relative_begin = pid_begin + pid_size;

	int warpid = threadIdx.x >> 5;
	int laneid = threadIdx.x - warpid << 5;
	__shared__ int voting[4];
	if (laneid == 0)
		voting[warpid] = i_constraints[2];
	__syncthreads();
	if (voting[warpid] != i_constraints[2])
		voting[warpid] = -1;
	__syncthreads();

	ShapeMatching(&f_constraints[relative_begin], P_Opt_x, P_Opt_y, &i_constraints[pid_begin], pid_size);

	for (int i = 0; i < pid_size - 1; i++) {
		float dx = param * (f_constraints[relative_begin + relative_size - 2] - f_constraints[relative_begin + 2 * i + 0]);
		float dy = param * (f_constraints[relative_begin + relative_size - 1] - f_constraints[relative_begin + 2 * i + 1]);

		float sum = dx, sum2 = dy;
		if (voting[warpid] > 0)
		{
			for (int i = 16; i >= 1; i >>= 1)
			{
				sum += __shfl_down(sum, i);
				sum2 += __shfl_down(sum2, i);
			}
			if (laneid == 0)
			{
				//child
				atomicAdd(&_right_hand_x[i_constraints[pid_begin + i]], -sum);
				atomicAdd(&_right_hand_y[i_constraints[pid_begin + i]], -sum2);
			}
		}
		else {
			//child
			atomicAdd(&_right_hand_x[i_constraints[pid_begin + i]], -dx);
			atomicAdd(&_right_hand_y[i_constraints[pid_begin + i]], -dy);
		}
		//parent
		atomicAdd(&_right_hand_x[i_constraints[pid_begin + pid_size - 1]], dx);
		atomicAdd(&_right_hand_y[i_constraints[pid_begin + pid_size - 1]], dy);
	}

}
/*For Circle Constraints*/
__forceinline__ __device__ int CircleKernel(
	void* __restrict__ constraints,
	float* __restrict__ P_Opt_x,
	float* __restrict__ P_Opt_y,
	float* __restrict__ _right_hand_x,
	float* __restrict__ _right_hand_y) {

	int*   i_constraints = (int  *)constraints;
	float* f_constraints = (float *)constraints;
	float differ_angle = f_constraints[0];
	int pid_size = i_constraints[1], pid_begin = 4;
	int r_lengths_size = i_constraints[2], r_lengths_begin = pid_begin + pid_size;
	int relative_size = 2 * pid_size, relative_begin = r_lengths_begin + r_lengths_size;
	float param = f_constraints[3];

	int warpid = threadIdx.x >> 5;
	int laneid = threadIdx.x - warpid << 5;
	__shared__ int voting[4];
	if (laneid == 0)
		voting[warpid] = i_constraints[2];
	__syncthreads();
	if (voting[warpid] != i_constraints[2])
		voting[warpid] = -1;
	__syncthreads();

	ShapeMatching(&f_constraints[relative_begin], P_Opt_x, P_Opt_y, &i_constraints[pid_begin], pid_size);

	for (int i = 0; i < pid_size; i++) {

		float dx = f_constraints[relative_begin + 2 * ((i + 1) % pid_size) + 0] - f_constraints[relative_begin + 2 * i + 0];
		float dy = f_constraints[relative_begin + 2 * ((i + 1) % pid_size) + 1] - f_constraints[relative_begin + 2 * i + 1];

		//printf("in loop with %d nodes, d=(%f,%f),i=%d,\n ", pid_size, dx, dy, i);

		dx = dx*param / f_constraints[r_lengths_begin + i];
		dy = dy*param / f_constraints[r_lengths_begin + i];

		float sum = dx, sum2 = dy;
		if (voting[warpid] > 0)
		{
			for (int i = 16; i >= 1; i >>= 1)
			{
				sum += __shfl_down(sum, i);
				sum2 += __shfl_down(sum2, i);
			}
			if (laneid == 0)
			{
				//child
				atomicAdd(&_right_hand_x[i_constraints[pid_begin + i]], -sum);
				atomicAdd(&_right_hand_y[i_constraints[pid_begin + i]], -sum2);
			}
		}
		else {
			//child
			atomicAdd(&_right_hand_x[i_constraints[pid_begin + i]], -dx);
			atomicAdd(&_right_hand_y[i_constraints[pid_begin + i]], -dy);
		}
		//parent
		atomicAdd(&_right_hand_x[i_constraints[pid_begin + (i + 1) % pid_size]], dx);
		atomicAdd(&_right_hand_y[i_constraints[pid_begin + (i + 1) % pid_size]], dy);

	}
}
/*For Edge Crossing Remove Constraints*/
__forceinline__ __device__ int CrossingRemovalKernel(
	void* __restrict__ constraints,
	float* __restrict__ P_Opt_x,
	float* __restrict__ P_Opt_y,
	float* __restrict__ _right_hand_x,
	float* __restrict__ _right_hand_y) {

	float *f_constraints = (float *)constraints;
	int*   i_constraints = (int  *)constraints;

	int edge_size = i_constraints[0], edge_begin = 2;
	int rest_length_size = i_constraints[0], rest_length_begin = edge_begin + 2 * edge_size;
	int edge_weights_begin = rest_length_begin + rest_length_size;
	float param = f_constraints[1];

	/*the idea direction is the sum direction of intersecting edges*/
	float  dire_x = 0, dire_y = 0;
	for (int i = 0; i < edge_size; i++) {
		float edge_x = P_Opt_x[i_constraints[edge_begin + i * 2 + 0]] - P_Opt_x[i_constraints[edge_begin + i * 2 + 1]];
		float edge_y = P_Opt_y[i_constraints[edge_begin + i * 2 + 0]] - P_Opt_y[i_constraints[edge_begin + i * 2 + 1]];
		/*keep the relative direction roughly*/
		if (dire_x*edge_x + dire_y*edge_y < 0) {
			dire_x += -edge_x;
			dire_y += -edge_y;
		}
		else {
			dire_x += edge_x;
			dire_y += edge_y;
		}
	}

	float norm = sqrt(dire_x * dire_x + dire_y * dire_y);
	if (isnan(norm)) norm = 1;
	dire_x = dire_x * edge_size * param / norm;
	dire_y = dire_y * edge_size * param / norm;
	for (int i = 0; i < edge_size; i++) {
		float edge_x = P_Opt_x[i_constraints[edge_begin + i * 2 + 0]] - P_Opt_x[i_constraints[edge_begin + i * 2 + 1]];
		float edge_y = P_Opt_y[i_constraints[edge_begin + i * 2 + 0]] - P_Opt_y[i_constraints[edge_begin + i * 2 + 1]];
		float dx = dire_x * f_constraints[edge_weights_begin + i] / f_constraints[rest_length_begin + i];
		float dy = dire_y * f_constraints[edge_weights_begin + i] / f_constraints[rest_length_begin + i];
		/*keep the relative direction roughly*/
		if (dire_x*edge_x + dire_y*edge_y < 0) {
			//source
			atomicAdd(&_right_hand_x[i_constraints[edge_begin + i * 2 + 0]], -dx);
			atomicAdd(&_right_hand_y[i_constraints[edge_begin + i * 2 + 0]], -dy);
			//target
			atomicAdd(&_right_hand_x[i_constraints[edge_begin + i * 2 + 1]], dx);
			atomicAdd(&_right_hand_y[i_constraints[edge_begin + i * 2 + 1]], dy);
		}
		else {
			//source
			atomicAdd(&_right_hand_x[i_constraints[edge_begin + i * 2 + 0]], dx);
			atomicAdd(&_right_hand_y[i_constraints[edge_begin + i * 2 + 0]], dy);
			//target
			atomicAdd(&_right_hand_x[i_constraints[edge_begin + i * 2 + 1]], -dx);
			atomicAdd(&_right_hand_y[i_constraints[edge_begin + i * 2 + 1]], -dy);
		}
	}
	return edge_size * 4 + 2;

}
/*For Node Noverlap Constraints*/
__forceinline__ __device__ int NoverlapKernel(
	void* __restrict__ constraints,
	float* __restrict__ P_Opt_x,
	float* __restrict__ P_Opt_y,
	float* __restrict__ _right_hand_x,
	float* __restrict__ _right_hand_y) {
	float *f_constraints = (float *)constraints;
	int   *i_constraints = (int   *)constraints;

	float x_dis = P_Opt_x[i_constraints[2]] - P_Opt_x[i_constraints[3]];
	float y_dis = P_Opt_y[i_constraints[2]] - P_Opt_y[i_constraints[3]];
	float width = f_constraints[0];
	float height = f_constraints[1];
	int source_node = i_constraints[2], target_node = i_constraints[3];

	float eucli_dis = sqrt(x_dis*x_dis + y_dis*y_dis);
	if (isnan(eucli_dis))
	{
		eucli_dis = 1;
	}
	float param = f_constraints[4] / eucli_dis;
	/*these two nodes are overlap*/
	if (abs(x_dis) < width && abs(y_dis) < height) {
		float dx = 0, dy = 0;
		/*move along x-axis costing less, more desiring move along y-axis*/
		if (2 * (width - abs(x_dis)) < height - abs(y_dis)) {
			dx = f_constraints[4] * width;
			dy = param * abs(y_dis);
			/*keep the relative direction roughly*/
			if (dx*x_dis + dy*y_dis < 0) {
				dx *= -1;
				dy *= -1;
			}
			atomicAdd(&_right_hand_x[source_node], -dx);
			atomicAdd(&_right_hand_y[source_node], -dy);
			atomicAdd(&_right_hand_x[target_node], dx);
			atomicAdd(&_right_hand_y[target_node], dy);
			//right_hand[0] = param * width * 35;
		}
		else {/*move along y-axis costing less*/
			dx = param * abs(x_dis);
			dy = f_constraints[4] * height;
			/*keep the relative direction roughly*/
			if (dx*x_dis + dy*y_dis < 0) {
				dx *= -1;
				dy *= -1;
			}
			atomicAdd(&_right_hand_x[source_node], -dx);
			atomicAdd(&_right_hand_y[source_node], -dy);
			atomicAdd(&_right_hand_x[target_node], dx);
			atomicAdd(&_right_hand_y[target_node], dy);
			//right_hand[1] = height * 20;
		}
	}
	else {
		atomicAdd(&_right_hand_x[source_node], param*x_dis / f_constraints[5]);//* x_dis
		atomicAdd(&_right_hand_y[source_node], param*y_dis / f_constraints[5]);// * y_dis
		atomicAdd(&_right_hand_x[target_node], -param*x_dis / f_constraints[5]);//* x_dis
		atomicAdd(&_right_hand_y[target_node], -param*y_dis / f_constraints[5]);// *y_dis
	}

	return 5;
}
/*For Symmetry Constraints*/
__forceinline__ __device__ int SymmetryKernel(
	void* __restrict__ constraints,
	float* __restrict__ P_Opt_x,
	float* __restrict__ P_Opt_y,
	float* __restrict__ _right_hand_x,
	float* __restrict__ _right_hand_y) {

	float *f_constraints = (float *)constraints;
	int   *i_constraints = (int   *)constraints;
	/*node number in the relative shape = number of pairs of closest nodes*/
	int pid_size = i_constraints[2];
	int model_node_begin = 5;
	int other_node_begin = model_node_begin + pid_size;
	int r_length_begin = other_node_begin + pid_size, r_length_size = i_constraints[3];
	float param = f_constraints[4];
	float rotate_angle = f_constraints[0];
	float move_dis = f_constraints[1];

	int index = 0;
	for (int i = 0; i < pid_size; i++) {
		for (int j = i + 1; j < pid_size; j++) {
			/*direction of model edge(model node i to model node j)*/
			float dire_x = P_Opt_x[i_constraints[model_node_begin + i]] - P_Opt_x[i_constraints[model_node_begin + j]];
			float dire_y = P_Opt_y[i_constraints[model_node_begin + i]] - P_Opt_y[i_constraints[model_node_begin + j]];

			float eucli_dist_model = sqrtf(dire_x*dire_x + dire_y*dire_y);
			if (isnan(eucli_dist_model)) eucli_dist_model = 1;
			_rotate(rotate_angle, dire_x, dire_y);
			moveCenterToZero(move_dis, dire_x, dire_y);
			dire_x *= -1;
			//dire_x = dire_x / eucli_dist_model;
			//dire_y = dire_y / eucli_dist_model;

			_rotate(-rotate_angle, dire_x, dire_y);
			moveCenterToZero(-move_dis, dire_x, dire_y);

			float dx = param * dire_x / f_constraints[r_length_begin + index];
			float dy = param * dire_y / f_constraints[r_length_begin + index];

			//source
			atomicAdd(&_right_hand_x[i_constraints[other_node_begin + i]], dx);
			atomicAdd(&_right_hand_y[i_constraints[other_node_begin + i]], dy);
			//target
			atomicAdd(&_right_hand_x[i_constraints[other_node_begin + j]], -dx);
			atomicAdd(&_right_hand_y[i_constraints[other_node_begin + j]], -dy);

			index++;
		}
	}

	return 5 + 2 * i_constraints[2] + i_constraints[3];
}
/*For Atom Constraints*/
__forceinline__ __device__ int AtomKernel(
	void * __restrict__ constraints,
	float * __restrict__ P_Opt_x,
	float  * __restrict__ P_Opt_y,
	float  * __restrict__ _right_hand_x,
	float   * __restrict__ _right_hand_y) {//
	float *f_constraints = (float *)constraints;
	int   *i_constraints = (int   *)constraints;
	int warpid = threadIdx.x >> 5;
	int laneid = threadIdx.x - warpid << 5;

	//sum in every threads
	__shared__ int voting[4];
	if (laneid == 0)
		voting[warpid] = i_constraints[2];
	__syncthreads();
	if (voting[warpid] != i_constraints[2])
		voting[warpid] = -1;
	__syncthreads();

	//if this constraint does not assign te direction for this edge
	if (f_constraints[0] == 0 && f_constraints[1] == 0) {

		//original euclidean edge vector
		float dx = 0, dy = 0;
		dx = P_Opt_x[i_constraints[2]] - P_Opt_x[i_constraints[3]];//2 3
		dy = P_Opt_y[i_constraints[2]] - P_Opt_y[i_constraints[3]];
		if (dx == 0 && dy == 0) {
			dx = 0.0001f;
			dy = 0.0001f;
		}
		//weight para / normalization
		float para = f_constraints[4] / sqrtf(dx*dx + dy*dy);//4

															 //normalize and weight the edge vector
		dx *= para;
		dy *= para;

		float sum = dx, sum2 = dy;
		if (voting[warpid] > 0)
		{
			for (int i = 16; i >= 1; i >>= 1)
			{
				sum += __shfl_down(sum, i);//left shift variables in threads
				sum2 += __shfl_down(sum2, i);
			}
			if (laneid == 0)
			{
				atomicAdd(&_right_hand_x[i_constraints[2]], sum);
				atomicAdd(&_right_hand_y[i_constraints[2]], sum2);
			}
		}
		else {
			//edge's source
			atomicAdd(&_right_hand_x[i_constraints[2]], dx);
			atomicAdd(&_right_hand_y[i_constraints[2]], dy);
		}
		//edge's target
		atomicAdd(&_right_hand_x[i_constraints[3]], -dx);
		atomicAdd(&_right_hand_y[i_constraints[3]], -dy);
	}
	//this constraint assigns the direction to this edge
	//(f_constraints[0],f_constraints[1]) is the unit direction vector
	//f_constraints[4] is the weight para
	else {
		atomicAdd(&_right_hand_x[i_constraints[2]], 3.0f*f_constraints[4] * f_constraints[0]);
		atomicAdd(&_right_hand_y[i_constraints[2]], 3.0f*f_constraints[4] * f_constraints[1]);
		atomicAdd(&_right_hand_x[i_constraints[3]], -3.0f*f_constraints[4] * f_constraints[0]);
		atomicAdd(&_right_hand_y[i_constraints[3]], -3.0f*f_constraints[4] * f_constraints[1]);
	}

	return 3;
}

__forceinline__ __device__ int StressKernel(
	void * __restrict__ constraints,
	float* __restrict__ P_Opt_x,
	float* __restrict__ P_Opt_y,
	float* __restrict__ _right_hand_x,
	float* __restrict__ _right_hand_y) { }
__device__ Device_Stub stubs[] = { CrossingRemovalKernel, StressKernel, CircleKernel, EqualAngleKernel, NoverlapKernel, SymmetryKernel, AtomKernel }; ///

__global__ void LocalSteps(cuDWORD * __restrict__ constraints, float* __restrict__ P_Opt_x, float* __restrict__ P_Opt_y, float* __restrict__ _right_hand_x, float* __restrict__ _right_hand_y,
	int *constraints_offset, int const_num) {//float *right_hand
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	for (; i<const_num; i += gridDim.x * blockDim.x)
	{
		cuDWORD* i_constraints = constraints + constraints_offset[i];
		stubs[*i_constraints](i_constraints + 2, P_Opt_x, P_Opt_y, _right_hand_x, _right_hand_y);//i**(constraints + 1) *(constraints + 2)
	}
}

/*For Stress Constraints*/
__global__ void local_stress(
	float* __restrict__ P_Opt_x, float* __restrict__ P_Opt_y,
	float* __restrict__ _right_hand_x, float* __restrict__ _right_hand_y,

	float *__restrict__ dist, int n
) {
	unsigned stepping = blockDim.x;
	unsigned j = blockIdx.x;
	static __shared__ float shared[64];

	float sum = 0, sum2 = 0;
	unsigned wid = threadIdx.x >> 5;	// warp ID
	unsigned lane = threadIdx.x - (wid << 5);
#pragma unroll
	for (; j < n; j += gridDim.x) {
		sum = 0; sum2 = 0;
#pragma unroll
		for (unsigned i = threadIdx.x; i < n; i += stepping)
		{
			if (i != j) {

				float t = rsqrt(pow(_right_hand_x[j] - _right_hand_x[i], 2.f) + pow(_right_hand_y[j] - _right_hand_y[i], 2.f));

				if (isinf(t) || isnan(t)) t = 1;
				sum += ((t / dist(j, i))) * (_right_hand_x[j] - _right_hand_x[i]);
				sum2 += ((t / dist(j, i))) * (_right_hand_y[j] - _right_hand_y[i]);
			}
		}

#pragma unroll
		for (unsigned offset = 16; offset > 0; offset >>= 1)
		{
			sum += __shfl_down(sum, offset);
			sum2 += __shfl_down(sum2, offset);
		}
		if (lane == 0)
		{
			shared[wid + wid] = sum;
			shared[wid + wid + 1] = sum2;
		}
		__syncthreads();

		sum = (threadIdx.x < 32) ? shared[lane + lane] : 0;// group to warp 1;
		sum2 = (threadIdx.x < 32) ? shared[lane + lane + 1] : 0;// group to warp 1;

		if (wid == 0) {
#pragma unroll
			for (unsigned offset = 16; offset > 0; offset >>= 1)
			{
				sum += __shfl_down(sum, offset);
				sum2 += __shfl_down(sum2, offset);
			}
		}
		if (threadIdx.x == 0)
		{
			P_Opt_x[j] = sum;
			P_Opt_y[j] = sum2;
		}
	}
	sum = 0; sum2 = 0;
}
__global__ void cuMemchecker(float *mem, float *mem2) {
	mem[0] = 0;
}
__global__ void jacobiPreconditioner(float *M, float *diag, int n, int innersize) {
	//	int i = threadIdx.x;
	int j = blockIdx.x;
	for (; j < n; j += gridDim.x) {
		for (int i = threadIdx.x; i < n; i += blockDim.x) {
			if (i != j) {
				M[j*innersize + i] /= M[j*innersize + j];
			}
		}
		__syncthreads();
		if (threadIdx.x == 0)
		{
			diag[j] = M[j*innersize + j];
			M[j*innersize + j] = 1;
		}
	}
}
__global__ void jacobiPrecondImpl(float *diag, float *x, float *y, int n) {
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	if (i < n) {
		x[i] /= diag[i];
		y[i] /= diag[i];

	}
	else if (i<2 * n) {

	}

}
DWORD WINAPI finalize(LPVOID);


class CudaCore {

private:
	viennacl::matrix<float> system_matrix;
	viennacl::vector<float> right_hand, right_hand_x, right_hand_y;
	viennacl::vector<float> result, result_x, result_y;
	float *d_right_x;
	float *d_right_y;
	float* res_x;
	float* res_y;
	bool runonce = true;
	viennacl::linalg::cg_solver<viennacl::vector<float>> *solver;
	void *d_mat;
	cudaStream_t nonBlocking;

	cuDWORD *d_constraints;

	int *d_const_idx;
	float *d_dist = 0;
	float *diag;
	//J
	float *d_J_val;
	int *d_J_RowPtr, *d_J_ColIdx;
	int n_data = -1, constraints_size, constraints_length, val_size;
	int const_num;

	int blocks, threads = 128;

public:
	cuDWORD *d_gc1;
	int *d_gc2;
	int MaxIter = 100;

	CudaCore()
	{}

	//
	void init(float* sysm, int n_data, cuDWORD* constraints, int* const_idx, int constraint_num, float **dists)
	{

		if (d_dist == 0 || this->n_data != n_data) {
			cudaStreamCreateWithFlags(&nonBlocking, cudaStreamNonBlocking);

			this->n_data = n_data; // n_data is the rows of L, = 2 * nodenum 

			system_matrix = viennacl::matrix<float>(n_data, n_data, viennacl::context(viennacl::CUDA_MEMORY));
			d_mat = system_matrix.handle().cuda_handle().get();
			cudaMalloc(&diag, n_data * sizeof(float));
			cudaMalloc(&d_dist, sizeof(float) * (n_data*(n_data - 1)) / 2);
			float *d_dist_ptr = d_dist;
			for (int i = 0; i < n_data; i++) {
				int to_cpy = n_data - i - 1;
				cudaMemcpyAsync(d_dist_ptr, dists[i] + i + 1, sizeof(float) * to_cpy, cudaMemcpyHostToDevice, nonBlocking);
				d_dist_ptr += to_cpy;
			}

			solver = new viennacl::linalg::cg_solver<viennacl::vector<float>>(viennacl::linalg::cg_tag(1e-5, 15));//, 80 for 12000

			right_hand_x = viennacl::vector<float>(n_data, viennacl::context(viennacl::CUDA_MEMORY));//d_right
			right_hand_y = viennacl::vector<float>(n_data, viennacl::context(viennacl::CUDA_MEMORY));//d_right

			d_right_x = (float *)right_hand_x.handle().cuda_handle().get();
			d_right_y = (float *)right_hand_y.handle().cuda_handle().get();

			result_x = viennacl::vector<float>(n_data, viennacl::context(viennacl::CUDA_MEMORY));//d_right;
			result_y = viennacl::vector<float>(n_data, viennacl::context(viennacl::CUDA_MEMORY));//d_right;

			res_x = (float*)result_x.handle().cuda_handle().get();
			res_y = (float*)result_y.handle().cuda_handle().get();
		}
		const_num = constraint_num;
		constraints_length = const_idx[const_num];
		cudaMalloc(&d_constraints, constraints_length * sizeof(cuDWORD));
		cudaMalloc(&d_const_idx, const_num * sizeof(int));
		cudaMemcpy(d_constraints, constraints, constraints_length * sizeof(cuDWORD), cudaMemcpyDefault);
		//	delete[] constraints;
		cudaMemcpy(d_const_idx, const_idx, sizeof(int)*const_num, cudaMemcpyDefault);
		delete[] const_idx;


		int internal_size = system_matrix.internal_size1();
		for (int i = 0; i < n_data; i++) {
			cudaMemcpyAsync(
				((cuDWORD *)d_mat) + i *internal_size,
				sysm + n_data * i, n_data * sizeof(float),
				cudaMemcpyHostToDevice// , nonBlocking
			);
		}
		//cudaStreamSynchronize(nonBlocking);
		//	jacobiPreconditioner<<<128,128>>>((float*)d_mat, diag, n_data, internal_size);
	}
	int t = 0;


	void Solve(float *data_x, float *data_y,
		float *out_x, float *out_y, int iters, int* es, float* shortest_path) {
		if (runonce) {
			cudaMemcpy(res_x, out_x, sizeof(float)*n_data, cudaMemcpyHostToDevice);
			cudaMemcpy(res_y, out_y, sizeof(float)*n_data, cudaMemcpyHostToDevice);
			runonce = false;
		}
		std::cout << "iterations: " << iters << std::endl;
		for (int i = 0; i < iters; i++) {
			if (i == 2)
				solver->tag_ = viennacl::linalg::cg_tag(1e-5, 13);
			else if (i == 10)//18-19
				solver->tag_ = viennacl::linalg::cg_tag(1e-5, 10);
			else if (i == 15)
				solver->tag_ = viennacl::linalg::cg_tag(1e-5, 8);

			cudaMemset(d_right_x, 0, sizeof(float)*(n_data));
			cudaMemset(d_right_y, 0, sizeof(float)*(n_data));
			local_stress << <128, 1024 >> >(d_right_x, d_right_y, res_x, res_y, d_dist, n_data);
			LocalSteps << <128, 128 >> >(d_constraints, res_x, res_y, d_right_x, d_right_y, d_const_idx, const_num); //

			(*solver)(system_matrix, right_hand_x, right_hand_y);

			cudaMemset(res_x, 0, sizeof(float)*(n_data));
			cudaMemset(res_y, 0, sizeof(float)*(n_data));
			local_stress << <128, 1024 >> >(res_x, res_y, d_right_x, d_right_y, d_dist, n_data);
			LocalSteps << <128, 128 >> >(d_constraints, d_right_x, d_right_y, res_x, res_y, d_const_idx, const_num); //

			(*solver).operator()<decltype(system_matrix), true>(system_matrix, result_x, result_y);


		}

		cudaMemcpyAsync(
			out_x,
			res_x,
			n_data * sizeof(float),
			cudaMemcpyDeviceToHost
		);
		cudaMemcpyAsync(
			out_y,
			res_y,
			n_data * sizeof(float),
			cudaMemcpyDeviceToHost
		);
		/*
		d_gc1 = d_constraints;
		d_gc2 = d_const_idx;

		CreateThread(0, 0, finalize, this, 0, 0);
		*//*d_constraints = 0;
		d_const_idx = 0;*/

		cudaFree(d_constraints);

		d_constraints = 0;

		cudaFree(d_const_idx);

		d_const_idx = 0;
	}


	~CudaCore() {
		if (!right_hand_x.empty())
		{

			if (!d_constraints)

				cudaFree(d_constraints);

			d_constraints = 0;

			if (!d_const_idx)

				cudaFree(d_const_idx);

			d_const_idx = 0;

			if (!d_dist)

				cudaFree(d_dist);

			d_dist = 0;

			delete solver;
		}
	}
};
DWORD WINAPI finalize(LPVOID param) {
	CudaCore* core = ((CudaCore*)param);
	cudaFree(core->d_gc1);
	cudaFree(core->d_gc2);
	return 0;
}
CudaSolver::CudaSolver() {
	core = new CudaCore();

}
CudaSolver::~CudaSolver() {
	delete core;
}

void CudaSolver::init(float *mat, int n_mat, cuDWORD* constraints, int* const_idx, int constraint_size, float **dists) {
	core->init(mat, n_mat, constraints, const_idx, constraint_size, dists);
}

void CudaSolver::Solve(float *right_hand_x, float *right_hand_y, float* P_Out_x, float* P_Out_y, int iters,
	int* es, float* shortest_path) {
	core->Solve(right_hand_x, right_hand_y, P_Out_x, P_Out_y, iters, es, shortest_path);
}



void CudaSolver::setMaxIter(int mi) {
	core->MaxIter = mi;
}
__inline__ __device__  void _rotate(float alpha, float &edge_X, float &edge_Y) {
	float thow = (sqrtf(edge_X * edge_X + edge_Y * edge_Y));
	float p_angle = 0; // x==0, y==0
	if (edge_X > 0) {
		p_angle = atanf(edge_Y / edge_X);
	}
	else if (edge_X < 0 && edge_Y >= 0) {
		p_angle = atanf(edge_Y / edge_X) + 3.14159;
	}
	else if (edge_X< 0 && edge_Y < 0) {
		p_angle = atanf(edge_Y / edge_X) - 3.14159;
	}
	else if (edge_X == 0 && edge_Y < 0) {
		p_angle = -3.14159 / 2;
	}
	else if (edge_X == 0 && edge_Y > 0) {
		p_angle = 3.14159 / 2;
	}
	edge_X = thow * cos(p_angle - alpha);
	edge_Y = thow * sin(p_angle - alpha);
}
__inline__ __device__  void moveCenterToZero(float center_x, float &edge_X, float &edge_Y) {
	edge_X -= center_x;
}
__inline__ __device__ void ShapeMatching(float *relative_shape,
	float *P_Opt_x, float *P_Opt_y, int *pid, int node_num) {

	float X0cmx = 0.0f, X0cmy = 0.0f, Xcmx = 0.0f, Xcmy = 0.0f;
	for (int i = 0; i < node_num; i++) {
		X0cmx += relative_shape[2 * i + 0];
		X0cmy += relative_shape[2 * i + 1];

		Xcmx += P_Opt_x[pid[i]];
		Xcmy += P_Opt_y[pid[i]];
	}
	X0cmx /= node_num;
	X0cmy /= node_num;
	Xcmx /= node_num;
	Xcmy /= node_num;

	float Apqa = 0.0f, Apqb = 0.0f, Apqc = 0.0f, Apqd = 0.0f;
	for (int i = 0; i < node_num; i++) {
		float qix = relative_shape[2 * i + 0] - X0cmx;
		float qiy = relative_shape[2 * i + 1] - X0cmx;
		float pix = P_Opt_x[pid[i]] - Xcmx;
		float piy = P_Opt_y[pid[i]] - Xcmy;

		Apqa += pix*qix;
		Apqb += pix*qiy;
		Apqc += piy*qix;
		Apqd += piy*qiy;
	}

	float Sa, Sb, Sc, Sd;
	Sa = Apqa*Apqa + Apqc*Apqc;
	Sb = Apqa*Apqb + Apqc*Apqd;
	Sc = Apqa*Apqb + Apqc*Apqd;
	Sd = Apqb*Apqb + Apqd*Apqd;

	//matrix sqrt
	float tao_real = Sa + Sd, tao_i = 0;
	float sigma_real = Sa*Sd - Sb*Sc, sigma_i = 0;
	float s_real = 0, s_i = 0;
	if (sigma_real >= 0) {
		s_real = sqrt(sigma_real);
	}
	else {
		s_i = sqrt(abs(sigma_real));
	}

	float t_real = tao_real + 2 * s_real, t_i = tao_i + 2 * s_i;
	if (t_real >= 0) {
		t_real = sqrt(t_real);
	}
	else {
		t_i = sqrt(abs(t_real));
	}

	float _t_real = t_real / (t_real*t_real - t_i*t_i), _t_i = -t_i / (t_real*t_real - t_i*t_i);
	Sa = _t_real*(Sa + s_real) - _t_i*s_i;
	Sb = _t_real*(Sb);
	Sc = _t_real*(Sc);
	Sd = _t_real*(Sd + s_real) - _t_i*s_i;



	float ISa, ISb, ISc, ISd;
	float S_trace = Sa*Sd - Sb*Sc;
	if (S_trace == 0.0f) {
		S_trace = 0.0001;
	}
	ISa = Sd / S_trace;
	ISb = -Sb / S_trace;
	ISc = -Sc / S_trace;
	ISd = Sa / S_trace;

	float Ra, Rb, Rc, Rd;
	Ra = Apqa*ISa + Apqb*ISc;
	Rb = Apqa*ISb + Apqb*ISd;
	Rc = Apqc*ISa + Apqd*ISc;
	Rd = Apqc*ISb + Apqd*ISd;

	for (int i = 0; i < node_num; i++) {
		float Xi0x = relative_shape[2 * i + 0];
		float Xi0y = relative_shape[2 * i + 1];

		float gix = Ra*(Xi0x - X0cmx) + Rb*(Xi0y - X0cmy) + Xcmx;
		float giy = Rc*(Xi0x - X0cmx) + Rd*(Xi0y - X0cmy) + Xcmy;

		relative_shape[2 * i + 0] = gix;
		relative_shape[2 * i + 1] = giy;
	}
}
