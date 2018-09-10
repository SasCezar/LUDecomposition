
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <chrono>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include "parser.hpp"

using namespace std;
using namespace aria::csv;


cudaError_t decomposeCuda(float *A, int size, int threads_per_block, int num_blocks);


float * read_csv(string path, int N) {
	/**
	 * Reads a csv files and creates a matrix
	*/

	std::ifstream f(path);
	CsvParser parser = CsvParser(f).delimiter(',');

	int i = 0, j = 0;
	float *A = (float *)malloc(sizeof(float) * N * N);
	for (auto& row : parser) {
		j = 0;
		for (auto& field : row) {
			A[i * N + j] = std::stof(field);
			j++;
		}
		i++;
	}

	return A;
}

float **matrix_create(size_t m, size_t n) {
	float **result = (float **)malloc(sizeof(float *)*m);

	for (int i = 0; i < m; i++) {
		result[i] = (float *)malloc(sizeof(float)*n);
	}
	return result;
}

float **matrix_difference(float **ma, float **mb, int n) {

	float **mc = matrix_create(n, n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			mc[i][j] = ma[i][j] - mb[i][j];
		}
	}

	return mc;
}

float **matrix_multiplication(float **ma, float **mb, int n) {
	float **mc = matrix_create(n, n);

	float sum = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < n; k++) {
				sum = sum + ma[i][k] * mb[k][j];
			}
			mc[i][j] = sum;
		}
	}
	return mc;
}


double frobenius_norm(float **A, int n) {
	float sum = 0;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sum += A[i][j] * A[i][j];
		}
	}

	double norm = sqrt(sum);
	return norm;
}


typedef struct {
	float **L;
	float **U;
} LU;


float **initialize_matrix(int rows, int cols) {

	float **matrix = matrix_create(rows, cols);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			matrix[i][j] = 0;
		}
	}
	return matrix;
}

LU split_lu(float **a, int n) {
	float **L = initialize_matrix(n, n);
	float **U = initialize_matrix(n, n);


	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i < j)
			{
				U[i][j] = a[i][j];
			}
			if (i > j) {
				L[i][j] = a[i][j];
			}
			if (i == j) {
				L[i][j] = 1;
				U[i][j] = a[i][j];
			}
		}
	}


	LU decomposition = LU();
	decomposition.L = L;
	decomposition.U = U;
	return decomposition;
}

double compute_error(float **a, float **b, int n) {
	float **difference = matrix_difference(a, b, n);
	double error = frobenius_norm(difference, n);
	return error;
}



__constant__ int MATRIX_SIZE;


__global__ void decompose_multipliers(float *A, int rows_per_thread, int i) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	int jstart = (i + 1) + tid * rows_per_thread;
	int jend = jstart + rows_per_thread;

	for (int j = jstart; j < jend && j < MATRIX_SIZE; j++) {
		A[j * MATRIX_SIZE + i] = A[j * MATRIX_SIZE + i] / A[i * MATRIX_SIZE + i]; // Computes the multipliers and updates L in A
	}
}

__global__ void decompose_elimination(float *A, int rows_per_thread, int i) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	int jstart = (i + 1) + tid * rows_per_thread;
	int jend = jstart + rows_per_thread;

	for (int j = jstart; j < jend && j < MATRIX_SIZE; j++) { // Iterates over the remaining rows
		for (int k = i + 1; k < MATRIX_SIZE; k++) { // iterates over the remaining columns
			A[j * MATRIX_SIZE + k] -= A[j * MATRIX_SIZE + i] * A[i * MATRIX_SIZE + k]; // Updates U in A
		}
	}

}

void print_matrix(float **matrix, int n) {
	/**
	 * Prints the matrix
	 * @param matrix The matrix to print
	*/


	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%f\t", matrix[i][j]);
		}
		printf("\n");
	}
}

float **matrix2d(float *A, int n) {
	float **original = matrix_create(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			original[i][j] = A[i*n + j];
		}
	}
	return original;
}



int main(int argc, char const *argv[])
{
	int size = atoi(argv[1]);
	//int size = 2000;
	char path[255];
	int out = snprintf(path, 255, "C:\\Users\\sasce\\Desktop\\Matrices\\matrix_%i.csv", size);

	float *A = read_csv(path, size);
	float **original = matrix2d(A, size);

	struct cudaDeviceProp properties;
	cudaGetDeviceProperties(&properties, 0);
	// cout << "using " << properties.multiProcessorCount << " multiprocessors" << endl;
	// cout << "max threads per processor: " << properties.maxThreadsPerMultiProcessor << endl;

	// Decomplse matrix in parallel.
	int threads_per_block = atoi(argv[2]);
	int num_blocks = atoi(argv[3]);
	// int threads_per_block = 512;
	// int num_blocks = 2;
	auto start = std::chrono::high_resolution_clock::now();
	cudaError_t cudaStatus = decomposeCuda(A, size, threads_per_block, num_blocks);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "decomposeCuda failed!");
		return 1;
	}

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	cout << size << ";" << num_blocks << ";" << threads_per_block << ";" << (int)elapsed.count() << "\n";

	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}


	/*
	float **result = matrix2d(A, size);
	LU decomposition = split_lu(result, size);
	float **recomposed = matrix_multiplication(decomposition.L, decomposition.U, size);
	double error = compute_error(recomposed, original, size);
	printf("Error %f", error);
	*/
	return 0;
}

// Helper function for using CUDA to decompose matrix in parallel.
cudaError_t decomposeCuda(float *A, int size, int threads_per_block, int num_blocks)
{
	float *dev_a;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed! Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers for three vectors (two input, one output).
	cudaStatus = cudaMalloc((void**)&dev_a, size * size * sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_a, A, size * size * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


	// Copy input size from host memory to GPU.
	cudaStatus = cudaMemcpyToSymbol(MATRIX_SIZE, &size, sizeof(size));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyToSymbol failed!");
		goto Error;
	}

	// Launch a kernel on the GPU with one thread for each element.

	int ops_per_thread = ceil((size) / (float)(threads_per_block*num_blocks));

	dim3 thread_block(threads_per_block, 1, 1);
	dim3 grid(num_blocks, 1);
	// printf("Ops per thread %i \n", ops_per_thread);

	for (int i = 0; i < size; i++) { // Iterates over the columns to remove
		decompose_multipliers << <grid, thread_block >> > (dev_a, ops_per_thread, i);
		decompose_elimination << <grid, thread_block >> > (dev_a, ops_per_thread, i);
	}

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "decomposeKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}


	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching decomposeKernel!\n", cudaStatus);
		goto Error;
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(A, dev_a, size * size * sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

Error:
	cudaFree(dev_a);

	return cudaStatus;
}