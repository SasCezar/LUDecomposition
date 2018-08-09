// #include <iostream>
#include <time.h>
#include<sys/time.h>
// #include "lud_omp.h"
// #include "lud.h"
#include "lud_mpi.h"
#include "ludutils.h"
//#include "ludcrout.h"


int main(int argc, char *argv[]){
	printf("Loading matrix\n");
	Matrix original = read_csv("/mnt/c/Users/sasce/Desktop/Matrices/matrix_200.csv");
	Matrix matrix = duplicate_matrix(original);
	printf("Matrix loaded\n");
	// print_matrix(matrix);
	printf("Decomposing matrix\n");

	struct timeval tv;
	gettimeofday(&tv,NULL);
  	double start=tv.tv_sec;
	LU result = decompose_mpi(matrix);
	gettimeofday(&tv,NULL);
    double end=tv.tv_sec;
	double time_spent = (double)(end - start);
	printf("Matrix decomposed - Elapsed %f\n", time_spent);

	// printf("Lower\n");
	// print_matrix(result.L);
	// printf("Upper\n");
	// print_matrix(result.U);

	printf("Recomposing original matrix\n");
	Matrix recomposed = matrix_multiplication(result.L, result.U);
	printf("Matrix recomposed\n");
	// printf("Recomposed\n");
	// print_matrix(recomposed);

	printf("Computing error\n");
	double error = compute_error(original, recomposed);
	printf("error %f\n", error);
	free(matrix.matrix);

	return 0;
}