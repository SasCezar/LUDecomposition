// #include <iostream>
#include <time.h>
#include<sys/time.h>
// #include "lud_omp.h"
// #include "lud.h"
#include "lud_mpi.h"
#include "../ludutils.h"
//#include "ludcrout.h"


int main(int argc, char *argv[]){
	int rank, nprocs;
    int i=0, j=0, k=0;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	Matrix original;
	Matrix matrix;
	struct timeval tv;
	double start, end;
	if(rank == 0){
		printf("Loading matrix\n");
		original = read_csv("/mnt/c/Users/sasce/Desktop/Matrices/matrix_2000.csv");
		matrix = duplicate_matrix(original);
		printf("Matrix loaded\n");
		// print_matrix(matrix);
	
		printf("Decomposing matrix\n");

		gettimeofday(&tv,NULL);
		start=tv.tv_sec;
	}
	LU result = decompose_mpi(matrix);
	if(rank == 0){
		gettimeofday(&tv,NULL);
		end=tv.tv_sec;
		double time_spent = (double)(end - start);
		printf("Matrix decomposed - Elapsed %f\n", time_spent);

		printf("Recomposing original matrix\n");
		Matrix recomposed = matrix_multiplication(result.L, result.U);
		printf("Matrix recomposed\n");

		printf("Computing error\n");
		double error = compute_error(original, recomposed);
		printf("error %f\n", error);
		free(matrix.matrix);
	}
	MPI_Finalize();
	
	return 0;
}