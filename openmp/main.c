// #include <iostream>
#include <time.h>
#include<sys/time.h>
// #include "lud_omp.h"
// #include "lud.h"
#include "lud_mpi.h"
#include "ludutils.h"
//#include "ludcrout.h"


int main(int argc, char *argv[]){
	int rank, nprocs;
    int i=0, j=0, k=0;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		printf("Loading matrix\n");
		Matrix original = read_csv("/mnt/c/Users/sasce/Desktop/Matrices/matrix_6000.csv");
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

	MPI_Finalize();
	
	return 0;
}