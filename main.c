#include <iostream>
#include <time.h>
 #include "lud_omp.h"
// #include "lud.h"
#include "ludutils.h"
//#include "ludcrout.h"


int main(int argc, char *argv[]){
	printf("Loading matrix\n");
	Matrix matrix = read_csv("Matrices/matrix_2000.csv");
	printf("Matrix loaded\n");
	// print_matrix(matrix);
	printf("Decomposing matrix\n");
	clock_t begin = clock();
	LU result = decompose_omp(matrix);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
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
	double error = compute_error(matrix, recomposed);
	printf("error %f\n", error);
	free(matrix.matrix);

	return 0;
}