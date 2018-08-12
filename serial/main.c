// #include <iostream>
#include <time.h>
#include <sys/time.h>
#include "lud.h"
#include "../ludutils.h"

int main(int argc, char *argv[]){
	printf("Loading matrix\n");
	int size = atoi(argv[1]); 
	char path[255];
	int out = snprintf(path, 255, "/mnt/c/Users/sasce/Desktop/Matrices/matrix_%i.csv", size);
	Matrix original = read_csv(path);
	Matrix matrix = duplicate_matrix(original);
	printf("Matrix loaded\n");
	printf("Decomposing matrix\n");

	struct timeval tv;
	gettimeofday(&tv,NULL);
	double start=tv.tv_sec;
	LU result = decompose_serial(matrix);
	gettimeofday(&tv,NULL);
	double end=tv.tv_sec;
	double time_spent = (double)(end - start);
	printf("Matrix decomposed - Elapsed %f\n", time_spent);

	printf("Recomposing original matrix\n");
	Matrix recomposed = matrix_multiplication(result.L, result.U);
	printf("Matrix recomposed\n");
	
	printf("Computing error\n");
	double error = compute_error(original, recomposed);
	printf("error %f\n", error);
	
	return 0;
}