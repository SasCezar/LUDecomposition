#include <time.h>
#include<sys/time.h>
#include "lud_omp.h"
#include "../ludutils.h"


int main(int argc, char *argv[]){
	// printf("Loading matrix\n");
	// int size = atoi(argv[1]); 
	char path[255];
	printf("size;num_threads;elapsed\n");
	for(int size = 1000; size <= 9000; size += 1000){
		int out = snprintf(path, 255, "/mnt/c/Users/sasce/Desktop/Matrices/matrix_%i.csv", size);
		Matrix original = read_csv(path, size);
		for(int n = 2; n <= 8; n ++){
			// for(int j = 0; j < 9; j++){
			{
				Matrix matrix = duplicate_matrix(original);
				// printf("Matrix loaded\n");
				// printf("Decomposing matrix\n");

				struct timeval tv;
				gettimeofday(&tv,NULL);
				double start=tv.tv_sec;
				LU result = decompose_omp(matrix, n);
				gettimeofday(&tv,NULL);
				double end=tv.tv_sec;
				int time_spent = (int)(end - start);
				printf("%i;%i;%i\n", size, n, time_spent);
				free(matrix.matrix);
				// printf("Matrix decomposed - Elapsed %f\n", time_spent);
			}
		}
	}

	// printf("Recomposing original matrix\n");
	// Matrix recomposed = matrix_multiplication(result.L, result.U);
	// printf("Matrix recomposed\n");

	// print_matrix(result.L);
	// printf("\n\n\n");
	// print_matrix(result.U);

	// printf("Computing error\n");
	// double error = compute_error(original, recomposed);
	// printf("error %f\n", error);
	
	return 0;

	return 0;
}