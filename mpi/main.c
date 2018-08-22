#include <time.h>
#include <sys/time.h>
#include "mpi.h"
#include "../ludutils.h"

float *matrix1d(float **matrix, int n){
	float *r = malloc(sizeof(float) * n * n);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			r[i * n + j] = matrix[i][j];
		}
	}
	return r;
}

float **matrix2d(float *matrix, int n){
	float **r = matrix_create(n, n);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			r[i][j] = matrix[i * n + j];
		}
	}
	return r;
}

int main(int argc, char *argv[])
{
	int nprocs, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n = atol(argv[1]);
    char path[255];
	int out = snprintf(path, 255, "/mnt/c/Users/sasce/Desktop/Matrices/matrix_%i.csv", n);
	printf("Loading matrix\n");
	
	Matrix original = read_csv(path);
    Matrix matrix = duplicate_matrix(original);
	float *l = matrix1d(identity(n),n);
	float *A = matrix1d(matrix.matrix, n);

	printf("Matrix loaded\n");
	printf("Decomposing matrix\n");
    
    struct timeval tv;
    double start, end;
	gettimeofday(&tv,NULL);
	start = tv.tv_sec;
	
	int *map = malloc(n * sizeof(int));
  	for(int i = 0; i < n; i++)
	{
		map[i] = i % nprocs;
	}
	
	for(int k = 0; k < n - 1 ; k++) {
		
		MPI_Bcast(&A[k * n + k], n - k, MPI_FLOAT, rank, MPI_COMM_WORLD);
		
		for(int i = k + 1; i < n; i++) {
			if(map[i] == rank){
				l[i * n + k] = A[i * n + k] / A[k * n + k];
			}	
		}
        for(int j = k + 1; j < n; j++){
			for(int i = k + 1; i < n; i++) {
				if(map[i] == rank){
					A[i * n + j] -= l[i * n + k] * A[k * n + j];	
				}
			}
       	}
    }
	
	if(rank == 0){
		gettimeofday(&tv,NULL);
    	end=tv.tv_sec;

		double time_spent = (double)(end - start);
		printf("Matrix decomposed - Elapsed %f\n", time_spent);

		float **r = matrix2d(A, n);
		Matrix L = {.matrix = matrix2d(l, n), .n = n};
		LU result = split_lu(r, n);
		printf("Recomposing original matrix\n");

		print_matrix(L);
		printf("\n\n\n");
		print_matrix(result.U);

		Matrix recomposed = matrix_multiplication(L, result.U);
		printf("Matrix recomposed\n");
		printf("Computing error\n");
		double error = compute_error(original, recomposed);
		printf("error %f\n", error);
	}	

	MPI_Finalize();

	return EXIT_SUCCESS;
}