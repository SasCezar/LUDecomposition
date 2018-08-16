#include <time.h>
#include <sys/time.h>
#include "mpi.h"
#include "../ludutils.h"

/*
 * forw_elim - forward Gauss elimination
 *
 * @origin row pointer by reference
 * @master_row row in which lays diagonal
 */
void forw_elim(float **origin, float *master_row, size_t dim)
{
   if (**origin == 0)
      return;

   float k = **origin / master_row[0];

   int i;
   for (i = 1; i < dim; i++) {
      (*origin)[i] = (*origin)[i] - k * master_row[i];
   }
   **origin = k;
}


int main(int argc, char *argv[])
{
	const int root_p = 0;
    int n = 0, p, id;

    n = atol(argv[1]);
    char path[255];
	int out = snprintf(path, 255, "/mnt/c/Users/sasce/Desktop/Matrices/matrix_%i.csv", n);
	printf("Loading matrix\n");
	Matrix original = read_csv(path);
    Matrix matrix = duplicate_matrix(original);
	float **A = matrix.matrix;
	printf("Matrix loaded\n");
	printf("Decomposing matrix\n");

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    struct timeval tv;
    double start, end;
	gettimeofday(&tv,NULL);
	start = tv.tv_sec;

  
    for (int i = 0; i < n - 1; i++) {
        float *diag_row = &A[i][i];
        for (int j = i + 1; j < n; j++) {
            if (j % p == id) {
                float *save = &A[j][i];
                forw_elim(&save, diag_row, n - i);
            }
        }

        for (int j = i + 1; j < n; j++) {
            float *save = &A[j][i];
            MPI_Bcast(save, n - i, MPI_FLOAT, j % p, MPI_COMM_WORLD);
        }
    }

	if(id == 0){
		gettimeofday(&tv,NULL);
    	end=tv.tv_sec;

		double time_spent = (double)(end - start);
		printf("Matrix decomposed - Elapsed %f\n", time_spent);

		LU result = split_lu(A, n);
		printf("Recomposing original matrix\n");
		Matrix recomposed = matrix_multiplication(result.L, result.U);
		printf("Matrix recomposed\n");
		printf("Computing error\n");
		double error = compute_error(original, recomposed);
		printf("error %f\n", error);
	}	

	MPI_Finalize();

	return EXIT_SUCCESS;
}