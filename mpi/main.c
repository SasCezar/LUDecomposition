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

  
	int *map = malloc(size * (sizeof *map));

	for(i=0; i<size; i++)
	{
		map[i]= i % nprocs;
	}

	for(j =0; j < size-1; j++) {
        if(map[j] == rank)
        {
        	for(i = j+1; i < size; i++) {
        		A[i][j] = A[i][j]/A[j][j];
        	}
        }
        MPI_Bcast (&A[j][j],size-j,MPI_DOUBLE,map[j],MPI_COMM_WORLD);

        		for(k = j+1; k < size; k++)
        		{
        			if(map[j] == rank)
        			{
        				for(i = j+1; i < size; i++) {
        					A[k][i]= A[k][i] - (A[k][j] * A[j][i]);
       				    }
       			}
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