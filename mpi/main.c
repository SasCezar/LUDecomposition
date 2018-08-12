#include <time.h>
#include <sys/time.h>
#include "mpi.h"
#include "../ludutils.h"


int main(int argc, char *argv[])
{
    int rank, nprocs;     // for storing this process' rank, and the number of processes

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	double start, end;
    Matrix original, matrix; 
	float **A;
	int n = atoi(argv[1]);  
    
	if(rank == 0){
        char path[255];
		int out = snprintf(path, 255, "/mnt/c/Users/sasce/Desktop/Matrices/matrix_%i.csv", n);
		printf("Loading matrix\n");
		original = read_csv(path);
        matrix = duplicate_matrix(original);
		A = matrix.matrix;
		printf("Matrix loaded\n");
		printf("Decomposing matrix\n");
		start = MPI_Wtime();
    }else{
        A = matrix_create(n, n);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int rem = n % nprocs;       // Elements remaining after division among processes
    int sum = 0;                // Sum of counts. Used to calculate displacements
    
    int *sendcounts = malloc(sizeof(int)*nprocs); // Array describing how many elements to send to each process
    int *displs = malloc(sizeof(int)*nprocs); // Array describing the displacements where each segment begins

    // Calculate send counts and displacements
    for (int i = 0; i < nprocs; i++) {
        sendcounts[i] = n / nprocs;

        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }
        sendcounts[i] *= n;

        displs[i] = sum;
        sum += sendcounts[i];
    }

	int rows = sendcounts[rank] / n;
	float  **rec_buf = matrix_create(rows, n);

    // Print calculated send counts and displacements for each process
    if (0 == rank) {
        for (int i = 0; i < nprocs; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    // Divide the data among processes as described by sendcounts and displs
    
    MPI_Scatterv(&A[0][0], sendcounts, displs, MPI_FLOAT, &rec_buf[0][0], sendcounts[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);
    
	int mystart = displs[rank] / n;
    int myend = mystart + rows - 1;
	
	if(rank == 0){
        printf("rank, i\tj\tk\n");
    }
	MPI_Barrier(MPI_COMM_WORLD);
	
	printf("Rank %i: start %i - end %i", rank, mystart, myend);
	for (int i = mystart; i < myend; i++) { // Iterates over the columns to remove
        for(int j = i + 1; j < rows; j++){ // Iterates over the remaining rows
            float m = rec_buf[j][i] / rec_buf[i][i];
            for(int k = i+1; k < rows; k++){ // Iterates over the remaining columns
                printf("%i\t%i\t%i\t%i\n",rank, i, j, k);
				rec_buf[j][k] -= m * rec_buf[i][k];
            }
            rec_buf[j][i] = m;
        }
    }
	
	printf("RANK END %i\n", rank);
	
    MPI_Gatherv(&rec_buf[0][0], sendcounts[rank], MPI_FLOAT, &A[0][0], sendcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

	printf("Gathered\n");

	if(rank == 0){
		end = MPI_Wtime();
		double time_spent = (double)(end - start);
		printf("Matrix decomposed - Elapsed %f\n", time_spent);

		LU result = split_lu(A, n);
		printf("Recomposing original matrix\n");
		Matrix recomposed = matrix_multiplication(result.L, result.U);
		printf("Matrix recomposed\n");
		//Matrix recomposed = {.matrix = A, .n=n};
		printf("Computing error\n");
		double error = compute_error(original, recomposed);
		printf("error %f\n", error);
	}	

	MPI_Finalize();

	return 0;
}