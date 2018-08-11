#ifndef LUDMPI_H
#define LUDMPI_H

#include "../ludutils.h"
#include "mpi.h"

int map(int i, int n) {
  return i % n;
}

LU decompose_mpi(Matrix matrix){
    int n = matrix.n;
    float **A = matrix.matrix;
    
    int rank, nprocs;
    int i=0, j=0, k=0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);	/* get number of processes */

    int nrows = n / nprocs;

	int process_row_min =pid * nrows;
	int process_row_max = process_row_min + nrows - 1;

    printf("Process %d\n", rank);

    for (int i = 0; i < n; i++) { // Iterates over the columns to remove

        for(int j = i + 1; j < n; j++){ // Iterates over the remaining rows
            float m = A[j][i] / A[i][i];
            for(int k = i+1; k < n; k++){ // iterates over the remaining columns
                A[j][k] -= m* A[i][k];
            }
            A[j][i] = m;
        }
    }



    LU decomposition = split_lu(A, n);

    return decomposition;
}

#endif