#ifndef LUDMPI_H
#define LUDMPI_H

#include "ludutils.h"
#include "mpi.h"

int map(int i, int n) {
  return i % n;
}

LU decompose_mpi(Matrix matrix){
    int n = matrix.n;
    float **A = matrix.matrix;
    
    int rank, nprocs;
    int i=0, j=0, k=0;
    MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);	/* get number of processes */


	for(j =0; j < n-1; j++) {
        if(map(j, nprocs) == rank)
        {
        	for(i = j+1; i < n; i++) {
        		A[i][j] = A[i][j]/A[j][j];
        	}
        }
        MPI_Bcast (&A[j][j],n-j,MPI_DOUBLE,map(j, nprocs),MPI_COMM_WORLD);

        for(k = j+1; k < n; k++){
        	if(map(j, nprocs) == rank){
        		for(i = j+1; i < n; i++) {
        			A[k][i]= A[k][i] - (A[k][j] * A[j][i]);
				}
  			}
        }
    }

    MPI_Finalize();

    LU decomposition = split_lu(A, n);

    return decomposition;
}

#endif