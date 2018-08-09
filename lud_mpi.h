#ifndef LUDMPI_H
#define LUDMPI_H

#include "ludutils.h"
#include <mpi.h>

inline int map(int i, int n) {
  return i % nprocs;
}

LU decompose_mpi(Matrix matrix, int argc, char *argv[]){
    int n = matrix.n;
    float **A = matrix.matrix;
    
    int rank, nproc;
    int i=0, j=0, k=0;
    MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);	/* get number of processes */


	for(j =0; j < n-1; j++) {
        if(map(j, nproc) == rank)
        {
        	for(i = j+1; i < n; i++) {
        		A[i][j] = A[i][j]/A[j][j];
        	}
        }
        MPI_Bcast (&A[j][j],n-j,MPI_DOUBLE,map(j, nproc),MPI_COMM_WORLD);

        		for(k = j+1; k < n; k++)
        		{
        			if(map(j, nproc) == rank)
        			{
        				for(i = j+1; i < n; i++) {
        					A[k][i]= A[k][i] - (A[k][j] * A[j][i]);
        				}
        			}
        		}
        }

    // for (int i = 0; i < n; i++) {
    //     for(int j = i + 1; j < n; j++){
    //         float m = A[j][i] / A[i][i];
    //         for(int k = i+1; k < n; k++){
    //             A[j][k] -= m* A[i][k];
    //         }
    //         A[j][i] = m;
    //     }
    // }

    MPI_Finalize();

    LU decomposition = split_lu(A, n);

    return decomposition;
}

#endif