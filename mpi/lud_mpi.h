#ifndef LUDMPI_H
#define LUDMPI_H

#include "../ludutils.h"
#include "mpi.h"

LU decompose_mpi(Matrix matrix){
    int rank, nprocs;     // for storing this process' rank, and the number of processes

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   
    float **A;
    int n;
    if(rank == 0){
        A = matrix.matrix;
        n = matrix.n;
    }else{
        n = 2000;
        A = matrix_create(n, n);
    }
    
    
    printf("I'm here %i\n", rank);

    int rem = n % nprocs;       // Elements remaining after division among processes
    int sum = 0;                // Sum of counts. Used to calculate displacements

    int r_b_size = ceil(n / nprocs);
    float **rec_buf = matrix_create(n, n);

    int *sendcounts = malloc(sizeof(int) * nprocs); // Array describing how many elements to send to each process
    int *displs = malloc(sizeof(int) * nprocs); // Array describing the displacements where each segment begins
    

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

    // Divide the data among processes as described by sendcounts and displs
    MPI_Scatterv(&A[0][0], sendcounts, displs, MPI_FLOAT, &rec_buf[0][0], sendcounts[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);

    int rows = sendcounts[rank] / n ;
    int mystart = displs[rank] / n;
    int myend = mystart + rows; 
    
    if(rank == 0){
        printf("i\tj\tk\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = mystart; i < myend; i++) { // Iterates over the columns to remove
        for(int j = i + 1; j < rows; j++){ // Iterates over the remaining rows
            float m = rec_buf[j][i] / rec_buf[i][i];
            for(int k = i+1; k < n; k++){ // Iterates over the remaining columns
                rec_buf[j][k] -= m* rec_buf[i][k];
            }
            rec_buf[j][i] = m;
        }
    }

    MPI_Gatherv(&rec_buf[0][0], sendcounts[rank], MPI_FLOAT, &B[0][0], sendcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    LU decomposition = split_lu(A, n);

    return decomposition;
}

#endif