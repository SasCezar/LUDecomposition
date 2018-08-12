#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "ludutils.h"



int main(int argc, char *argv[])
{
    int rank, nprocs;     // for storing this process' rank, and the number of processes

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    Matrix original;
    float **matrix; 

    int n = atoi(argv[1]);  
    if(rank == 0){
        char path[255];
		int out = snprintf(path, 255, "/mnt/c/Users/sasce/Desktop/Matrices/matrix_%i.csv", n);
		original = read_csv(path);
        matrix = original.matrix;
    }else{
        matrix =  matrix_create(n, n);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int rem = n % nprocs;       // elements remaining after division among processes
    int sum = 0;                // Sum of counts. Used to calculate displacements
    
    int *sendcounts = malloc(sizeof(int)*nprocs); // array describing how many elements to send to each process
    int *displs = malloc(sizeof(int)*nprocs); // array describing the displacements where each segment begins

    // calculate send counts and displacements
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

    // print calculated send counts and displacements for each process
    if (0 == rank) {
        for (int i = 0; i < nprocs; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    // divide the data among processes as described by sendcounts and displs
    
    MPI_Scatterv(&matrix[0][0], sendcounts, displs, MPI_FLOAT, &rec_buf[0][0], sendcounts[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    // print what each process received
    printf("%d: ", rank);
    for (int i = 0; i < rows; i++) {
        for(int j = 0; j < n; j++){
            printf("%i\t", (int)rec_buf[i][j]);
        }
        printf("\n");
    }
    printf("\n");

     for (int i = 0; i < rows; i++) {
        for(int j = 0; j < n; j++){
            rec_buf[i][j] = -rec_buf[i][j];
        }
    }

    MPI_Gatherv(&rec_buf[0][0], sendcounts[rank], MPI_FLOAT, &matrix[0][0], sendcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if(0 == rank){
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++){
                printf("%i\t", (int)matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    }

    MPI_Finalize();

    free(sendcounts);
    free(displs);

    return 0;
}