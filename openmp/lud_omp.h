#ifndef LUDOMP_H
#define LUDOMP_H

#include "ludutils.h"
#include <omp.h>


LU decompose_omp(Matrix matrix){
    int n = matrix.n;
    float **A = matrix.matrix;

    int numThreads = 16;

	omp_set_num_threads(numThreads);

    for(int i = 0; i < n; i++) {

        // pivoting()

        #pragma omp parallel for shared(A,n,i)
        for(int j = i + 1; j < n; j++) {
            // int tid = omp_get_thread_num();
            // printf("Hello from thread = %d\n", tid);
            float m = A[j][i] / A[i][i];
            for(int k = i + 1; k < n; k++) {
                A[j][k] -= m * A[i][k];
            }
            A[j][i] = m;
        }

    }

    LU decomposition = split_lu(A, n);
    return decomposition;

}

#endif