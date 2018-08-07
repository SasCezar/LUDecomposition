#ifndef LUDMP_H
#define LUDMP_H

#include "ludutils.h"
#include <omp.h>


LU decompose_omp(Matrix matrix){
    int i, j;

    int n = matrix.n;
    float **A = matrix.matrix;

    int numThreads = 8;

	omp_set_num_threads(numThreads);

    for(i = 0; i < n - 1; ++i) {
        // For the vectoriser
        for(j = i + 1; j < n; j++) {
            A[j][i] /= A[i][i];
        }

        #pragma omp parallel for shared(A,n,i) private(j)
        for(j = i + 1; j < n; j++) {
            int k;
            int tid = omp_get_thread_num();
            printf("Hello from thread = %d\n", tid);
            const double Aji = A[j][i];
            for(k = i + 1; k < n; k++) {
                A[j][k] -= Aji * A[i][k];
            }
        }
    }

    LU decomposition = split_lu(A, n);
    return decomposition;

}

#endif