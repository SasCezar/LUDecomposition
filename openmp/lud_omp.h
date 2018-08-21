#ifndef LUDOMP_H
#define LUDOMP_H

#include "../ludutils.h"
#include <omp.h>


LU decompose_omp(Matrix matrix){
    int n = matrix.n;
    float **A = matrix.matrix;

    int num_threads = 16;

	omp_set_num_threads(num_threads);

    int tid;
    for(int i = 0; i < n; i++) {
        #pragma omp parallel for shared(A,n,i)
        for(int j = i + 1; j < n; j++) {
            A[j][i] = A[j][i] / A[i][i];
        }
        #pragma omp parallel for shared(A,n,i)
        for(int j = i + 1; j < n; j++) {
            tid = omp_get_thread_num();
            for(int k = i + 1; k < n; k++) {
                printf("Thread=%d did i=%d - j=%d - k=%d\n",tid, i, j, k);
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }

    LU decomposition = split_lu(A, n);
    return decomposition;

}

#endif