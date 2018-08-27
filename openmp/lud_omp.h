#ifndef LUDOMP_H
#define LUDOMP_H

#include "../ludutils.h"
#include <omp.h>


LU decompose_omp(Matrix matrix, int num_threads){
    int n = matrix.n;
    float **A = matrix.matrix;

	omp_set_num_threads(num_threads);

    int tid;
    for(int i = 0; i < n; i++) {
        #pragma omp parallel for shared(A,n,i)
        for(int j = i + 1; j < n; j++) {
            A[j][i] = A[j][i] / A[i][i];
        }
        #pragma omp parallel for shared(A,n,i)
        for(int j = i + 1; j < n; j++) {
            for(int k = i + 1; k < n; k++) {
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }

    LU decomposition = split_lu(A, n);
    return decomposition;

}

#endif