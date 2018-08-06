#ifndef LUDMP_H
#define LUDMP_H

#include "ludutils.h"


void lup_od_omp(Matrix matrix){
    int i,k;

    int n = matrix.n;
    float **A = matrix.matrix;

    for(k = 0; k < n - 1; ++k) {
        // For the vectorization
        for(i = k + 1; i < n; i++) {
            A[i][k] /= A[k][k];
        }

        #pragma omp parallel for shared(A,n,k) private(i) schedule(static, 64)
        for(i = k + 1; i < n; i++) {
            int j;
            const double Aik = A[i][k];
            for(j = k + 1; j < n; j++) {
                A[i][j] -= Aik * A[k][j];
            }
        }
    }
}

#endif