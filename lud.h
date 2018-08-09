#ifndef LUDS_H
#define LUDS_H

#include "ludutils.h"

LU decompose_serial(Matrix matrix){
    int n = matrix.n;
    float **A = matrix.matrix;

    for (int i = 0; i < n; i++) {
        for(int j = i + 1; j < n; j++){
            float m = A[j][i] / A[i][i];
            for(int k = i+1; k < n; k++){
                A[j][k] -= m* A[i][k];
            }
            A[j][i] = m;
        }
    }

    LU decomposition = split_lu(A, n);

    return decomposition;
}

#endif