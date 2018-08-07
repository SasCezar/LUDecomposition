#ifndef LUDS_H
#define LUDS_H

#include "ludutils.h"

LU decompose_serial(Matrix mat){
    int n = mat.n;
    float **matrix = mat.matrix;

    Matrix ma = duplicate_matrix(mat);
    float **a = ma.matrix;

    for (int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++){
            float m = a[j][i] / a[i][i];
            for(int k = i+1; k < n; k++){
                a[j][k] -= m* a[i][k];
            }
            a[j][i] = m;
        }
    }

    LU decomposition = split_lu(a, n);

    return decomposition;
}

#endif