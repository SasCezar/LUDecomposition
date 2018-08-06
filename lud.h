#ifndef LUDS_H
#define LUDS_H

#include "ludutils.h"

LU split_lu(float **a, int n){
    float **L = initalize_matrix(n, n);
    float **U = initalize_matrix(n, n);

    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i < j)
            {
                U[i][j] = a[i][j];
            }
            if (i > j){
                L[i][j] = a[i][j];
            }
            if (i == j){
                L[i][j] = 1;
                U[i][j] = a[i][j];
            }
        }
    }


    Matrix mL = { .matrix = L, .n = n};
    Matrix mU = { .matrix = U, .n = n};
    LU decomposition = { .L = mL, .U = mU};
    return decomposition;
}

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