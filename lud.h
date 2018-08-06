#ifndef LUDS_H
#define LUDS_H

#include "ludutils.h"

LU decompose_serial(Matrix mat){
    int n = mat.n;
    float **matrix = mat.matrix;
    float **L = initalize_matrix(n, n);
    float **U = initalize_matrix(n, n);

    for (int i = 0; i < n; i++) {
 
        // Upper Triangular
        for (int k = i; k < n; k++) {
 
            // Summation of L(i, j) * U(j, k)
            int sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
 
            // Evaluating U(i, k)
            U[i][k] = matrix[i][k] - sum;
        }
 
        // Lower Triangular
        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal as 1
            else {
 
                // Summation of L(k, j) * U(j, i)
                int sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
 
                // Evaluating L(k, i)
                L[k][i] = (matrix[k][i] - sum) / U[i][i];
            }
        }
    }

    Matrix mL = { .matrix = L, .n = n};
    Matrix mU = { .matrix = U, .n = n};
    LU decomposition = { .L = mL, .U = mU};

    return decomposition;

}

#endif