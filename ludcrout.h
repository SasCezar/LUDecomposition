#ifndef LUDC_H
#define LUDC_H

#include "ludutils.h"

LU decompose_crout(Matrix mat){
    int n = mat.n;
    float **A = mat.matrix;
    
    float **L = initalize_matrix(n, n);
    float **U = initalize_matrix(n, n);

    int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {
				printf("det(L) close to 0!\n Can't divide by 0...\n");
				exit(EXIT_FAILURE);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}

    Matrix mL = { .matrix = L, .n = n};
    Matrix mU = { .matrix = U, .n = n};
    LU decomposition = { .L = mL, .U = mU};

    return decomposition;
}

#endif