#ifndef MATERR_H
#define MATERR_H

#include <stdlib.h>

int equals(float **a, float **b, float threshold){
    /**
     * Compares two matrixes with a given error threshold
     * @param **a The first matrix
     * @param **b The second matrix
     * @param threshold The maximum difference acceptable between two elements of the matrixes
    */
    int rowsa = sizeof(a) / sizeof(a[0])+1;
    int rowsb = sizeof(a) / sizeof(a[0])+1;

    int colsa = sizeof(a[0])/sizeof(a[0][0]);
    int colsb = sizeof(b[0])/sizeof(b[0][0]);
    

    if (rowsa != rowsb || colsa != colsb) {
        return 0;
    }
    

    for(int i = 0; i < rowsa; i++)
    {
        for(int j = 0; j < colsa; j++)
        {
            float diff = abs(a[i][j] - b[i][j]);
            
            if (diff > threshold) {
                return 0;
            }
            
        }
    }
    return 1;
}

#endif