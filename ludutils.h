#ifndef LUDUTILS_H
#define LUDUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float **initalize_matrix(int rows, int cols) {
    float **matrix = (float **)malloc(rows * sizeof(float *)); ;

    for(int i = 0; i < rows; i++)
    {
        matrix[i] = (float *)malloc(cols * sizeof(float));
        for(int j = 0; j < cols; j++)
        {
            matrix[i][j] = 0;
        }
    }

    return matrix;
}


typedef struct{
    float **matrix;
    int n;
} Matrix;

typedef struct{
    Matrix L;
    Matrix U;
} LU;

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

/**
 * Class that implements an invalid file path exception
*/
class not_valid_path : public std::runtime_error {
public:
	/**
     * Secondary constructor, takes a message
     * @param msg Message of error
	*/
	not_valid_path(const char *msg) :std::runtime_error(msg) {}
};

/**
 * Class that implements an invalid file
*/
class not_valid_file : public std::runtime_error {
public:
	/**
     * Secondary constructor, takes a message
     * @param msg Message of error
	*/
	not_valid_file(const char *msg) :std::runtime_error(msg) {}
};

char* concat(const char *s1, const char *s2)
{
    /**
     * Helper function to concatenate two strings
     * @param s1 Pointer to the first string
     * @param s2 Pointer to the second string  
     * @return *result Concatenated strings     
    */
    char *result = (char*)malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

void partialPivot(int n, double** a, double* b, int j){
	int i,k,m,rowx;
	double xfac, temp, temp1, amax;

	amax = (double) fabs(a[j][j]);
    m = j;
    for (i=j+1; i<n; i++){   /* Find the row with largest pivot */
    	xfac = (double) fabs(a[i][j]);
    	if(xfac > amax) {amax = xfac; m=i;}
    }

    if(m != j) {  /* Row interchanges */
    	rowx = rowx+1;
    	temp1 = b[j];
    	b[j]  = b[m];
    	b[m]  = temp1;
    	for(k=j; k<n; k++) {
    		temp = a[j][k];
    		a[j][k] = a[m][k];
    		a[m][k] = temp;
    	}
    }
}

Matrix read_csv(char* path){
    /**
     * Reads a csv files and creates a matrix
     * @param path Path to the csv files
     * @param *matrix Pointer to the matrix 
     * @throw not_valid_fie
     * @throw not_valid_path
    */

    FILE *fp;
    fp = fopen (path, "r");
    if (!fp) {
        throw not_valid_path(concat("Unable to find file:", path));
    }

    int idx = 0;
    int j = 0;
    char *buffer = NULL;
    size_t len = 0;
    ssize_t read;

    getline (&buffer, &len, fp);
    int n = atoi(buffer);
    float **matrix = (float **)malloc(n * sizeof(float *)); ;
    float cell;
    while ((read = getline (&buffer, &len, fp)) != -1) {
        matrix[idx] = (float *)malloc(n * sizeof(float));

        char *token = strtok(buffer, ","); // Split the row using the ';' char
        int j = 0;
        while (token != NULL) // Iterate over the tokens extracted by strtok function
        {
            sscanf(token, "%f", &cell);
            matrix[idx][j] = cell;
            j++;
            token = strtok(NULL, ",");
        }
        if (j != n) {
                throw not_valid_file("Invalid informations in file.\nThe number of matrix cols is not equal to the cols defined in file.");
        }
        idx++;
    }
    
    if (idx != n) {
            throw not_valid_file("Invalid informations in file.\nThe number of matrix rows is not equal to the rows defined in file.");
    }



    fclose (fp);

    Matrix result = {.matrix = matrix, .n = n};
    return result;
}


void print_matrix(Matrix mat){
    /**
     * Prints the matrix
     * @param matrix The matrix to print 
    */

    int n = mat.n;
    float **matrix = mat.matrix;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

Matrix matrix_multiplication(Matrix a,Matrix b){
    int n = a.n;
    float **ma = a.matrix;
    float **mb = b.matrix;
    float **mc = initalize_matrix(n, n);


    float sum = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            sum = 0;
            for(int k = 0; k < n; k++){
                sum = sum + ma[i][k] * mb[k][j];
            }
            mc[i][j] = sum;        
        }
    }
    Matrix res = {.matrix = mc, .n = n};
    return res;
}

Matrix matrix_difference(Matrix a, Matrix b){
    int n = a.n;
    float **ma = a.matrix;
    float **mb = b.matrix;
    float **mc = initalize_matrix(n, n);


    float sum = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            mc[i][j] = ma[i][j] - mb[i][j];
        }
    }
    Matrix res = {.matrix = mc, .n = n};
    return res;
}


double frobenius_norm(Matrix matrix){
    float **A = matrix.matrix;
    float sum = 0;
    
    for(int i = 0; i < matrix.n; i++)
    {
        for(int j = 0; j < matrix.n; j++)
        {
            sum += A[i][j] * A[i][j];
        }        
    }
    
    double norm = sqrt(sum);
    return norm;
}

double compute_error(Matrix a, Matrix b){
    Matrix difference = matrix_difference(a, b);
    double error = frobenius_norm(difference);
    return error;
}

Matrix duplicate_matrix(Matrix mat){
    int n = mat.n;

    float **matrix = (float **)malloc(n * sizeof(float *)); ;

    float sum = 0;
    for(int i = 0; i < n; i++)
    {
        matrix[i] = (float *)malloc(n * sizeof(float));
        for(int j = 0; j < n; j++)
        {
            matrix[i][j] = mat.matrix[i][j];      
        }
    }
    
    Matrix res = {.matrix = matrix, .n = n};
    return res;
}

#endif