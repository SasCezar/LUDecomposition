#include <iostream>
#include "ludutils.h"
#include "lud.h"


int main(int argc, char *argv[]){
	Matrix matrix = read_csv("test.csv");
	print_matrix(matrix);
	LU result = decompose_serial(matrix);
	printf("Lower\n");
	print_matrix(result.L);
	printf("Upper\n");
	print_matrix(result.U);

	Matrix recomposed = matrix_multiplication(result.L, result.U);
	printf("Recomposed\n");
	print_matrix(recomposed);

	free(matrix.matrix);

	return 0;
}