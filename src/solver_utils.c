#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <solver_utils.h>

/* 
	Computes the theoretical optimal relaxation factor omega for the SOR iteration method.
	Works for Laplace and Poisson equations in a rectangular domain.
	In this case the domain is a n-by-n grid with x,y = [0, 1].
*/
double _omega(const unsigned int n){
	
	double rho = cos(M_PI/(n-1));
	double omega = 2.0/(1.0+sqrt(1.0-pow(rho,2)));
	
	return omega;
}


void _init_values(double** T, double (*func)(double x, double y), const unsigned int n, double* x_vals, double* y_vals, int i0, int i1, int j0, int j1){
	for(int i = i0; i < i1; i++){
		for(int j = j0; j < j1; j++){
			T[i][j] = func(x_vals[i], y_vals[j]);
		}
	}
}

void _copy_neumann_border(double** T, const unsigned int n, int dest_x0, int dest_x1, int dest_y0, int dest_y1, int src_x0, int src_x1, int src_y0, int src_y1){

	int dest_dx = dest_x1 - dest_x0 + 1;
	int dest_dy = dest_y1 - dest_y0 + 1;
	int src_dx = src_x1 - src_x0 + 1;
	int src_dy = src_y1 - src_y0 + 1;

	if (dest_dx != src_dx || dest_dy != src_dy){
		printf("Error in copy border: Dimensions do not match\n" );
		exit(1);
	}

	for(int i = 0; i < dest_dx; i++){
		for(int j = 0; j < dest_dy; j++){
			T[dest_x0 + i][dest_y0 + j] = T[src_x0 + i][src_y0 + j];
		}
	}
}

double* create_linear_array(const double from, const double to, const unsigned int count){

	double h = (to - from)/((double) count - 1);

	double *arr = malloc( count * sizeof(double) );

	arr[0] = from;
	arr[count-1] = to;


	for(int i = 1; i < (count-1); i++){
		arr[i] = arr[i-1] + h;
	}

	return arr;

}