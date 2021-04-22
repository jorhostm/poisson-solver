#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <solver_utils.h>

/* 
	Computes the theoretical optimal relaxation factor omega for the SOR iteration method.
	Works for Laplace and Poisson equations in a rectangular domain.
	In this case the domain is a n-by-n grid with x,y = [0, 1].
	rho(A) = cos(PI*h), h = 1/(n-1)
*/
double _omega(const unsigned int n){
	
	double rho = cos(M_PI/(n-1)); 					// Spectral radius of the coefficients-matrix A, rho(A) for Ax = b
	
	double omega = 2.0/(1.0+sqrt(1.0-pow(rho,2)));	// Omega as a function of the spectral radius rho(A)
	
	return omega;
}


void _init_values(double** T, double (*func)(double x, double y), const unsigned int n, double* x_vals, double* y_vals, int i0, int i1, int j0, int j1){
	
	if (func == NULL)
	{
		return;
	}
	

	for(int i = i0; i < i1; i++){
		for(int j = j0; j < j1; j++){
			T[i][j] = func(x_vals[i], y_vals[j]);
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

int is_admissable(const unsigned int n, const unsigned int multigrid_min){
	
	if ( n <= multigrid_min)
	{
		return 1;
	}
	else if ( (n-1)%2 )
	{
		return 0;
	}
	else
	{
		return is_admissable(n/2 + 1, multigrid_min);
	}
	

}

unsigned int find_admissable(unsigned int n, const unsigned int multigrid_min){

	unsigned int result = n;

	while ( ! is_admissable(result, multigrid_min))
	{
		result++;
	}
	
	return result;
}