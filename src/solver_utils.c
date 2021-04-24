/**
 * @file solver_utils.c
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This file implements helper functions used by the Poisson equation solver
 * @version 1.0
 * @date 2021-04-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <solver_utils.h>

/**
 * @brief Computes the theoretical optimal relaxation factor omega for the SOR iteration method.
 *	Works for Laplace and Poisson equations in a rectangular domain.
 *	In this case the domain is a n-by-n grid with x,y = [0, 1].
 *	The specral radius rho(A) = cos(PI*h), where h is the step size = 1/(n-1)
 * 
 * @param n The amount of discretized points along an axis.
 * @return double The theoretical optimal relaxation factor
 */
double get_omega(const unsigned int n){
	
	double rho = cos(M_PI/(n-1)); 					// Spectral radius of the coefficients-matrix A, rho(A) for Ax = b
	
	double omega = 2.0/(1.0+sqrt(1.0-pow(rho,2)));	// Omega as a function of the spectral radius rho(A)
	
	return omega;
}


void add_values(double** T, double (*func)(double x, double y), const unsigned int n, double* x_vals, double* y_vals, const unsigned int i_start, const unsigned int i_end, const unsigned int j_start, const unsigned int j_end){
	
	if (func == NULL)
	{
		return;
	}
	

	for(unsigned int i = i_start; i < i_end; i++){
		for(unsigned int j = j_start; j < j_end; j++){
			T[i][j] = func(x_vals[i], y_vals[j]);
		}
	}
}

void sor_iterate(double **T, double **b, const unsigned int n, const double omega, double *Tsum, double *dTsum, const unsigned int i_start, const unsigned int i_end, const unsigned int j_start, const unsigned int j_end){
	for(unsigned int i = i_start; i < i_end; i++){
		for(unsigned int j = j_start; j < j_end; j++){
			
			/* Compute the residual */
                double R = a0*T[i][j-1] + a1*T[i-1][j] + a2*T[i][j+1] + a3*T[i+1][j] - 4*T[i][j] - h2*b[i][j];
                /* Apply relaxation factor */
				double dT = omega*0.25*R;
                /* Add the residual to the solution */
				T[i][j] += dT;
			
				T_sum += fabs(T[i][j]);
				dT_sum += fabs(dT);

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