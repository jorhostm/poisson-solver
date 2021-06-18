/**
 * @file solver_utils.c
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This file implements helper functions used by the Poisson equation solver
 * @version 1.0
 * @date 2021-04-23
 * 
 * 
 */

#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include <solver_utils.h>

/**
 * @brief Computes the theoretical optimal relaxation factor omega for the SOR iteration method using Chebyshev acceleration.
 *	Works for Laplace and Poisson equations in a rectangular domain.
 *	In this case the domain is a n-by-n grid with x,y = [0, 1].
 *	The specral radius the Jacobi matrix rho(C_J) = cos(PI*h), where h is the step size = 1/(n-1)
 * 
 * @param n The amount of discretized points along an axis.
 * @return double The theoretical optimal relaxation factor
 */
double get_next_omega(double omega, const unsigned int n){
	
	double rho = cosl(M_PI/(n-1)); 		// Spectral radius of the Jacobi iteration matrix, rho(C_J)
	
	omega = 1.0/(1-rho*rho*omega/4);	// Chebyshev acceleration
	
	return omega;
}

/**
 * @brief Adds the values to a 2D array by evaluating the given function for the corresponding x- and y-values.
 * 
 * @param T 		The 2D array
 * @param func 		The given function to evaluate
 * @param n 		The amount of discretized points along an axis.
 * @param x_vals 	The array of x-values
 * @param y_vals 	The array of y-values
 * @param i_start 	Start index for x-values
 * @param i_end 	End index for x-values
 * @param j_start 	Start index for y-values
 * @param j_end 	End index for y-values
 */
void add_values(double** T, double (*func)(double x, double y), const unsigned int n, double* x_vals, double* y_vals, const unsigned int i_start, const unsigned int i_end, const unsigned int j_start, const unsigned int j_end){
	
	if (func == NULL)
	{
		return;
	}
	

	for(int i = i_start; i < i_end; i++){
		for(int j = j_start; j < j_end; j++){
			T[i][j] = func(x_vals[i], y_vals[j]);
		}
	}
}

/**
 * @brief Create a linear array from one value to another
 * 
 * @param from The start value
 * @param to 	The end value
 * @param count The total number of values in the array
 * @return double* pointer to the linear array
 */
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