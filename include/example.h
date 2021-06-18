/**
 * @file example.h
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This header file contains example functions that can be used with the general Poisson solver API
 * @version 1.0
 * @date 2021-04-23
 * 
 * 
 */

#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <stdio.h>
#include <sys/time.h>

#include <poissonsolver.h>

#define MIN_N 10    // Arbitrary minimum number for discretization

const char *help_msg = " -h \t\t Display help message\n -n [int] \t Set descretisaton number\n -m \t\t Use multigrid\n -s \t\t Save solution to file\n -r [double]\t Relative tolerance (default: 1e-6)\n -p [int]\t Number of problems to solve (default: 5)\n";

/**
 * @brief Parameters to be passed to a thread running solve(). Used to create, solve, correct and save the solution of a BVP
 */
typedef struct solver_params {
    char *name;
    unsigned int n;
    double (*phi)(double x, double y);
    double (*g)(double x, double y);
    unsigned int save_data;
    unsigned int use_multigrid;
    double reltol;
    int nm_flags;
    double value_at_xyz[3];
} solver_params;

/**
 * @brief Example of a Dirichlet boundary function. 
 * Used for the Mixed boundaries example, where phi(x,1) = 1 and phi(1,y) = 0, Neumann for x = 0 and y = 0
 * 
 * @param x The x-coordinate
 * @param y The y-coordinate
 * @return double The result of phi at boundary coordinate (x,y)
 */
double mix(double x, double y){
    if(y == 1.0){
        return 1.0;
    }

    return 0.0;
}

/**
 * @brief Example of a Dirichlet boundary function. 
 * Used for the Mixed boundaries example, where phi(x,1) = 1 and phi(1,y) = 0, Neumann for x = 0 and y = 0
 * 
 * @param x The x-coordinate
 * @param y The y-coordinate
 * @return double The result of phi at boundary coordinate (x,y)
 */
double phi(double x, double y){

	return 0.25*(x*x+y*y);
	
}

/**
 * @brief Example of a function g(x,y) on the right side of the equation
 * 
 * @param x The x-coordinate
 * @param y The y-coordinate
 * @return double The result of g(x,y)
 */
double g1(double x, double y){

	return 12.0-12.0*x-12.0*y;

}

/**
 * @brief Example of a function g(x,y) on the right side of the equation
 * 
 * @param x The x-coordinate
 * @param y The y-coordinate
 * @return double The result of g(x,y)
 */
double g2(double x, double y){

	return (6.0-12.0*x)*(3.0*y*y-2.0*y*y*y) + (3.0*x*x-2.0*x*x*x)*(6.0-12.0*y);

}

/**
 * @brief An example function that create, solve, correct and save a BVP using function calls to the solver API
 * 
 * @param parameters The parameters that describe the BVP, how it should be solved and saved
 */
void solve( void *parameters){
    
    solver_params *params = (solver_params*) parameters;
    
    char *name = params->name;
	
    /* Create a poisson boundary value problem */
	bvp_t bvp = bvp_create(params->n, params->phi, params->g, params->nm_flags);
	
    /* Start timer */
    struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
	
    /* Solve the BVP */
	int i = bvp_solve(bvp, params->use_multigrid, params->reltol);
	
    /* End timer */
	clock_gettime(CLOCK_MONOTONIC, &end);
	double secs = (double) (end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1000000000;

    /* If we only have neumann boundaries, then we need to correct the solution */
    if (params->nm_flags == NM_ONLY)
    {
        double *values = params->value_at_xyz;
        double x = values[0];
        double y = values[1];
        double z = values[2];
        
        double delta_t = bvp_get_value_at(bvp,x,y);
	    bvp_shift_solution(bvp, z-delta_t);
    }

	printf("%s finished after %f seconds using %d iterations\n", name, secs, i);
	
	if (params->save_data){
		bvp_print_solution_to_file(bvp, name);
		bvp_create_gnuplot_data(bvp, name);
	};

	bvp_destroy(bvp);

}
#endif