/**
 * @file example.c
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This file works as an example on how to use the poisson-solver functions to create,
 *  solve and save solutions to Poisson Boundary Value Problems
 * @version 1.0
 * @date 2021-04-22
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <ctype.h>
#include <getopt.h>

#include <example.h>

#define MIN_N 10	// Arbitrary minimum number for discretization

int main(int argc, char **argv){
	
	// Default flags and values
	unsigned int n = MULTIGRID_MIN;	// Number of discretisation points for the n * n grid
	unsigned int use_multigrid = 0;	// Whether to use multigrids or not
	unsigned int save_data = 0;		// Whether to save the data to file
	unsigned int np = 5;			// Number of problems to solve
	double reltol = 1e-6;			// Relative tolerance

  	int c;
	// Process the options
	while ((c = getopt(argc, argv, "hmsn:r:p:")) != -1)
		switch (c)
		{
		case 'h':
			printf(help_msg);
			return EXIT_SUCCESS;
			break;
		case 'm':
			use_multigrid = 1;
			break;
		case 's':
			save_data = 1;
			break;
		case 'n':
			n = atoi(optarg);
			if(n < MIN_N){
					fprintf(stderr,"Invalid argument %s for n. Needs to be a positive integer >= %d\n", optarg, MIN_N);
					return EXIT_FAILURE;
				}
			break;
		case 'r':
			reltol = atof(optarg);
			if(reltol <= 0.0){
				fprintf(stderr,"Invalid argument %s for reltol. Needs to be a double > 0\n", optarg);
				return EXIT_FAILURE;
			}
			break;
		case 'p':
			np = atoi(optarg);
			if(np < 1 || np > 5){
					fprintf(stderr,"Invalid argument %s for np. Needs to be a positive integer > 0 and < 5\n", optarg);
					return EXIT_FAILURE;
				}
			break;
		case '?':
			if (optopt == 'n' || optopt == 'r')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,"Unknown option character `\\x%x'.\n", optopt);
			return EXIT_FAILURE;
		default:
			abort ();
		}
	
/**
 * @brief Testing the solver.
 *	Solve the poisson equation with mixed boundary conditions and g(x,y) = 0
 *	Neumann at x = 0 and y = 0, Dirichlet at x = 1 and y = 1 where phi(x,1) = 1 and phi(y,1) = 0
 */
	solver_params mixed = {"Mixed", n, &mix, NULL, save_data, use_multigrid, reltol, NM_X0 | NM_Y0};
/**
 * @brief Solve the poisson equation with only Dirichlet boundary conditions and g(x,y) = g1
 */
	solver_params dirichlet1 = {"Dirichlet1", n, &phi, &g1, save_data, use_multigrid, reltol, NM_NONE};

/**
 * @brief Solve the poisson equation with only Dirichlet boundary conditions and g(x,y) = g2
 */
	solver_params dirichlet2 = {"Dirichlet2", n, &phi, &g2, save_data, use_multigrid, reltol, NM_NONE};

/**
 * @brief Solve the poisson equation with only Neumann boundary conditions and g(x,y) = g1,
 * given phi(0,0) = 0
 */
	solver_params neumann1 = {"Neumann1", n, NULL, &g1, save_data, use_multigrid, reltol, NM_ONLY, {0,0,0}};
/**
 * @brief Solve the poisson equation with only Neumann boundary conditions and g(x,y) = g2,
 * given phi(0,0) = 0
 */
	solver_params neumann2 = {"Neumann2", n, NULL, &g2, save_data, use_multigrid, reltol, NM_ONLY, {0,0,0}};
	
/**
 * Solve each problem in their own separate thread
 */
	solver_params *problems[5] = {&mixed, &dirichlet1, &dirichlet2, &neumann1, &neumann2};
	pthread_t thread[5];
	
	for(int i = 0; i < np; i++){
		pthread_create( &thread[i], NULL, &solve, (void*) problems[i]);
	}
	for(int i = 0; i < np; i++){
		pthread_join( thread[i], NULL);
	}

	return EXIT_SUCCESS;
}