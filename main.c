#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "include/poissonsolver.h"
#include "include/functions.h"


int main(int argc, char *argv[]){
	
	if(argc == 1){
		fprintf(stderr, "Invalid number of arguments");
		return EXIT_FAILURE;
	}

	// Get the number of discretised points along an axis
	int n = atoi(argv[1]);
	if(n <= 0){
		fprintf(stderr,"Invalid argument %s", argv[1]);
		return EXIT_FAILURE;
	}

	time_t t0 = time(NULL);
	
	/* 
		Pre-allocate the solution array for an n*n grid
		Reusable to save memory
	*/
	double *result = malloc(n*n*sizeof(double));
	
	/*	
		Testing of the solver
		Solve the poisson equation with mixed boundary condition and g(x,y) = 0
		Neumann at x = 0 and y = 0, Dirichlet phi(x,1) = 1 and phi(1,y) = 0
		Then save the solution to file
	*/
	solve_poisson(result,n, &mix, &zero, NM_TRUE,NM_FALSE,NM_TRUE,NM_FALSE);
	print_solution_to_file(result, n, "solutions/Mixed.txt");


	/*
		Solve the poisson equation with Dirichlet boundary condition and g(x,y) = g1
		Then save the solution to file
	*/
	solve_poisson(result,n, &phi, &g1, NM_FALSE,NM_FALSE,NM_FALSE,NM_FALSE);
	print_solution_to_file(result, n, "solutions/Dirichlet1.txt");
	
	
	/*
		Solve the poisson equation with Dirichlet boundary condition and g(x,y) = g2
		Then save the solution to file
	*/
	solve_poisson(result, n, &phi, &g2, NM_FALSE,NM_FALSE,NM_FALSE,NM_FALSE);
	print_solution_to_file(result, n, "solutions/Dirichlet2.txt");
	

	/*
		Solve the poisson equation with Neumann boundary condition and g(x,y) = g1
		Correct the solution so that phi(0,0) = 0
		Then save the solution to file
	*/
	double delta_t;
	solve_poisson(result,n, &zero, &g1, NM_TRUE,NM_TRUE,NM_TRUE,NM_TRUE);
	delta_t = get_value_at(result, n, 0, 0);
	shift_solution(result,n,-delta_t);
	print_solution_to_file(result, n, "solutions/Neumann1.txt");
	

	/*
		Solve the poisson equation with Neumann boundary condition and g(x,y) = g2
		Correct the solution so that phi(0,0) = 0
		Then save the solution to file
	*/
	solve_poisson(result,n, &zero, &g2, NM_TRUE,NM_TRUE,NM_TRUE,NM_TRUE);
	delta_t = get_value_at(result, n, 0, 0);
	shift_solution(result,n,-delta_t);
	print_solution_to_file(result, n, "solutions/Neumann2.txt");
	
	free(result);
	
	
	time_t t = time(NULL) - t0;
	printf("Time elapsed: %d seconds\n", t);

	return EXIT_SUCCESS;
}