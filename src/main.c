#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include <poissonsolver.h>
#include <functions.h>

static int n;

static void solve0();
static void solve1();
static void solve2();
static void solve3();
static void solve4();

int main(int argc, char *argv[]){
	
	if(argc == 1){
		fprintf(stderr, "Invalid number of arguments");
		return EXIT_FAILURE;
	}

	// Get the number of discretised points along an axis
	n = atoi(argv[1]);
	if(n <= 0){
		fprintf(stderr,"Invalid argument %s", argv[1]);
		return EXIT_FAILURE;
	}

	time_t t0 = time(NULL);
	
	void *solve[5] = {&solve0, &solve1, &solve2, &solve3, &solve4};
	pthread_t thread[5];
	
	for(int i = 0; i < 5; i++){
		pthread_create( &thread[i], NULL, solve[i], NULL);
	}
	for(int i = 0; i < 5; i++){
		pthread_join( thread[i], NULL);
	}
	
	time_t t = time(NULL) - t0;
	printf("Time elapsed: %d seconds\n", t);

	return EXIT_SUCCESS;
}

/*	
	Testing the solver
	Solve the poisson equation with mixed boundary conditions and g(x,y) = 0
	Neumann at x = 0 and y = 0, Dirichlet phi(x,1) = 1 and phi(1,y) = 0
	Then save the solution to file
*/
static void solve0(){


	bvp_t *bvp;
	bvp_create(&bvp, n, &mix, &zero, NM_TRUE,NM_FALSE,NM_TRUE,NM_FALSE);
	solve_poisson_bvp(bvp);
	print_solution_to_file(bvp,"./solutions/Mixed.txt");
	bvp_free(bvp);
}

/*
	Solve the poisson equation with Dirichlet boundary conditions and g(x,y) = g1
	Then save the solution to file
*/
static void solve1(){


	bvp_t *bvp;
	bvp_create(&bvp, n, &phi, &g1, NM_FALSE,NM_FALSE,NM_FALSE,NM_FALSE);
	solve_poisson_bvp(bvp);
	print_solution_to_file(bvp,"solutions/Dirichlet1.txt");
	bvp_free(bvp);
}

/*
	Solve the poisson equation with Dirichlet boundary conditions and g(x,y) = g2
	Then save the solution to file
*/
static void solve2(){
	

	bvp_t *bvp;
	bvp_create(&bvp, n, &phi, &g2, NM_FALSE,NM_FALSE,NM_FALSE,NM_FALSE);
	solve_poisson_bvp(bvp);
	print_solution_to_file(bvp,"solutions/Dirichlet2.txt");
	bvp_free(bvp);
}

/*
	Solve the poisson equation with Neumann boundary conditions and g(x,y) = g1
	Correct the solution so that phi(0,0) = 0
	Then save the solution to file
*/
static void solve3(){

	bvp_t *bvp;
	bvp_create(&bvp, n, &zero, &g1, NM_TRUE,NM_TRUE,NM_TRUE,NM_TRUE);
	solve_poisson_bvp(bvp);
	double delta_t = get_value_at(bvp,0,0);
	shift_solution(bvp, -1*delta_t);
	print_solution_to_file(bvp,"solutions/Neumann1.txt");
	bvp_free(bvp);
}

/*
	Solve the poisson equation with Neumann boundary conditions and g(x,y) = g2
	Correct the solution so that phi(0,0) = 0
	Then save the solution to file
*/
static void solve4(){
	
	bvp_t *bvp;
	bvp_create(&bvp, n, &zero, &g2, NM_TRUE,NM_TRUE,NM_TRUE,NM_TRUE);
	solve_poisson_bvp(bvp);
	double delta_t = get_value_at(bvp,0,0);
	shift_solution(bvp, -1*delta_t);
	print_solution_to_file(bvp,"solutions/Neumann2.txt");
	bvp_free(bvp);
}