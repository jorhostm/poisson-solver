#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "poissonsolver.h"


static double phi(double x, double y){

	return 0.25*(x*x+y*y);
	
}
static double g0(double x, double y){

	return 0.0;

}

static double g1(double x, double y){

	return 12-12*x-12*y;

}
static double g2(double x, double y){

	return (6-12*x)*(3*y*y-2*y*y*y) + (3*x*x-2*x*x*x)*(6-12*y);

}

int main(int argc, char *argv[]){
	
	if (argc > 2){
		int n = atoi(argv[1]);



		time_t t0 = time(NULL);

		//double* r = solve_poisson(n, &phi, &g1, NM_FALSE,NM_FALSE,NM_FALSE,NM_FALSE);
		double* r = solve_poisson(n, &phi, &g1, NM_TRUE,NM_TRUE,NM_TRUE,NM_TRUE);

		double x = 0.0;
		double y = 0.0;
		double temp = get_value_at(r, n, x,y);
		printf("value at %f,%f is %f\n",x,y, temp);

		shift_solution(r, n, -1*temp);

		time_t t = time(NULL) - t0;
		printf("Time elapsed: %d seconds\n", t);

		print_solution_to_file(r, n, argv[2]);

		free(r);

		if(n<30){
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					printf("%f\n", r[i*n+j]);
				}
			printf("\n");
			}
		}
	}
	
	return 0;
}