/**
 * @file poissonsolver.c
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This file constitute the implementation of the general Poisson equation solver API
 * @version 1.0
 * @date 2021-04-22
 * 
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#include <poissonsolver.h>
#include <solver_utils.h>


/**
 * @brief Iterate using Succesive Over-Relaxation for the Red-Black-Gauss-Seidel method
 * 
 * @param bvp The Boundary Value Problem
 * @param reltol The relative tolerence for the relative residual used to determine convergence. 
 * 				Smaller => greater accuracy ,but more iterations
 * @return int The number of iterations used to reach convergence
 */
static int red_black_sor(bvp_t restrict bvp, const double reltol){

	const unsigned int n = bvp->n;
	double **T = bvp->result;
	double **b = bvp->b;
	double h = 1.0/( n - 1.0);
	double h2 = h*h;
	double omega = 1.0;

	const unsigned int i_start 	= bvp->nm_flags & NM_X0 ? 0 : 1;
	const unsigned int i_end 	= bvp->nm_flags & NM_X1 ? n : n-1;
	const unsigned int j_start 	= bvp->nm_flags & NM_Y0 ? 0 : 1;
	const unsigned int j_end 	= bvp->nm_flags & NM_Y1 ? n : n-1;

	double rel_res;
	int iterations = 0;
	 
	do {

	 	double T_sum = 0.0;
		double dT_sum = 0.0;

		int isw = 0;
		
		/* Odd-even ordering */
		for (int pass = 1; pass <= 2; pass++){

			int jsw = isw;
			
			for (int i = i_start; i < i_end; i++){
				
				int a1 = 1;	//	T(i-1,j)
				int a3 = 1;	// 	T(i+1,j)
				
				/* Check if there is a Neumann boundary at x = 0 */
				if(i == 0){
					a1 = 0;
					a3 = 2;
				}
				/* Check if there is a Neumann boundary at x = 1 */
				else if(i == n-1){
					a1 = 2;
					a3 = 0;
				}

				for (int j = j_start + jsw; j < j_end; j+=2){

					int a0 = 1;	// 	T(i,j-1)
					int a2 = 1; // 	T(i,j+1)

					/* Check if there is a Neumann boundary at y = 0 */
					if(j == 0){
						a0 = 0;
						a2 = 2;
					}
					/* Check if there is a Neumann boundary at y = 1 */
					else if(j == n-1){
						a0 = 2;
						a2 = 0;
					}
					
					/* Compute the residual */
					double R = a0*T[i][j-1] + a1*T[i-1][j] + a2*T[i][j+1] + a3*T[i+1][j] - 4*T[i][j] - h2*b[i][j];
					/* Apply relaxation factor */
					double dT = omega*0.25*R;
					/* Add the residual to the solution */
					T[i][j] += dT;
				
					T_sum += fabs(T[i][j]);
					dT_sum += fabs(dT);

				}

				jsw = 1 - jsw;

			}

			isw = 1 - isw;
			omega = get_next_omega(omega,n);

		}

		rel_res = dT_sum/T_sum;	
		iterations++;
    
    } while (rel_res >= reltol);

    return iterations;
}

/**
 * @brief Solve the given Poisson Boundary Value Problem
 * 
 * @param bvp The Boundary Value Problem
 * @param use_multigrid Whether or not to us multigrid preconditioning
 * @param reltol The relative tolerence for the relative residual used to determine convergence. 
 * 				Smaller => greater accuracy ,but more iterations
 * @return int The number of iterations used to reach convergence
 */
int bvp_solve(bvp_t restrict bvp, const unsigned int use_multigrid, const double reltol){
	int n = bvp->n;
	int iterations = 0;
	
	if(n >= MULTIGRID_MIN && use_multigrid){
		bvp_t coarse_bvp = bvp_create(n/2, bvp->phi,bvp->g,bvp->nm_flags);;
		iterations += bvp_solve(coarse_bvp, use_multigrid, reltol);
		bvp_copy(coarse_bvp, bvp);
		bvp_destroy(coarse_bvp);
	}
	
	iterations += red_black_sor(bvp, reltol);

	return iterations;
}


/**
 * @brief Get the value at phi(x,y) using bilinear interpolation
 * 
 * @param bvp The Boundary Value Problem
 * @param x The x-coordinate
 * @param y The y-coordinate
 * @return double 
 */
double bvp_get_value_at(bvp_t restrict bvp, const double x, const double y){

	double **result = bvp->result;
	int n = bvp->n;
	double h = 1.0/(n-1.0);

	// Domain check
	if (x < 0 || x > 1 || y < 0 || y > 1){
		return NAN;
	}

	double *x_values = bvp->x_val;
	double *y_values = bvp->y_val;

	unsigned int i_1 = x*(n-1);
	unsigned int i_2 = ceil(x*(n-1));
	unsigned int j_1 = y*(n-1);
	unsigned int j_2 = ceil(y*(n-1));

	double t11 = result[i_1][j_1];
	double t21 = result[i_2][j_1];
	double t12 = result[i_1][j_2];
	double t22 = result[i_2][j_2];

	double x_norm = (x-x_values[i_1])/h;
	double y_norm = (y-y_values[j_1])/h;

	double t_xy = t11*(2.0 - x_norm - y_norm) + t21*(1.0 + x_norm - y_norm) + t12*(1.0 - x_norm + y_norm) + t22*(x_norm + y_norm);

	return 0.25*t_xy;
}

/**
 * @brief Shifts the entire solution by delta_t. Used to correct a solution with only Neumann boundaries
 * 
 * @param bvp The Boundary Value Problem
 * @param delta_t The value to shift the solution by
 */
bvp_t bvp_shift_solution(bvp_t restrict bvp, const double delta_t){

	if (bvp != NULL){

		double *result = bvp->result[0];
		int n = bvp->n;
		int n2 = n*n;
		
		for (int i = 0; i < n2; i++){

			result[i] += delta_t;

		}
	}

	return bvp;
}

/**
 * @brief Create and initialize a Poisson Boundary Value Problem
 * 
 * @param n The number of discretised points along an axis
 * @param phi The function phi(x,y) used to prescribe Dirichlet boundaries. If NULL, is threated as phi(x,y) = 0
 * @param g The function g(x,y) on the right side of the Poisson equation. If NULL, is threated as g(x,y) = 0
 * @param nm_flags Flags that is used to specify Neumann boundaries, if any. If there isn't any, set it as 0
 * @return bvp_t The Boundary Value Problem with an initial n-by-n grid with prescribed Dirichlet boundaries, if any
 */
bvp_t bvp_create(unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int nm_flags){
	
	bvp_t bvp = malloc(sizeof(struct bvp_t));

	bvp->n = n;
	bvp->phi = phi;
	bvp->g = g;
	bvp->nm_flags = nm_flags;

	double *x = create_linear_array(0.0, 1.0, n);
	double *y = x;

	bvp->x_val = x;
	bvp->y_val = y;

	// Initialize the RHS
	double *rhs = calloc(n*n,sizeof(double));
	bvp->b = calloc(n, sizeof(double*));

	for (int i = 0; i < n; i++)
	{
		bvp->b[i] = rhs + i*n;
	}
	add_values(bvp->b,g,n,x,y,0,n,0,n);

	// Initialize the 2D solution array, including "ghost cells" on the sides, to prevent out-of-bound access
	double *datapoints = calloc((n+2)*n,sizeof(double));
	
	bvp->result = calloc(n+2, sizeof(double*));
	for (int i = 0; i < n+2; i++)
	{
		bvp->result[i] = datapoints + i*n;
	}
	(bvp->result)++;
	
	// Initialize Dirichlet boundaries
	if(!(bvp->nm_flags & NM_X0)) add_values(bvp->result,phi,n,x,y,0,n,0,1);
	if(!(bvp->nm_flags & NM_X1)) add_values(bvp->result,phi,n,x,y,0,n,n-1,n);
	if(!(bvp->nm_flags & NM_Y0)) add_values(bvp->result,phi,n,x,y,0,1,0,n);
	if(!(bvp->nm_flags & NM_Y1)) add_values(bvp->result,phi,n,x,y,n-1,n,0,n);
	
	return bvp;
}

/**
 * @brief Copy result of a bvp over to another bvp. If destination BVP is NULL, create a new one. Uses bilinear interpolation
 * 
 * @param src_bvp The source BVP solution, non-null
 * @param dest_bvp The destination BVP solution, NULL if a new one is to be created
 * @return ibvp_t The destination BVP if successful, NULL otherwise.
 */
bvp_t bvp_copy(const bvp_t restrict src_bvp, bvp_t restrict dest_bvp){
	
	if (src_bvp == NULL)
		return NULL;
	
	else if (dest_bvp == NULL)
		dest_bvp = bvp_create(src_bvp->n , src_bvp->phi , src_bvp->g , src_bvp->nm_flags);

	
	const int dest_n = dest_bvp->n;
	double **dest_r = dest_bvp->result;

	for(int i = 0; i < dest_n; i++){
		
		double x = dest_bvp->x_val[i];
		
		for(int j = 0; j < dest_n; j++){
			
			double y = dest_bvp->y_val[j];
			dest_r[i][j] = bvp_get_value_at(src_bvp, x, y);
		
		}
	
	}

	return dest_bvp;
}

/**
 * @brief Free up memory occupied by the Boundary Value Problem
 * 
 * @param bvp The Boundary Value Problem
 */
void bvp_destroy(bvp_t restrict bvp){
	--(bvp->result);
	free((bvp->result[0]));
	free(bvp->result);

	free(bvp->x_val);
	free(bvp->b[0]);
	free(bvp->b);

	free(bvp);
}

/**
 * @brief Saves the solution to file in the following format:
 * n n
 * phi(0,0)
 * phi(0,1*h)
 * .
 * .
 * phi(0,(n-1)*h)
 * phi(1*h,0)
 * .
 * .
 * phi(1,1)
 * @param bvp The Boundary Value Problem
 * @param filename The full filename. Can include path
 */
void bvp_print_solution_to_file(bvp_t restrict bvp, const char *name){
	int n = bvp->n;
	double **result = bvp->result;

	char filename[50];
	sprintf(filename,"solutions/opengl/%s.dat", name);

	FILE *f;
    f = fopen(filename, "w");

    if (f == NULL)
	{
    	printf("Error opening file!\n");
    	exit(1);
	}

    fprintf(f, "%d %d", n, n);
    for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                
            	fprintf(f, "\n%f",result[i][j]);
                
            }
        }

    fclose(f);

}

/**
 * @brief Saves the solution to file in a format suitable for gnuplot: x y z, where z = phi(x,y)
 * 
 * @param bvp The Boundary Value Problem
 * @param filename The full filename. Can include path
 */
void bvp_create_gnuplot_data(bvp_t restrict bvp, const char *name){
	const unsigned int n = bvp->n;
	double **result = bvp->result;
	const double *x_vals = bvp->x_val;
	const double *y_vals = bvp->y_val;
	

	char filename[50];
	sprintf(filename,"solutions/gnuplot/%s.dat", name);
	
	FILE *f;
    f = fopen(filename, "w");

    if (f == NULL)
	{
    	printf("Error opening file!\n");
    	exit(1);
	}
	
    for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                
            	fprintf(f, "%f %f %f\n",x_vals[i], y_vals[j], result[i][j]);
                
            }
			fprintf(f, "\n");
        }
	fclose(f);

}