#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <poissonsolver.h>
#include <solver_utils.h>

struct bvp_t{
    double *result;
    unsigned int n;
	double *x_val;
	double *y_val;
    double (*phi)(double x, double y);
    double (*g)(double x, double y);
    double *b;
    int nm_x0;
    int nm_x1;
    int nm_y0;
    int nm_y1;
};



// SOR iterative solution
static int sor(bvp_t *bvp){
	unsigned int n = bvp->n;
	double *T = bvp->result;
	double *b = bvp->b;
	double h = 1.0/((double) n - 1);
	register double h2 = h*h;

	double *h2b = malloc(n*n*sizeof(double));
	for (int i = 1; i < n-1; i++)
	{
		for (int j = 1; j < n-1; j++)
		{
			h2b[i*n+j] = h2*b[i*n+j];
		}
		
	}
	

	register double omega = _omega(n);

	register double reltol = 0.000001;

    double rel_res = 1.0;

    int iterations = 0;

	register int nm_x0 = bvp->nm_x0;
	register int nm_x1 = bvp->nm_x1;
	register int nm_y0 = bvp->nm_y0;
	register int nm_y1 = bvp->nm_y1;

    
    while (rel_res > reltol){
        
		double Tmax = 0.0;
		double dTmax = 0.0;
        
		for (int i = 1; i < n-1; i++){

			int a1 = 1;	//	T(i-1,j)
			int a3 = 1;	// 	T(i+1,j)
			
			/* Check if there is a Neumann boundary at x = 0 */
			if(i == 1 && nm_x0){
				a1 = 0;
				a3 = 2;
			}
			/* Check if there is a Neumann boundary at x = 1 */
			else if(i == n-2 && nm_x1){
				a1 = 2;
				a3 = 0;
			}

            for (int j = 1; j < n-1; j++){

				int index = i*n + j;

				int a0 = 1;	// 	T(i,j-1)
            	int a2 = 1; // 	T(i,j+1)

				/* Check if there is a Neumann boundary at y = 0 */
                if(j == 1 && nm_y0){
                	a0 = 0;
                	a2 = 2;
                }
				/* Check if there is a Neumann boundary at y = 1 */
                else if(j == n-2 && nm_y1){
                	a0 = 2;
                	a2 = 0;
                }
				/* Compute the residual */
                double R = a0*T[i*n+j-1] + a1*T[(i-1)*n+j] - 4*T[index] + a2*T[i*n+j+1] + a3*T[(i+1)*n+j] - h2b[index];
                /* Apply relaxation coefficient */
				double dT = omega*R;
                /* Add the residual to the solution */
				T[index] += dT;
                
				/* Update the maximum value and residual */
				Tmax = fmax(fabs(T[index]), Tmax);
				dTmax = fmax(fabs(dT), dTmax);

            }
        }
		
        rel_res = dTmax/Tmax;
       	
        iterations++;
    } 
    
	free(h2b);
    return iterations;  
}

int solve_poisson_bvp(bvp_t *bvp){
	int it = sor(bvp);
	int n = bvp->n;
	
	if (bvp->nm_x0) _copy_neumann_border(bvp->result, n, 0,0,0,n-1, 1,1,0,n-1);
	if (bvp->nm_x1) _copy_neumann_border(bvp->result, n, n-1,n-1,0,n-1, n-2,n-2,0,n-1);
	if (bvp->nm_y0) _copy_neumann_border(bvp->result, n, 0,n-1,0,0, 0,n-1,1,1);
	if (bvp->nm_y1) _copy_neumann_border(bvp->result, n, 0,n-1,n-1,n-1, 0,n-1, n-2,n-2);

	return it;
}

void print_solution_to_file(bvp_t *bvp, const char *name){
	int n = bvp->n;
	double *result = bvp->result;

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
    for (int j = 0; j < n; j++){
            for (int i = 0; i < n; i++){
                
            	fprintf(f, "\n%f",result[j*n+i]);
                
            }
        }

    fclose(f);

}

void create_gnuplot_data(bvp_t *bvp, const char *name){
	int n = bvp->n;
	double *result = bvp->result;
	double *x = bvp->x_val;
	double *y = bvp->y_val;
	
	char filename[50];
	sprintf(filename,"solutions/gnuplot/%s.dat", name);
	
	FILE *f;
    f = fopen(filename, "w");

    if (f == NULL)
	{
    	printf("Error opening file!\n");
    	exit(1);
	}

    for (int j = 0; j < n; j++){
            for (int i = 0; i < n; i++){
                
            	fprintf(f, "%f %f %f\n",x[j], y[i], result[j*n+i]);
                
            }
        }

    fclose(f);

}

// Finds the value of T(x,y) using bilinear interpolation
double get_value_at(bvp_t *bvp, const double x, const double y){

	double *result = bvp->result;
	int n = bvp->n;

	// Domain check
	if (x < 0 || x > 1 || y < 0 || y > 1){
		printf("( %f , %f ) is outside of the domain\n", x, y);
    	exit(1);
	}

	double *x_values = create_linear_array(0,1,n);
	double *y_values = x_values;

	int i_left = 0;
	int i_right = n-1;
	int j_left = 0;
	int j_right = n-1;

	int delta_i = i_right - i_left;
	int delta_j = j_right - j_left;

	while (delta_i > 1 || delta_j > 1 ){

		delta_i = i_right - i_left;
		delta_j = j_right - j_left;
		
		int i_mid = (i_left + i_right) / 2;
		int j_mid = (j_left + j_right) / 2;

		if (x <= x_values[i_mid]){
			
			i_right = i_mid;

		}

		else{

			i_left = i_mid;

		}

		if (y <= y_values[j_mid]){
			
			j_right = j_mid;

		}

		else{

			j_left = j_mid;

		}
	
	}

	double t11 = result[i_left*n + j_left];
	double t21 = result[i_right*n + j_left];
	double t12 = result[i_left*n + j_right];
	double t22 = result[i_right*n + j_right];

	double h = x_values[1] - x_values[0];

	double x_normalised = (x-x_values[i_left])/h;
	double y_normalised = (y-y_values[j_left])/h;

	double a11 = t11;
	double a21 = t21 - t11;
	double a12 = t12 - t11;
	double a22 = t22 + t11 - (t21 + t12);

	double t_xy = a11 + a21 * x_normalised + a12 * y_normalised + a22 * x_normalised * y_normalised;

	free(x_values);
	//free(y_values);
	return t_xy;
}

void shift_solution(bvp_t *bvp, const double delta_t){

	double *result = bvp->result;
	int n = bvp->n;
	n = n*n;
	
	for (int i = 0; i < n; i++){

		result[i] += delta_t;

	}
}

bvp_t *bvp_create(const unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int neumann_x0, int neumann_x1, int neumann_y0, int neumann_y1){
	bvp_t *bvp = malloc(sizeof(bvp_t));

	bvp->result = calloc(n*n,sizeof(double));
	bvp->n = n;
	bvp->phi = phi;
	bvp->g = g;
	bvp->b = malloc(n*n*sizeof(double));
	bvp->nm_x0 = neumann_x0;
	bvp->nm_x1 = neumann_x1;
	bvp->nm_y0 = neumann_y0;
	bvp->nm_y1 = neumann_y1;

	double *x = create_linear_array(0, 1, n);
	double *y = x;

	bvp->x_val = x;
	bvp->y_val = y;

	// Initialize the RHS
	_init_values(bvp->b,g,n,x,y,0,n,0,n);
	
	// Initialize the boundaries
	if(!neumann_x0) _init_values(bvp->result,phi,n,x,y,0,n,0,1);
	if(!neumann_x1) _init_values(bvp->result,phi,n,x,y,0,n,n-1,n);
	if(!neumann_y0) _init_values(bvp->result,phi,n,x,y,0,1,0,n);
	if(!neumann_y1) _init_values(bvp->result,phi,n,x,y,n-1,n,0,n);

	return bvp;
}

void bvp_destroy(bvp_t *bvp){
	free(bvp->result);
	bvp->result = NULL;

	free(bvp->x_val);
	bvp->x_val = NULL;

	free(bvp->b);
	bvp->b = NULL;

	free(bvp);
}