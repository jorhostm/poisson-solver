#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <poissonsolver.h>
#include <solver_utils.h>

struct bvp_t{
    double **result;
    unsigned int n;
	double *x_val;
	double *y_val;
    double (*phi)(double x, double y);
    double (*g)(double x, double y);
    double **b;
    int nm_x0;
    int nm_x1;
    int nm_y0;
    int nm_y1;
};

/* 
	Copy result of a coarser bvp over to a bvp with a finer grid 
	and interpolate the rest
*/
static int bvpcpy(const bvp_t *src_bvp, bvp_t *dest_bvp){
	const int src_n = src_bvp->n;
	const int dest_n = dest_bvp->n;
	double **src_r = src_bvp->result;
	double **dest_r = dest_bvp->result;

	if (src_n != dest_n/2)
	{
		return -1;
	}
	
	/* 
		The destination grid has 4 times the amount of points
		Copy over the first batch of datapoints,
		then interpolate the rest
	*/
	for (int i = 0; i < src_n; i++)
	{
		for (int j = 0; j < src_n; j++)
		{
			dest_r[2*i][2*j] = src_r[i][j];
		}
		
	}

	/* 
		In case of Neumann boundary at x = 0, interpolate
	*/
	if (dest_bvp->nm_x0)
	{
		for (int j = 0; j < src_n-1; j++)
		{
			dest_r[0][2*j+1] = 0.5*(dest_r[0][2*j] + dest_r[0][2*j+2]);
		}
	}

	/* 
		In case of Neumann boundary at y = 0, interpolate
	*/
	if (dest_bvp->nm_y0)
	{
		for (int i = 0; i < src_n-1; i++)
		{
			dest_r[2*i+1][0] = 0.5*(dest_r[2*i][0] + dest_r[2*i+2][0]);
		}
	}

	/* 
		In case of Neumann boundary at x = 1, interpolate
	*/
	if (dest_bvp->nm_x1)
	{
		const int i = src_n-1;
		for (int j = 0; j < src_n-1; j++)
		{
			dest_r[2*i][2*j+1] = 0.5*(dest_r[2*i][2*j] + dest_r[2*i][2*j+2]);
		}
		unsigned int n = dest_n;
		_copy_neumann_border(dest_r, n, n-1,n-1,0,n-1, n-2,n-2,0,n-1);
	}

	/* 
		In case of Neumann boundary at y = 1, interpolate
	*/
	if (dest_bvp->nm_y1)
	{
		const int j = src_n-1;
		for (int i = 0; i < src_n-1; i++)
		{
			dest_r[2*i+1][2*j] = 0.5*(dest_r[2*i][2*j] + dest_r[2*i+2][2*j]);
		}
		unsigned int n = dest_n;
		_copy_neumann_border(dest_r, n, 0,n-1,n-1,n-1, 0,n-1, n-2,n-2);
	}

	/*
		Second batch
		Interpolate between the points from the first batch \
	*/
	for (int i = 0; i < src_n-1; i++)
	{
		for (int j = 0; j < src_n-1; j++)
		{
			double corner_sum = dest_r[2*i][2*j];	// Top-Left
			corner_sum += dest_r[2*i+2][2*j];		// Top-Right
			corner_sum += dest_r[2*i][2*j+2];		// Bottom-Left
			corner_sum += dest_r[2*i+2][2*j+2];		// Bottom-Right
			dest_r[2*i+1][2*j+1] = 0.25*corner_sum;	// Middle = Average of the 4-corner-neighbourhood
		}
		
	}
	/*
		Third batch
		Interpolate between the points from the first and second batch
	*/
	for (int i = 0; i < src_n-1; i++)
	{
		for (int j = 0; j < src_n-1; j++)
		{
			double cross_sum = dest_r[2*i+2][2*j];	// Top
			cross_sum += dest_r[2*i+2][2*j+2];		// Botton
			cross_sum += dest_r[2*i+1][2*j+1];		// Left
			cross_sum += dest_r[2*i+3][2*j+1];		// Right
			dest_r[2*i+2][2*j+1] = 0.25*cross_sum;	// Middle = Average of the 4-cross-neighbourhood
		}
		
	}
	
	/*
		Fourth batch
		Interpolate between the points from the first, second and third batch
	*/
	for (int i = 0; i < src_n-1; i++)
	{
		for (int j = 0; j < src_n-1; j++)
		{
			double sum = dest_r[2*i][2*j+1];	// Top-Left
			sum += dest_r[2*i+1][2*j+1];		// Top
			sum += dest_r[2*i+2][2*j+1];		// Top-Right
			sum += dest_r[2*i][2*j+2];			// Left
			sum += dest_r[2*i+2][2*j+2];		// Right
			sum += dest_r[2*i][2*j+3];			// Bottom-Left
			sum += dest_r[2*i+1][2*j+3];		// Bottom
			sum += dest_r[2*i+2][2*j+3];		// Bottom-Right
			dest_r[2*i+1][2*j+2] = 0.125*sum;	// Middle = Average of 8-neighbourhood
		}
	}

	return 0;
}

/*
 Solve the BVP using the Successive Over-Relaxation iterative method
*/
static int sor(bvp_t *bvp){
	
	const unsigned int n = bvp->n;
	double **T = bvp->result;
	double **b = bvp->b;
	double h = 1.0/((double) n - 1);
	double h2 = h*h;
	double omega = _omega(n);

	/* Relative tolerence and relative residual*/
	const double reltol = 1e-5;
    double rel_res = 1.0;

    int iterations = 0;

	const unsigned int nm_x0 = bvp->nm_x0;
	const unsigned int nm_x1 = bvp->nm_x1;
	const unsigned int nm_y0 = bvp->nm_y0;
	const unsigned int nm_y1 = bvp->nm_y1;

    
    while (rel_res > reltol){
        
		double T_sum = reltol; // T_sum != 0
		double dT_sum = 0.0;
        
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
                double R = a0*T[i][j-1] + a1*T[i-1][j] + a2*T[i][j+1] + a3*T[i+1][j] - 4*T[i][j] - h2*b[i][j];
                /* Apply relaxation factor */
				double dT = omega*0.25*R;
                /* Add the residual to the solution */
				T[i][j] += dT;
                
				
				T_sum += fabs(T[i][j]);
				dT_sum += fabs(dT);

            }
        }
		
        rel_res = dT_sum/T_sum;
       	
        iterations++;
    } 
    return iterations;  
}

int solve_poisson_bvp(bvp_t *bvp, unsigned int useMultiGrid){
	int n = bvp->n;
	int iterations = 0;
	
	if(n >= 256 && n % 2 == 0 && useMultiGrid){
		bvp_t *coarse_bvp = bvp_create(n/2, bvp->phi,bvp->g,bvp->nm_x0, bvp->nm_x1,bvp->nm_y0,bvp->nm_y1);
		iterations += solve_poisson_bvp(coarse_bvp, --useMultiGrid);
		bvpcpy(coarse_bvp, bvp);
		bvp_destroy(coarse_bvp);
	}
	iterations += sor(bvp);
	
	if (bvp->nm_x0) _copy_neumann_border(bvp->result, n, 0,0,0,n-1, 1,1,0,n-1);
	if (bvp->nm_x1) _copy_neumann_border(bvp->result, n, n-1,n-1,0,n-1, n-2,n-2,0,n-1);
	if (bvp->nm_y0) _copy_neumann_border(bvp->result, n, 0,n-1,0,0, 0,n-1,1,1);
	if (bvp->nm_y1) _copy_neumann_border(bvp->result, n, 0,n-1,n-1,n-1, 0,n-1, n-2,n-2);

	return iterations;
}


void print_solution_to_file(bvp_t *bvp, const char *name){
	int n = bvp->n;
	double **result = bvp->result;

	mkdir("solutions");
	mkdir("solutions/opengl");

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

void create_gnuplot_data(bvp_t *bvp, const char *name){
	const unsigned int n = bvp->n;
	double **result = bvp->result;
	const double *x_vals = bvp->x_val;
	const double *y_vals = bvp->y_val;
	
	mkdir("solutions");
	mkdir("solutions/gnuplot");
	mkdir("solutions/gnuplot/images");
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
        }

    fclose(f);
}

// Finds the value of T(x,y) using bilinear interpolation
double get_value_at(bvp_t *bvp, const double x, const double y){

	double **result = bvp->result;
	int n = bvp->n;

	// Domain check
	if (x < 0 || x > 1 || y < 0 || y > 1){
		printf("( %f , %f ) is outside of the domain\n", x, y);
    	exit(1);
	}

	double *x_values = bvp->x_val;
	double *y_values = bvp->y_val;

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

	double t11 = result[i_left][j_left];
	double t21 = result[i_right][j_left];
	double t12 = result[i_left][j_right];
	double t22 = result[i_right][j_right];

	double h = x_values[1] - x_values[0];

	double x_normalised = (x-x_values[i_left])/h;
	double y_normalised = (y-y_values[j_left])/h;

	double a11 = t11;
	double a21 = t21 - t11;
	double a12 = t12 - t11;
	double a22 = t22 + t11 - (t21 + t12);

	double t_xy = a11 + a21 * x_normalised + a12 * y_normalised + a22 * x_normalised * y_normalised;

	return t_xy;
}

void shift_solution(bvp_t *bvp, const double delta_t){

	double *result = bvp->result[0];
	int n = bvp->n;
	int n2 = n*n;
	
	for (int i = 0; i < n2; i++){

		result[i] += delta_t;

	}
}

bvp_t *bvp_create(const unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int neumann_x0, int neumann_x1, int neumann_y0, int neumann_y1){
	bvp_t *bvp = malloc(sizeof(bvp_t));
	
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
	double *datapoints = calloc(n*n,sizeof(double));
	bvp->b = calloc(n, sizeof(double*));
	for (int i = 0; i < n; i++)
	{
		bvp->b[i] = datapoints + i*n;
	}
	_init_values(bvp->b,g,n,x,y,0,n,0,n);

	// Initialize the 2D solution array
	datapoints = calloc(n*n,sizeof(double));
	bvp->result = calloc(n, sizeof(double*));
	for (int i = 0; i < n; i++)
	{
		bvp->result[i] = datapoints + i*n;
	}
	
	
	// Initialize the boundaries
	if(!neumann_x0) _init_values(bvp->result,phi,n,x,y,0,n,0,1);
	if(!neumann_x1) _init_values(bvp->result,phi,n,x,y,0,n,n-1,n);
	if(!neumann_y0) _init_values(bvp->result,phi,n,x,y,0,1,0,n);
	if(!neumann_y1) _init_values(bvp->result,phi,n,x,y,n-1,n,0,n);

	return bvp;
}

void bvp_destroy(bvp_t *bvp){
	free(bvp->result[0]);
	free(bvp->result);
	bvp->result = NULL;

	free(bvp->x_val);
	bvp->x_val = NULL;
	bvp->y_val = NULL;

	free(bvp->b[0]);
	free(bvp->b);
	bvp->b = NULL;

	free(bvp);
}