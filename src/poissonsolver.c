#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <poissonsolver.h>
#include <solver_utils.h>

struct bvp_t{
	double *datapoints;
    double **result;
    unsigned int n;
	double *x_val;
	double *y_val;
    double (*phi)(double x, double y);
    double (*g)(double x, double y);
    double **b;
    int nm_flags;
};

/* 
	Copy result of a coarser bvp over to a bvp with a finer grid 
	and interpolate the remaining datapoints
*/
static int bvpcpy(const bvp_t src_bvp, bvp_t dest_bvp){
	
	const int src_n = src_bvp->n;
	const int dest_n = dest_bvp->n;
	double **src_r = src_bvp->result;
	double **dest_r = dest_bvp->result;


	if (src_n != dest_n/2 + 1)
	{
		return -1;
	}
	
	/* 
		The destination grid has about 4 times the amount of points as the source grid,
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
		In case of Neumann boundary at x = 0, 
		interpolate the remaining datapoints on that boundary
	*/
	if (dest_bvp->nm_flags & NM_X0)
	{
		for (int j = 0; j < src_n-1; j++)
		{
			dest_r[0][2*j+1] = 0.5*(dest_r[0][2*j] + dest_r[0][2*j+2]);
		}
	}

	/* 
		In case of Neumann boundary at y = 0, 
		interpolate the remaining datapoints on that boundary
	*/
	if (dest_bvp->nm_flags & NM_Y0)
	{
		for (int i = 0; i < src_n-1; i++)
		{
			dest_r[2*i+1][0] = 0.5*(dest_r[2*i][0] + dest_r[2*i+2][0]);
		}
	}

	/* 
		In case of Neumann boundary at x = 1, 
		interpolate the remaining datapoints on that boundary
	*/
	if (dest_bvp->nm_flags & NM_X1)
	{
		const int i = src_n-1;
		for (int j = 0; j < src_n-1; j++)
		{
			dest_r[2*i][2*j+1] = 0.5*(dest_r[2*i][2*j] + dest_r[2*i][2*j+2]);
		}
	}

	/* 
		In case of Neumann boundary at y = 1, 
		interpolate the remaining datapoints on that boundary
	*/
	if (dest_bvp->nm_flags & NM_Y1)
	{
		const int j = src_n-1;
		for (int i = 0; i < src_n-1; i++)
		{
			dest_r[2*i+1][2*j] = 0.5*(dest_r[2*i][2*j] + dest_r[2*i+2][2*j]);
		}
	}

	/*
		Second batch
		Interpolate between the points from the first batch

				◯   ◯
		x 	:=
				◯	◯
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

				  ◯ 
		△ 	:=	x	x
				  ◯
	*/
	// Account for odd numbered grid size, so it does not go out of bounds
	for (int i = 0; i < src_n-2; i++)
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

				△ x △
		□ 	:=	◯	◯
				△ x △
	*/
	for (int i = 0; i < src_n-1; i++)
	{
		// Account for odd numbered grid size, so it does not go out of bounds
		for (int j = 0; j < src_n-2; j++)
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
static int sor(bvp_t bvp, const double reltol){
	
	const unsigned int n = bvp->n;
	double **T = bvp->result;
	double **b = bvp->b;
	double h = 1.0/((double) n - 1);
	double h2 = h*h;
	double omega = _omega(n);

	/* Relative residual*/
	
    double rel_res = 1.0;

    int iterations = 0;

	const unsigned int i_start 	= bvp->nm_flags & NM_X0 ? 0 : 1;
	const unsigned int i_end 	= bvp->nm_flags & NM_X1 ? n : n-1;
	const unsigned int j_start 	= bvp->nm_flags & NM_Y0 ? 0 : 1;
	const unsigned int j_end 	= bvp->nm_flags & NM_Y1 ? n : n-1;

    
    while (rel_res > reltol){
        
		double T_sum = 0.0;
		double dT_sum = 0.0;
        
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

            for (int j = j_start; j < j_end; j++){

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
				
				//T_sum = fmax(fabs(T[i][j]), T_sum);
				//dT_sum = fmax(fabs(dT), dT_sum);
				T_sum += fabs(T[i][j]);
				dT_sum += fabs(dT);

            }
        }
		
        rel_res = dT_sum/T_sum;
       	
        iterations++;
    } 
    return iterations;  
}

int solve_poisson_bvp(bvp_t bvp, unsigned int use_multigrid, const double reltol){
	int n = bvp->n;
	int iterations = 0;
	
	if(n >= MULTIGRID_MIN && use_multigrid){
		bvp_t coarse_bvp = bvp_create(n/2 + 1, bvp->phi,bvp->g,bvp->nm_flags);
		iterations += solve_poisson_bvp(coarse_bvp, use_multigrid, reltol);
		bvpcpy(coarse_bvp, bvp);
		bvp_destroy(coarse_bvp);
	}
	
	iterations += sor(bvp, reltol);

	return iterations;
}


void print_solution_to_file(bvp_t bvp, const char *name){
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

void create_gnuplot_data(bvp_t bvp, const char *name){
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

// Finds the value of T(x,y) using bilinear interpolation
double get_value_at(bvp_t bvp, const double x, const double y){

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

void shift_solution(bvp_t bvp, const double delta_t){

	double *result = bvp->result[0];
	int n = bvp->n;
	int n2 = n*n;
	
	for (int i = 0; i < n2; i++){

		result[i] += delta_t;

	}
}

bvp_t bvp_create(unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int nm_flags){
	bvp_t bvp = malloc(sizeof(struct bvp_t));
	n = find_admissable(n, MULTIGRID_MIN);
	
	bvp->n = n;
	bvp->phi = phi;
	bvp->g = g;
	bvp->nm_flags = nm_flags;

	double *x = create_linear_array(0, 1, n);
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
	_init_values(bvp->b,g,n,x,y,0,n,0,n);

	// Initialize the 2D solution array, including "ghost cells" around the boundary, to prevent out-of-bound access
	double *datapoints = calloc((n+2)*(n+2),sizeof(double));
	bvp->datapoints = datapoints;
	bvp->result = calloc(n+2, sizeof(double*));
	for (int i = 0; i < n+2; i++)
	{
		bvp->result[i] = datapoints + i*(n)+1;
	}
	(bvp->result)++;
	
	// Initialize the boundaries
	if(!(bvp->nm_flags & NM_X0)) _init_values(bvp->result,phi,n,x,y,0,n,0,1);
	if(!(bvp->nm_flags & NM_X1)) _init_values(bvp->result,phi,n,x,y,0,n,n-1,n);
	if(!(bvp->nm_flags & NM_Y0)) _init_values(bvp->result,phi,n,x,y,0,1,0,n);
	if(!(bvp->nm_flags & NM_Y1)) _init_values(bvp->result,phi,n,x,y,n-1,n,0,n);
	
	return bvp;
}

void bvp_destroy(bvp_t bvp){
	free(bvp->datapoints);
	free(--(bvp->result));
	bvp->result = NULL;

	free(bvp->x_val);
	bvp->x_val = NULL;
	bvp->y_val = NULL;

	free(bvp->b);
	bvp->b = NULL;

	free(bvp);
}