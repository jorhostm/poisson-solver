#include <stdlib.h>
#include <math.h>
#include <stdio.h>


static double _omega(const unsigned int n){
	
	double pi = acos(-1);
	double rho = cos(pi/(n));
	
	return 2.0/(1+sqrt(1-pow(rho,2)));
}


static void _init_values(double* T, double (*f)(double x, double y), const unsigned int n, double* x, double* y, int i0, int i1, int j0, int j1){
	for(int i = i0; i < i1; i++){
		for(int j = j0; j < j1; j++){
			T[i*n+j] = (*f)(x[i],y[j]);
		}
	}
}

static void _copy_neumann_border(double* T, const unsigned int n, int dest_x0, int dest_x1, int dest_y0, int dest_y1, int src_x0, int src_x1, int src_y0, int src_y1){

	int dest_dx = dest_x1 - dest_x0 + 1;
	int dest_dy = dest_y1 - dest_y0 + 1;
	int src_dx = src_x1 - src_x0 + 1;
	int src_dy = src_y1 - src_y0 + 1;

	if (dest_dx != src_dx || dest_dy != src_dy){
		printf("Error in copy border: Dimensions do not match\n" );
		exit(1);
	}

	for(int i = 0; i < dest_dx; i++){
		for(int j = 0; j < dest_dy; j++){
			T[(dest_x0 + i)*n+(dest_y0 + j)] = T[(src_x0 + i)*n+(src_y0 + j)];
		}
	}
}

static double* create_linear_array(const double from, const double to, const unsigned int count){

	double h = (to - from)/((double) count - 1);

	double *arr = malloc( count * sizeof(double) );

	arr[0] = from;
	arr[count-1] = to;


	for(int i = 1; i < (count-1); i++){
		arr[i] = arr[i-1] + h;
	}

	return arr;

}

static int _sor_iterate(double* T, double* b, const unsigned int n, const unsigned int start, const unsigned int end, const int nm_x0, const int nm_x1, const int nm_y0, const int nm_y1){
	// SOR iterative solution
	int n2 = n*n;
	int count = end-start;
	double h = 1.0/((double) n - 1);
	double h2 = h*h;
	double omega = _omega(n);

	double reltol = 0.0001;

    double rel_res = 1.0;

    int iterations = 0;
    
    
    while (rel_res > reltol){
        
		double Tmax = 0.0;
		double dTmax = 0.0;
        
		for (int i = start; i < end; i++){
            for (int j = start; j < end; j++){

            	int a0 = 1;
            	int a1 = 1;
            	int a2 = 1;
            	int a3 = 1;

				if(i == 1 && nm_x0){
                	a1 = 0;
                	a3 = 2;
                }
                else if(i == n-2 && nm_x1){
                	a1 = 2;
                	a3 = 0;
                }

                if(j == 1 && nm_y0){
                	a0 = 0;
                	a2 = 2;
                }
                else if(j == n-2 && nm_y1){
                	a0 = 2;
                	a2 = 0;
                }

                double R = a0*T[i*n+j-1] + a1*T[(i-1)*n+j] - 4*T[i*n+j] + a2*T[i*n+j+1] + a3*T[(i+1)*n+j] - h2*b[i*n+j];
                double dT = 0.25*omega*R;
                T[i*n+j] += dT;
                
				Tmax = fmax(fabs(T[i*n+j]), Tmax);
				dTmax = fmax(fabs(dT), dTmax);

            }
        }

        rel_res = dTmax/Tmax;
       	
        iterations++;
    } 
    

    return iterations;  
}

void solve_poisson(double *T, const unsigned int n, double (*phi)(double x, double y), double (*g)(double x, double y), int neumann_x0, int neumann_x1, int neumann_y0, int neumann_y1){

	double x0, y0 = 0;
	double xEnd, yEnd = 1;
	double h = 1.0/((double) n - 1);

	/* Amount of gridpoints*/
	int n2 = n*n;

	double *x = create_linear_array(0, 1, n);
	double *y = x;

	
	double* b = malloc(n2*sizeof(double));

	_init_values(b,(*g),n,x,y,0,n,0,n);
	
	_init_values(T,(*phi),n,x,y,0,n,0,1);
	_init_values(T,(*phi),n,x,y,0,n,n-1,n);
	_init_values(T,(*phi),n,x,y,0,1,0,n);
	_init_values(T,(*phi),n,x,y,n-1,n,0,n);
        
    int it = _sor_iterate(T, b, n, 1, n-1, neumann_x0, neumann_x1, neumann_y0, neumann_y1);
   
    printf("Finnished after %d iterations\n", it);
	
	if (neumann_x0) _copy_neumann_border(T, n, 0,0,0,n-1, 1,1,0,n-1);
	if (neumann_x1) _copy_neumann_border(T, n, n-1,n-1,0,n-1, n-2,n-2,0,n-1);
	if (neumann_y0) _copy_neumann_border(T, n, 0,n-1,0,0, 0,n-1,1,1);
	if (neumann_y1) _copy_neumann_border(T, n, 0,n-1,n-1,n-1, 0,n-1, n-2,n-2);
	

	free(x);
//	free(y);
	free(b);


}

void print_solution_to_file(double *T, const unsigned int n, const char *filename){

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
                
            	fprintf(f, "\n%lf",T[j*n+i]);
                
            }
        }

    fclose(f);

}

// Finds the value of T(x,y) using bilinear interpolation
double get_value_at(double *T, const unsigned int n, const double x, const double y){

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

	double t11 = T[i_left*n + j_left];
	double t21 = T[i_right*n + j_left];
	double t12 = T[i_left*n + j_right];
	double t22 = T[i_right*n + j_right];

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

	/*Debug stuff*/
	/*
	printf("i_left = %d, i_right = %d, j_left = %d, j_right = %d\n", i_left, i_right, j_left, j_right);
	printf("t11 = %f, t21 = %f, t12 = %f, t22 = %f\n", t11, t21, t12, t22);
	*/
	return t_xy;
}

void shift_solution(double *T, const unsigned int n, const double delta_t){

	const int nT = n*n;
	
	for (int i = 0; i < nT; i++){

		T[i] += delta_t;

	}


}
