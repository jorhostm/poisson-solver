#ifndef SOLVER_UTILS_H
#define SOLVER_UTILS_H

/**
 * @brief Computes the theoretical optimal relaxation factor omega for the SOR iteration method.
 *	Works for Laplace and Poisson equations in a rectangular domain.
 *	In this case the domain is a n-by-n grid with x,y = [0, 1].
 *	The specral radius rho(A) = cos(PI*h), where h is the step size = 1/(n-1)
 * 
 * @param n The amount of discretized points along an axis.
 * @return double The theoretical optimal relaxation factor
 */
double get_omega(const unsigned int n);

void add_values(double **T, double (*fun)(double x, double y), const unsigned int n, double* x, double* y, const unsigned int i_start, const unsigned int i_end, const unsigned int j_start, const unsigned int j_end);

void sor_iterate(double **T, double **b, const unsigned int n, const double omega, double *Tsum, double *dTsum, const unsigned int i_start, const unsigned int i_end, const unsigned int j_start, const unsigned int j_end);

double* create_linear_array(const double from, const double to, const unsigned int count);

int is_admissable(const unsigned int n, const unsigned int multigrid_min);

unsigned int find_admissable(unsigned int n, const unsigned int multigrid_min);

#endif