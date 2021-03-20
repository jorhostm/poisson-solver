#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#define NM_TRUE  1
#define NM_FALSE  0

void solve_poisson(double *T, const unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int neumann_x0, int neumann_x1, int neumann_y0, int neumann_y1);

void print_solution_to_file(double *T, const unsigned int n, const char *filename);

double get_value_at(double *T, const unsigned int n, const double x, const double y);

void shift_solution(double *T, const unsigned int n, const double delta_t);

#endif