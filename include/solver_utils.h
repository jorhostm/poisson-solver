#ifndef SOLVER_UTILS_H
#define SOLVER_UTILS_H

double _omega(const unsigned int n);

void _init_values(double **T, double (*fun)(double x, double y), const unsigned int n, double* x, double* y, int i0, int i1, int j0, int j1);

double* create_linear_array(const double from, const double to, const unsigned int count);

int is_admissable(const unsigned int n, const unsigned int multigrid_min);

unsigned int find_admissable(unsigned int n, const unsigned int multigrid_min);

#endif