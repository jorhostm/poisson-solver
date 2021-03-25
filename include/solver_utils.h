#ifndef SOLVER_UTILS_H
#define SOLVER_UTILS_H

double _omega(const unsigned int n);

void _init_values(double* T, double (*fun)(double x, double y), const unsigned int n, double* x, double* y, int i0, int i1, int j0, int j1);

void _copy_neumann_border(double* T, const unsigned int n, int dest_x0, int dest_x1, int dest_y0, int dest_y1, int src_x0, int src_x1, int src_y0, int src_y1);

double* create_linear_array(const double from, const double to, const unsigned int count);

#endif