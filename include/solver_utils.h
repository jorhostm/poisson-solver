/**
 * @file solver_utils.h
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This header file decleares the helper functions for the poissonsolver. 
 * See solver_utils.c for implementation details
 * @version 1.0
 * @date 2021-05-04
 * 
 */
#ifndef SOLVER_UTILS_H
#define SOLVER_UTILS_H

#include <pthread.h>

struct bvp_t{

    double **result;
    unsigned int n;
    double *x_val;
    double *y_val;
    double (*phi)(double x, double y);
    double (*g)(double x, double y);
    double **b;
    int nm_flags;
};

double bvp_utils_get_next_omega(double omega, const unsigned int n);

void bvp_utils_add_values(double **T, double (*func)(double x, double y), const unsigned int n, double* x, double* y, const unsigned int i_start, const unsigned int i_end, const unsigned int j_start, const unsigned int j_end);

double* bvp_utils_create_linear_array(const double from, const double to, const unsigned int count);

#endif