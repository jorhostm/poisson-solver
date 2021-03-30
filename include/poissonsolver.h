#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#define NM_TRUE  1
#define NM_FALSE  0

// Poisson Boundary Value Problem type
typedef struct bvp_t bvp_t;

double get_value_at(bvp_t *bvp, const double x, const double y);

void shift_solution(bvp_t *bvp, const double delta_t);

void print_solution_to_file(bvp_t *bvp, const char *filename);

void create_gnuplot_data(bvp_t *bvp, const char *filename);

int solve_poisson_bvp(bvp_t *bvp, unsigned int useMultiGrid);

bvp_t *bvp_create(const unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int neumann_x0, int neumann_x1, int neumann_y0, int neumann_y1);

void bvp_destroy(bvp_t *bvp);

#endif