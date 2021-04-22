#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#define NM_X0           0x1
#define NM_X1           0x2
#define NM_Y0           0x4
#define NM_Y1           0x8
#define NM_ONLY          NM_X0 | NM_X1 | NM_Y0 | NM_Y1

#define MULTIGRID_MIN   200

// Poisson Boundary Value Problem type
typedef struct bvp_t* bvp_t;

double get_value_at(bvp_t bvp, const double x, const double y);

void shift_solution(bvp_t bvp, const double delta_t);

void print_solution_to_file(bvp_t bvp, const char *filename);

void create_gnuplot_data(bvp_t bvp, const char *filename);

int solve_poisson_bvp(bvp_t bvp, unsigned int use_multigrid, const double reltol);

bvp_t bvp_create(const unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int nm_flags);

void bvp_destroy(bvp_t bvp);

#endif