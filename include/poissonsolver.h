/**
 * @file poissonsolver.h
 * @author Jørgen Høstmark (jorgenhostm@gmail.com)
 * @brief This header file defines the API functions of the Poisson equation solver
 * @version 1.0
 * @date 2021-04-23
 * 
 * 
 */

#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H


//Bit flags used to set Neumann boundaries
#define NM_NONE         0x0
#define NM_X0           0x1
#define NM_X1           0x2
#define NM_Y0           0x4
#define NM_Y1           0x8
#define NM_ONLY         (NM_X0 | NM_X1 | NM_Y0 | NM_Y1)

#ifndef MULTIGRID_MIN
#define MULTIGRID_MIN   100
#endif

/**
 * @brief Poisson Boundary Value Type
 */
typedef struct bvp_t* bvp_t;


/**
 * @brief Create and initialize a Poisson Boundary Value Problem
 * 
 * @param n The number of discretised points along an axis
 * @param phi The function phi(x,y) used to prescribe Dirichlet boundaries. If NULL, is threated as phi(x,y) = 0
 * @param g The function g(x,y) on the right side of the Poisson equation. If NULL, is threated as g(x,y) = 0
 * @param nm_flags Flags that is used to specify Neumann boundaries, if any. If there isn't any, set it as 0
 * @return bvp_t The Boundary Value Problem with an initial n-by-n grid with prescribed Dirichlet boundaries, if any
 */
bvp_t bvp_create(const unsigned int n, double (*phi)(double x, double y),double (*g)(double x, double y), int nm_flags);


/**
 * @brief Solve the given Poisson Boundary Value Problem
 * 
 * @param bvp The Boundary Value Problem
 * @param use_multigrid Whether or not to us multigrid preconditioning
 * @param reltol The relative tolerence for the relative residual convergence. Smaller => greater accuracy but more iterations
 * @return int The number of iterations used to reach convergence
 */
int bvp_solve(bvp_t restrict bvp, const unsigned int use_multigrid, const double reltol);


/**
 * @brief Iterate using Succesive Over-Relaxation for the Red-Black-Gauss-Seidel method
 * 
 * @param bvp The Boundary Value Problem
 * @param reltol The relative tolerence for the relative residual used to determine convergence. 
 *              Smaller => greater accuracy ,but more iterations
 * @return int The number of iterations used to reach convergence
 */
int red_black_sor(bvp_t restrict bvp, const double reltol);


/**
 * @brief Get the value at phi(x,y) using bilinear interpolation
 * 
 * @param bvp The Boundary Value Problem
 * @param x The x-coordinate
 * @param y The y-coordinate
 * @return double 
 */
double bvp_get_value_at(const bvp_t restrict bvp, const double x, const double y);


/**
 * @brief Shifts the entire solution by dz. Used to correct a solution with only Neumann boundaries
 * 
 * @param bvp The Boundary Value Problem
 * @param dz The value to shift the solution by
 */
bvp_t bvp_shift_solution(bvp_t restrict bvp, const double dz);


/**
 * @brief Copy result of a bvp over to another bvp. If destination BVP is NULL, creates a new one. Uses bilinear interpolation
 * 
 * @param src_bvp The source BVP solution
 * @param dest_bvp The destination BVP solution
 * @return bvp_t The destination BVP if successful, NULL otherwise
 */
bvp_t bvp_copy(const bvp_t restrict src, bvp_t restrict dst);


/**
 * @brief Free up memory occupied by the Boundary Value Problem
 * 
 * @param bvp The Boundary Value Problem
 */
void bvp_destroy(bvp_t restrict bvp);


/**
 * @brief Saves the solution to file
 * @param bvp The Boundary Value Problem
 * @param filename
 */
void bvp_print_solution_to_file(const bvp_t restrict bvp, const char *filename);


/**
 * @brief Saves the solution to file in a format suitable for gnuplot: x y z, where z = phi(x,y)x y z, where z = phi(x,y)
 * 
 * @param bvp The Boundary Value Problem
 * @param filename
 */
void bvp_create_gnuplot_data(const bvp_t restrict bvp, const char *filename);

#endif