# poisson-solver
A general solver for the poisson equation, supporting Dirichlet and Neumann boundaries. Example solutions:

![](https://raw.githubusercontent.com/jorhostm/poisson-solver/master/build/solutions/gnuplot/images/Mixed.png)
![](https://raw.githubusercontent.com/jorhostm/poisson-solver/master/build/solutions/gnuplot/images/Dirichlet1.png)
![](https://raw.githubusercontent.com/jorhostm/poisson-solver/master/build/solutions/gnuplot/images/Neumann1.png)
## Compilation requirements
- GCC
- Gnu make
- FreeGLUT
- GLU
- GNU-getopt
- pthreads

## Compilation
To make all, run:
```
make all
```

## Usage
Usage of the solver is demonstrated in the provided example program.
To run the example program for a 100-by-100 grid, relative tolerance of 10^-5 and saving the solution:
```
./example.out -n 100 -r 1e-5 -s
```
To use multigrids, use the -m flag:
```
./example.out -n 1000 -r 1e-6 -s -m
```
To only solve a limited number of the example BVPs, use the -p flag:
```
./example.out -n 10000 -r 1e-6 -m -p 1
```
To use the visualisation demo program, provide it with the path to the solution file:
```
./demo.out solutions/opengl/Mixed.dat
```
Rotate the solution using WASD. To zoom press W and S while holding down left shift. Press space to swap between 3D and 2D. Press ESC to quit.
The testing program needs several inputs: From and to value for n, step size of n, whether to use multigrids, which problem to solve, the relative tolerance and the name of the file to store the results. Example:
```
./test.out -f 10 -t 3000 -s 10 -m -p 0 -r 1e-5 -n multigridMixed.dat
```
