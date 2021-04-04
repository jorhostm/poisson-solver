set term png
set output sprintf("solutions/gnuplot/images/%s.png", ARG1)
set title ARG1
set ticslevel 0
set pm3d interpolate 1,1
data = sprintf("solutions/gnuplot/%s.dat", ARG1)
splot data with pm3d