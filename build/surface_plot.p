set title ARG1
set pm3d interpolate 1,1
set dgrid3d 50,50 qnorm 4
data = sprintf("solutions/gnuplot/%s.dat", ARG1)
splot data with pm3d
pause mouse close