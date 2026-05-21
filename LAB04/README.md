in order to compile: gcc main.c -o ./test -lm

********************************************


This code implements a multigrid relaxation solver for the Laplace equation 
on a 2D grid to compute the electrostatic potential. 
It iteratively refines the solution across different grid resolutions, 
starting from a coarse mesh and interpolating to a fine mesh to achieve 
convergence efficiently.
