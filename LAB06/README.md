To compile: gcc -O main.c -o ./test -lm

Please note that computing might take a while.
The code does not create the ./u folder automatically, so you have to make sure it exists before executing the program.

******************************************************************************

This C program models the 2D transport equation to simulate how a scalar quantity moves through a velocity field over time.
 It calculates the initial velocity vectors based on an imported stream function, applies an implicit relaxation scheme 
to solve the time-dependent advection-diffusion, and outputs both spatial data for animation and diagnostic data 
(total mass conservation and average displacement) to verify numerical stability.
