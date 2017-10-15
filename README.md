# Airy-mKdV
Numerical implementation to solve Airy flow  and the modified Korteweg-de Vries equation. 

This also as a soliton visualizator.

##DESCRIPTION OF THE PROBLEM
This code computes the numerical solution of the two dimensional, periodic Initial Value Problem for the so-called Airy flow (a two dimensional dispersive equation).  The formulation is based in the theta-L formulation and solution theory, in which the dynamical variables (theta,sL) are used instead of (x,y) classical coordinates. The variable theta simplifies the formulation, so the arc-length sl must be constant over time. The code calculates the evolution for the theta-equation, converts back to (x,y) for visualization. 

##NUMERICAL METHOD

-The code is initialized with a 2 dimensional shape calling the function
 initiall, as part of
-It extracts the initial value for theta_0
 library.f
###TIME DISCRETIZATION
-The first step in the evolution of theta is computed with forward Euler
-The third and subsequent steps are computed using Adams Bashforth method
###SPACE DISCRETIZATION
it is computed with pseudo spectral methods (multiplication of nonlinear terms in physical space)
#SOFTWARE REQUIREMENTS
Before running the code you will need to install
- gfortran for doing the main computations
- Matlab 2012 or later for visualization and data analysis

#HOW THE RUN THE CODE
The code receives input data from

in.d 

containing the following parameters

 t        -time
 N      -space discretization size h=1/N, must be a poet of 2 for Fourier series to work
 dt      -time discretization size for computation
 ncut  -to measure error 
 tp      -time discretization size for printing (no need to print each computation)

you have to run

make 

in the terminal to update your files, don’t forget to initialize a smooth closed curve using the function initiall as part of  library.f


Then you need to type 

./hsclosed<in.d

in the terminal to update your initial parameters.



The main computation is done with 

airyflow.f,

which uses nl.f to handle the non-linear terms

#OUTPUT OF THE CODE
On file:
The code generates
- xp, yp files containing the evolution  of (x,y) positions 
- thetap files containing the evolution of the theta variable
- slp file containing the evolution of the arc-length (must be approx. constant)
- H1p,H2p,H3p containing information of numerical conservation of the underlying equation
-tp containing the vector of time
On screen:
It keeps track of
-x,y positions, 
- average radius
- Maximum and Minimum curvatures for tangent and normal vectors
- Length
- Area
- Center of mass
- Perturbation size

 VISUALIZE THE RESULTS 
You can use matlab to visualize the evolution of the curve that solves Airy flow and the evolution of the curvature that according to the modified Korteweg de Vries equation from the curvature of the curve. The function is called

visualize.m

and receives as a parameter the initial arc length (length of the initial curve)
