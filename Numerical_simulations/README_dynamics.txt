DYNAMICS PROGRAM
================

OVERVIEW
--------

The file dynamics.c numerically integrates the three-species
resource-consumer-predator model for one fixed set of parameters and initial
conditions. It generates a consumer-predator phase-plane trajectory that can be
used to inspect the long-term dynamics of the system. The output command can be
easily adapted to generate time series for R, C and P or phase portraits for
other pairs of state variables.


FILE
----

  dynamics.c   C source code containing the model and numerical integrator


MODEL VARIABLES
---------------

  R = x[0]   resource density
  C = x[1]   consumer density
  P = x[2]   predator density

The program implements

  dR/dt = R(1 - D - R) - x_c y_c R C/(R + R_0),

  dC/dt = x_c C[y_c R/(R + R_0) - 1]
          - Phi y_p C/(C + C_0),

  dP/dt = Phi y_p C/(C + C_0) - x_p P,

where

  Phi = x_p(1 - sigma)P + x_i sigma P^2.


PARAMETERS IN THE SUPPLIED SOURCE
---------------------------------

  x_c   = 0.4
  y_c   = 2.099
  R_0   = 0.16129
  x_p   = 0.08
  x_i   = 0.08
  y_p   = 2.876
  C_0   = 0.5
  D     = 0.001
  sigma = 0.0

The control parameters D and sigma are defined inside the function var0(). To
study another ecological regime, edit these values and recompile the program.


INITIAL CONDITIONS
------------------

  R(0) = 0.55
  C(0) = 0.35
  P(0) = 0.8


NUMERICAL METHOD
----------------

The equations are integrated with an embedded seventh/eighth-order
Runge-Kutta-Fehlberg method using adaptive time stepping.

  final integration time:       time_f = 10,000
  local error tolerance:        1e-15
  minimum permitted time step:  HMIN = 1e-5
  maximum permitted time step:  HMAX = 0.1

Population densities below the quasi-extinction threshold of 1e-12 are set to
zero, making extinction irreversible within the simulation.


OUTPUT
------

The program creates:

  CP.txt

Each row contains two space-separated columns:

  column 1: consumer density, C
  column 2: predator density, P

The output therefore represents the trajectory projected onto the (C,P)
phase plane. Time and resource density are not written to this file. The final
values of R, C and P are printed to standard output when the integration ends.

The fprintf() instruction can be readily modified to obtain other outputs. For
example, writing time together with R, C and P produces the three population
time series, whereas writing any two state variables produces the corresponding
phase portrait.


COMPILATION
-----------

A C compiler and the standard mathematics library are required. With GCC:

  gcc -O3 -Wall -o dynamics dynamics.c -lm


RUNNING
-------

Run the compiled program from a writable directory:

  ./dynamics

The file CP.txt will be created in the current working directory. An existing
file with that name will be overwritten.


PLOTTING THE PHASE PLANE
------------------------

For example, CP.txt can be plotted with gnuplot using:

  plot "CP.txt" using 1:2 with lines title "trajectory"

Label the horizontal axis C and the vertical axis P.
