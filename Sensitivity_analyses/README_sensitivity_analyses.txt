SENSITIVITY-ANALYSIS PROGRAMS
=============================

This directory contains two C programs used to evaluate the robustness of the
three-species resource-consumer-predator model to variation in model parameters
and initial conditions. Each sampled parameter set is integrated numerically
and classified according to the species persisting at the end of the run.


FILES
-----

1. sensitivity-break.c

   Performs a Monte Carlo sensitivity analysis over selected combinations of
   habitat loss (D) and the fraction of cooperative hunters (sigma). Random
   parameter values and initial conditions are generated independently for
   each simulation. The random-number generator is initialised from the system
   time, so independent executions need not produce identical samples.

   In the supplied version:

     - D is evaluated from 0.9 to values below 1.0 in increments of 0.1;
     - sigma is evaluated from 0.0 to 1.0 in increments of 0.1;
     - 1,000,000 parameter/initial-condition combinations are simulated for
       each (D, sigma) pair;
     - each trajectory is integrated up to t = 30,000;
     - the local numerical-error tolerance is 1e-15.

   The program writes the estimated outcome probabilities to:

     D_sigma_RCP1e7_D1.txt   full coexistence (R, C and P)
     D_sigma_RC1e7_D1.txt    resource-consumer coexistence (R and C)
     D_sigma_R1e7_D1.txt     resource persistence only (R)
     D_sigma_01e7_D1.txt     extinction of all three populations

   Each output row contains D, sigma and the corresponding estimated
   probability. Note that the filenames are retained from the original
   analysis workflow; the supplied source currently sets ITERS = 1e6.


2. sensitivity-break_fixed_seed.c

   Performs the same type of Monte Carlo classification for one fixed pair of
   control parameters. In the supplied version, D = 0.9 and sigma = 1.0.
   The calls that would initialise the random-number generator from the system
   time are commented out. Consequently, drand48() starts from its deterministic
   default state on POSIX-compatible systems, allowing the same random sample
   to be regenerated on repeated runs on the same platform.

   In the supplied version:

     - 5,000,000 parameter/initial-condition combinations are simulated;
     - each trajectory is integrated up to t = 20,000;
     - the local numerical-error tolerance is 1e-15;
     - the four estimated outcome probabilities are printed to standard output.

   To study another combination, edit the constants D and sigma near the
   beginning of the source file and recompile.


MODEL VARIABLES AND PARAMETERS
------------------------------

State variables:

  R   resource density (x[0])
  C   consumer density (x[1])
  P   predator density (x[2])

Sampled parameters:

  x_c   consumer metabolic and per-capita mortality scale (e[0])
  y_c   maximum consumer ingestion rate relative to metabolism (e[1])
  R_0   half-saturation constant for consumer resource uptake (e[2])
  x_p   predator metabolic, solitary-production and mortality scale (e[3])
  x_i   coefficient of density-dependent cooperative predator production (e[4])
  y_p   maximum predator ingestion rate relative to metabolism (e[5])
  C_0   half-saturation constant for predator consumption (e[6])

Control parameters:

  D       fraction of habitat loss
  sigma   fraction of predators cooperating in hunting

The cooperative/solitary predator term implemented in both programs is

  Phi = x_p (1 - sigma) P + x_i sigma P^2.


NUMERICAL METHOD
----------------

The three ordinary differential equations are integrated using an embedded
seventh/eighth-order Runge-Kutta-Fehlberg method with adaptive time stepping.
The supplied bounds on the time step are HMIN = 1e-5 and HMAX = 0.1.


OUTCOME CLASSIFICATION
----------------------

The programs classify simulations into four ecological outcomes:

  RCP   persistence of the resource, consumer and predator
  RC    persistence of the resource and consumer, with predator loss
  R     persistence of the resource only
  0     extinction of all three populations

Population densities below the numerical cutoff are set to zero, making
extinction irreversible within a simulation. This represents effective
extinction at densities that would be highly vulnerable to stochastic loss in
a natural population.

The two supplied source files were produced at different stages of the
analysis and currently use slightly different numerical cutoffs and sampling
ranges. Therefore, the constants and ranges in the selected source file should
be checked against those reported in the final manuscript before running the
full analysis.


COMPILATION
-----------

A C compiler and the standard mathematics library are required. With GCC:

  gcc -O3 -Wall -o sensitivity-break sensitivity-break.c -lm

  gcc -O3 -Wall -o sensitivity-break_fixed_seed \
      sensitivity-break_fixed_seed.c -lm

The programs use drand48(), getpid() and other POSIX/GNU functions and are
therefore intended for Linux, macOS or another POSIX-compatible environment.


RUNNING
-------

Run the parameter scan with:

  ./sensitivity-break

Run the reproducible single-combination analysis with:

  ./sensitivity-break_fixed_seed

To preserve the console output from either program, redirect it to a file, for
example:

  ./sensitivity-break_fixed_seed > sensitivity_fixed_seed_results.txt

The requested sample sizes are computationally demanding. For preliminary
tests, reduce ITERS and norm by the same factor, then restore their reported
values for the final analysis.


REPRODUCIBILITY
---------------

For an exact reproducible run, record:

  - the source-file version;
  - the values of D, sigma, ITERS, norm and time_f;
  - all parameter and initial-condition sampling ranges;
  - the extinction cutoff;
  - the compiler and compiler version;
  - the operating system; and
  - the random seed or deterministic initial state.

For archival publication, preserve the exact source files used to generate the
reported results together with their raw output files and analysis scripts.
