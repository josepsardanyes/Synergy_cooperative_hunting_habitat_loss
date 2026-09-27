README
======

Program
-------
Unimodal_Map_Supplementary_Material.nb

Purpose
-------
This Wolfram Mathematica notebook numerically constructs the one-dimensional
return map used in Figure 3 from the dynamics of a three-level
resource-consumer-predator model.

The central idea is to integrate the continuous-time ecological model, extract
successive local minima of the predator population P(t), and plot each minimum
against the following one,

    P_(n+1)^min versus P_n^min.

The resulting set of points provides a numerical first-return map. For the
parameter values used in the notebook, this map has the unimodal shape employed
in the symbolic-dynamics analysis associated with Figure 3.

Software
--------
The notebook was created with Wolfram Mathematica 15.0.

No external Mathematica packages or external data files are required. All
calculations use standard built-in Mathematica functions.

Model
-----
The state variables are

    R(t) : resource population,
    C(t) : consumer population,
    P(t) : predator population.

In the Mathematica code these variables are denoted by RR[t], CC[t], and PP[t].

The predator function is defined by

    F(P) = P [(1 - sigma) xp + sigma xi P].

The differential-equation system implemented in the notebook is

    dR/dt = R (1 - D/k - R/k)
            - xc yc C R/(R + R0),

    dC/dt = xc C [yc R/(R + R0) - 1]
            - F(P) yp C/(C + C0),

    dP/dt = F(P) yp C/(C + C0)
            - xp P.

The initial conditions are

    R(0) = 0.55,
    C(0) = 0.35,
    P(0) = 0.8.

Parameters used in the supplied calculation
-------------------------------------------
The parameter values explicitly assigned inside the main calculation are

    C0    = 0.50
    xc    = 0.4
    yc    = 2.009
    xp    = 0.08
    yp    = 2.876
    R0    = 0.16129
    xi    = 0.08
    k     = 0.9995
    sigma = 0.1
    D     = 0.015
    tmax  = 5000

The Mathematica variable DD represents the parameter D.

The notebook contains outer For loops for C0, k, and DD. In the supplied
version, however, their upper and lower limits are equal, so each loop executes
only once:

    C0 : 0.50
    k  : initially 1.0 in the outer loop
    DD : 0.015

Inside the DD loop, k is subsequently set explicitly to 0.9995 and C0 to 0.50.
Consequently, the actual numerical integration shown in the notebook uses

    k = 0.9995, C0 = 0.50, D = 0.015, sigma = 0.1.

The increments appearing in the loop definitions would become relevant only if
the corresponding upper limits were enlarged.

How the program works
---------------------

1. Clear previous definitions
   --------------------------
   Clear[...] removes existing definitions of the principal model variables and
   parameters before the calculation begins.

2. Initialize the parameter loops
   ------------------------------
   The program enters nested For loops for C0, k, and DD. In the present
   supplementary-material version these loops select a single parameter set.

   At the beginning of the DD calculation, the program prints

       DD = 0.015

   and, after assigning sigma, prints

       sigma = 0.1.

3. Define the predator function
   ----------------------------
   Mathematica defines

       FF[PP_] := PP ((1 - sigma) xp + sigma xi PP).

   This function enters both the consumer and predator differential equations
   and incorporates the cooperation parameter sigma.

4. Define the ODE system
   ---------------------
   The three differential equations and their initial conditions are collected
   in the Mathematica variable

       hindequat.

5. Numerically integrate the model
   --------------------------------
   NDSolve solves the system over

       0 <= t <= 5000.

   The relevant solver settings are

       WorkingPrecision -> MachinePrecision
       Method           -> StiffnessSwitching
       MaxSteps         -> 200000

   The numerical solution is stored in

       sol1.

   StiffnessSwitching allows Mathematica to switch numerical integration
   strategies if stiffness is detected during the trajectory.

6. Sample the numerical trajectory
   --------------------------------
   The NDSolve solution is evaluated at integer times

       t = 0, 1, 2, ..., 5000.

   Each sampled entry has the form

       {R(t), C(t), P(t)}.

   The complete sampled trajectory is stored in

       tab66.

7. Extract the predator time series
   --------------------------------
   The program takes the third component of every entry in tab66 and stores it
   in

       listaP.

   Thus,

       listaP = {P(0), P(1), ..., P(5000)}

   evaluated from the numerical solution.

   Notice that, in this particular notebook, no transient interval is explicitly
   removed before the minima are detected. All sampled predator values are
   passed to the local-minimum detection step.

8. Detect successive local minima of P(t)
   ---------------------------------------
   The code scans listaP and identifies a sampled point P_(i+1) as a local
   minimum whenever

       P_(i+1) - P_i     < 0

   and

       P_(i+2) - P_(i+1) > 0.

   In other words, the sampled trajectory must be decreasing immediately before
   the point and increasing immediately after it.

   Each detected minimum is appended to

       listmin.

   The command

       tabz = listmin

   then creates the sequence

       tabz = {P_1^min, P_2^min, P_3^min, ...}.

9. Construct the first-return map
   ------------------------------
   Consecutive minima are paired according to

       {P_n^min, P_(n+1)^min}.

   Mathematica performs this operation for all successive elements of tabz and
   stores the resulting points in

       listamaptotal.

   Hence,

       listamaptotal =
       {
         {P_1^min, P_2^min},
         {P_2^min, P_3^min},
         {P_3^min, P_4^min},
         ...
       }.

   This is the numerical first-return map generated from the minima of the
   predator variable.

10. Plot the return map
    -------------------
    ListPlot is used to display listamaptotal.

    The horizontal coordinate represents

        P_n^min

    and the vertical coordinate represents

        P_(n+1)^min.

    The principal plotting options are

        AspectRatio -> 1
        Axes        -> False
        Frame       -> True
        AbsolutePointSize[3]
        PlotRange   -> {{0.58, 0.78}, {0.58, 0.78}}
        PlotStyle   -> black.

    The resulting graphic is stored in

        gtotal.

11. Construct the diagonal P_(n+1) = P_n
    -------------------------------------
    The notebook also defines the identity function y = x through

        gdiag = Plot[x, {x, 0.55, 0.75}, ...].

    In the final displayed figure, however, the diagonal is inserted explicitly
    using Graphics and Line rather than by directly combining gdiag with the
    return-map plot.

12. Create the final Figure 3 return-map graphic
    ---------------------------------------------
    Show combines gtotal with a dashed diagonal line extending from

        (0.55, 0.55)

    to

        (0.78, 0.78).

    The diagonal therefore represents

        P_(n+1)^min = P_n^min.

    The frame ticks are specified at approximately

        0.55, 0.60, 0.65, 0.70, 0.75,

    and the frame is formatted in black with font size 17.

    The complete graphic is stored in

        rrrr

    and displayed with

        Print[rrrr].

Interpretation of the return map
--------------------------------
The numerical points in the final graph represent the transformation

    P_n^min  ->  P_(n+1)^min,

where successive local minima are taken from the predator trajectory.

A concentration of the points near a one-dimensional curve indicates that the
sequence of predator minima can be approximately described by an iterated map

    P_(n+1)^min = f(P_n^min).

For the parameter set used here, the curve has the unimodal structure relevant
to the analysis presented in Figure 3. Such a return map provides the
one-dimensional representation subsequently used for symbolic-dynamics and
related dynamical analyses.

The dashed identity line

    P_(n+1)^min = P_n^min

is useful for visually identifying fixed points of the return map: intersections
between the numerical map and this diagonal correspond to values satisfying

    f(P*) = P*.

Important Mathematica variables
-------------------------------
    DD             model parameter D
    sigma          cooperation parameter
    FF             predator-dependent function F(P)
    hindequat      differential equations plus initial conditions
    sol1           numerical NDSolve solution
    tab66          sampled {R,C,P} trajectory
    listaP         sampled predator time series
    listmin        detected local minima of P(t)
    tabz           sequence of predator minima
    listamaptotal  pairs {P_n^min,P_(n+1)^min}
    gtotal         ListPlot of the numerical return map
    gdiag          separately defined plot of y = x
    rrrr           final return-map figure including the dashed diagonal

How to run the notebook
-----------------------
1. Open

       Unimodal_Map_Figure_3_Supplementary_Material.nb

   in Wolfram Mathematica.

2. Select the main input cell and evaluate it with Shift+Enter, or use
   Evaluation -> Evaluate Notebook.

3. Mathematica numerically integrates the ODE system, samples the solution,
   detects the local minima of P(t), constructs the return map, and displays
   the final figure.

4. The notebook also prints the values of DD and sigma used in the calculation.

No external input files are needed.

No external output file is automatically created; the figure appears directly
as Mathematica notebook output.

Changing the parameters
-----------------------
The calculation can be adapted by changing, for example,

    DD
    sigma
    k
    C0
    xc, yc, xp, yp, R0, xi
    the initial conditions
    tmax
    the sampling interval
    the PlotRange.

To perform an actual parameter scan, change the upper limit of the relevant
For loop so that it exceeds its initial value. The existing increments are

    C0 increment : 0.002
    k increment  : 0.00055
    DD increment : 0.000464914.

Care should be taken with k because the current program sets k = 0.9995 inside
the DD loop. This internal assignment overrides the value supplied by the
outer k loop. If a genuine scan over k is intended, that internal assignment
should be removed or modified.

Numerical considerations
------------------------
1. Discrete sampling of minima
   The NDSolve interpolating solution is sampled only at integer values of time.
   Therefore, listmin contains local minima of the discretely sampled series,
   not the exact continuous-time minima of the interpolating solution.

2. Sampling resolution
   If the oscillations become sufficiently fast relative to the unit sampling
   interval, the detected positions and values of the minima can differ from
   their continuous-time values. A smaller sampling step can be used if higher
   temporal resolution is required.

3. Transients
   The supplied version uses the entire sampled trajectory, including the
   initial transient. If the purpose is to reconstruct only the asymptotic
   attractor, a transient interval may be discarded before detecting the local
   minima.

4. Number of minima
   The return map requires at least two detected local minima. Parameter sets
   approaching an equilibrium, extinction state, or weakly oscillatory regime
   may yield too few minima for a meaningful return map.

5. Plot range
   The final plot is restricted to

       0.58 <= P_n^min     <= 0.78
       0.58 <= P_(n+1)^min <= 0.78.

   If parameters are changed substantially, the return-map points may fall
   outside this window and the PlotRange should then be adjusted.

6. Efficiency
   The notebook uses AppendTo inside For loops. This is straightforward and
   transparent for supplementary code, although Mathematica operations such as
   Part, Differences, Select, Partition, or Reap/Sow can be more efficient for
   very large parameter scans.

Relation to the symbolic-dynamics calculation
----------------------------------------------
This notebook constructs the geometric object underlying the symbolic analysis:
the unimodal first-return map based on consecutive predator minima.

The associated symbolic-dynamics supplementary program can then use the
critical point of this return map to partition its domain and encode successive
iterates by symbols according to their positions relative to that critical
point.

Thus, conceptually, the workflow is

    continuous three-species ODE model
              |
              v
       predator trajectory P(t)
              |
              v
      successive minima P_n^min
              |
              v
    first-return pairs
    (P_n^min, P_(n+1)^min)
              |
              v
       unimodal return map
              |
              v
       symbolic dynamics.

Summary
-------
The notebook integrates a three-level ecological model for the parameter set
used in Figure 3 and extracts successive local minima of the predator
population. Consecutive minima are paired to reconstruct a numerical
one-dimensional first-return map. Mathematica then plots this map together with
the identity line. The resulting unimodal map provides a reduced description of
the oscillatory predator dynamics and forms the basis for the corresponding
symbolic-dynamics analysis.
