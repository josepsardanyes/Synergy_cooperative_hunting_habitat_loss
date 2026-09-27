README
======

File
----
Variation_3_Lyapunov_Exponents_Supplementary_Material.nb

Purpose
-------
This Wolfram Mathematica notebook computes the three Lyapunov exponents of the
three-dimensional resource-consumer-predator model used in Figure 3, while the
habitat-loss/disturbance parameter D is varied.

The calculation is based on the variational equations associated with the
nonlinear ODE system. The tangent dynamics are integrated together with the
original trajectory and the tangent vectors are repeatedly orthogonalized and
renormalized. The accumulated logarithmic stretching factors provide numerical
estimates of the Lyapunov exponents.

The final result consists of three data sets,

    {D, lambda_1}, {D, lambda_2}, {D, lambda_3},

and three corresponding plots showing how the Lyapunov spectrum changes with D.

Software
--------
The notebook was created with Wolfram Mathematica 15.0.

It uses standard Mathematica functionality, including NDSolve, differentiation,
vector operations, Projection, Norm, and ListPlot. No external data files are
required.

Model
-----
Let

    x1 = R : resource,
    x2 = C : consumer,
    x3 = P : predator.

The vector field F = (F1,F2,F3) implemented in the notebook is

    F1 = R (1 - D/k - R/k)
         - xc yc C R/(R + R0),

    F2 = xc C [yc R/(R + R0) - 1]
         - P [(1-sigma) xp + sigma xi P] yp C/(C + C0),

    F3 = P [(1-sigma) xp + sigma xi P] yp C/(C + C0)
         - xp P.

Thus the ecological model is

    dX/dt = F(X),

with

    X = (R,C,P).

Parameter values
----------------
The fixed parameters in the supplied notebook are

    xc    = 0.4
    yc    = 2.009
    xp    = 0.08
    yp    = 2.876
    R0    = 0.16129
    xi    = 0.08
    C0    = 0.50
    sigma = 0.1
    k     = 1

The initial condition for the nonlinear trajectory is

    R(0) = 0.55
    C(0) = 0.35
    P(0) = 0.8.

The parameter D is scanned over

    D = 0.0, 0.000125, 0.000250, ..., 0.04.

Hence the notebook evaluates the Lyapunov spectrum at 321 parameter values
(including both endpoints, subject to ordinary machine-number loop arithmetic).

Lyapunov-exponent method
------------------------
For a nonlinear autonomous system

    X' = F(X),

an infinitesimal perturbation delta X evolves according to the variational
equation

    delta X' = J(X(t)) delta X,

where

    J = DF

is the Jacobian matrix of the vector field.

Because the present system has dimension three, three independent tangent
vectors are evolved. Their asymptotic average rates of expansion or contraction
give the three Lyapunov exponents.

If the stretching factors accumulated during successive renormalization
intervals are s_j^(n), the numerical estimate after N intervals has the form

    lambda_j(N)
       = [sum_(n=1)^N log(s_j^(n))]/(N Delta t),

where Delta t is the duration of each integration interval.

In this notebook,

    Delta t = tmax = 100

and

    K = 1000

renormalization intervals are used for every value of D. Thus the accumulated
Lyapunov estimate corresponds to a total integration time of up to

    K tmax = 100000

for each D value.

How the program works
---------------------

1. Initialization
   --------------
   Clear[...] removes definitions left from previous Mathematica evaluations.

   The ecological parameters are assigned and the output lists are initialized:

       ListaLyap1 = {}
       ListaLyap2 = {}
       ListaLyap3 = {}.

   These lists will eventually contain the three Lyapunov exponents as functions
   of D.

2. Parameter scan
   ---------------
   The main For loop varies DD (the Mathematica variable representing D) from

       0.0 to 0.04

   with increment

       0.000125.

   For each parameter value, the notebook prints

       DD=<current value>

   so that the progress of the calculation can be followed.

3. Definition of the nonlinear vector field
   -----------------------------------------
   The function

       F[{x1_,x2_,x3_}]

   defines the right-hand side of the three ecological differential equations.

   Here x1, x2, and x3 correspond respectively to R, C, and P.

4. Construction of the Jacobian
   ----------------------------
   The helper function

       JacobianMatrix[funs_List, vars_List] :=
           Outer[D, funs, vars]

   differentiates the components of the vector field with respect to the three
   state variables.

   The resulting 3 x 3 matrix J is the Jacobian

       J(X) = DF(X).

   This matrix determines the linearized evolution of infinitesimal
   perturbations around the current nonlinear trajectory.

5. Nonlinear and variational variables
   ------------------------------------
   The Mathematica variables y_1(t), y_2(t), and y_3(t) represent the nonlinear
   state variables.

   Nine additional variables represent three tangent vectors. Therefore,
   NDSolve simultaneously integrates a total of 12 scalar variables:

       3 nonlinear state variables
       +
       9 variational variables.

   The three tangent vectors can be regarded as the columns/rows of a 3 x 3
   fundamental perturbation matrix.

6. Original differential equations
   --------------------------------
   The list EQ3 contains the three equations

       y_i'(t) = F_i(y_1(t),y_2(t),y_3(t)),
       i = 1,2,3.

7. Variational equations
   ---------------------
   The tangent matrix Y is formed from the nine perturbation variables.

   Multiplication of the Jacobian J by the tangent matrix gives the linearized
   dynamics. The corresponding nine scalar equations are collected in EQ9.

   In matrix notation this part of the calculation is

       Y'(t) = J(X(t)) Y(t).

8. Initial tangent vectors
   -----------------------
   YI9 initializes the three perturbation directions as mutually independent
   unit directions (an identity-type tangent basis).

   The nonlinear initial conditions are stored in YI3:

       y_1(0) = 0.55
       y_2(0) = 0.35
       y_3(0) = 0.8.

9. Repeated integration
   --------------------
   The program sets

       K    = 1000
       tmax = 100.

   A Do loop is then performed K times.

   During each iteration NDSolve integrates the combined nonlinear and
   variational system over

       0 <= t <= 100

   using

       WorkingPrecision -> MachinePrecision
       Method           -> StiffnessSwitching
       MaxSteps         -> 80000.

10. Recover the tangent vectors
    ---------------------------
    At t = tmax, the evolved tangent variables are extracted from the numerical
    solution and assembled into the matrix PhiT.

    PhiT therefore describes how the tangent basis has been deformed during the
    current integration interval.

11. Gram-Schmidt orthogonalization
    -------------------------------
    If tangent vectors were simply evolved indefinitely, they would tend to
    align with the most unstable direction. This would prevent the three
    Lyapunov exponents from being determined independently.

    The notebook therefore applies a Gram-Schmidt procedure after every
    integration interval.

    A direct call to GramSchmidt is present in the code as a commented-out
    alternative. The active implementation constructs the orthogonal vectors
    explicitly using Mathematica's Projection function.

    Schematically,

       w1 = phi1

       w2 = phi2 - Proj_w1(phi2)

       w3 = phi3 - Proj_w1(phi3) - Proj_w2(phi3).

    This produces three mutually orthogonal tangent directions.

12. Measure tangent-vector stretching
    ---------------------------------
    The command

       norms = Map[Norm,W]

    computes the lengths of the three orthogonalized vectors.

    These lengths are the stretching/contraction factors accumulated during the
    current interval.

    They are appended to

       s.

    Thus s stores the successive three-component vectors of stretching factors.

13. Renormalize the tangent basis
    -----------------------------
    The command

       PhiT = W/norms

    normalizes the orthogonal vectors.

    The normalized tangent basis is then used to define YI9 for the next
    integration interval.

    This repeated

       integrate -> orthogonalize -> measure -> normalize

    cycle is the core of the Lyapunov-exponent algorithm.

14. Continue the nonlinear trajectory
    ---------------------------------
    The endpoint of the nonlinear solution from the current interval becomes
    the initial state YI3 for the next interval.

    Consequently, the trajectory is not restarted at (0.55,0.35,0.8) on every
    Gram-Schmidt step. Instead, the integration proceeds continuously through
    phase space in consecutive blocks of length 100.

15. Compute finite-time Lyapunov estimates
    --------------------------------------
    After all K iterations, the code evaluates

       lce =
         Rest[FoldList[Plus,0,Log[s]]] /
         (tmax Range[K]).

    Log[s] contains the logarithms of the stretching factors.

    FoldList forms their cumulative sums, and division by

       tmax {1,2,...,K}

    converts the accumulated logarithmic stretching into average rates per unit
    time.

    Therefore lce[[n]] contains the three finite-time Lyapunov estimates after
    n integration/renormalization intervals.

16. Separate the three convergence histories
    -----------------------------------------
    The program creates

       ListaIndLyap1
       ListaIndLyap2
       ListaIndLyap3

    by extracting each component of lce.

    These lists represent the convergence histories of the individual Lyapunov
    estimates for the current D value.

17. Store the final exponents for the current D
    -------------------------------------------
    Only the last cumulative estimate is used in the D-dependent spectrum.

    The notebook appends

       {DD, Last[lce][[1]]}  to ListaLyap1
       {DD, Last[lce][[2]]}  to ListaLyap2
       {DD, Last[lce][[3]]}  to ListaLyap3.

    The entire calculation is then repeated for the next D value.

18. Plot the Lyapunov spectrum
    --------------------------
    After the D scan has finished, three joined ListPlot graphics are created:

       GLyap1 : lambda_1 versus D
       GLyap2 : lambda_2 versus D
       GLyap3 : lambda_3 versus D.

    In the supplied notebook they are displayed using red, green, and blue,
    respectively. Each graph has

       AspectRatio -> 0.3
       Frame       -> True
       PlotRange   -> All.

Interpretation
--------------
For an autonomous three-dimensional continuous-time dynamical system, the
Lyapunov exponents quantify the asymptotic sensitivity of nearby trajectories.

In general:

    lambda > 0
        indicates exponential separation in the corresponding tangent
        direction;

    lambda < 0
        indicates exponential contraction;

    lambda approximately 0
        is expected along the flow direction for a regular non-equilibrium
        autonomous trajectory, up to finite-time and numerical errors.

A positive maximum Lyapunov exponent is the standard numerical signature of
sensitive dependence on initial conditions and is therefore used as evidence
for chaotic dynamics.

The complete three-exponent spectrum also contains information that is not
available from the maximum exponent alone, including contraction in the
remaining phase-space directions.

Important Mathematica variables
-------------------------------
    DD              scanned parameter D
    F               nonlinear vector field
    J               3 x 3 Jacobian matrix DF
    Y               tangent/variational matrix
    EQ3             three nonlinear differential equations
    EQ9             nine variational differential equations
    YI3             nonlinear initial conditions for each block
    YI9             tangent-vector initial conditions for each block
    sol             NDSolve solution for one integration block
    PhiT            evolved tangent-vector matrix
    W               Gram-Schmidt orthogonalized tangent vectors
    norms           norms/stretching factors of W
    s               stretching factors accumulated over all blocks
    K               number of renormalization blocks (=1000)
    tmax            duration of each block (=100)
    lce             cumulative finite-time Lyapunov estimates
    ListaIndLyap1   convergence history of exponent 1 for current D
    ListaIndLyap2   convergence history of exponent 2 for current D
    ListaIndLyap3   convergence history of exponent 3 for current D
    ListaLyap1      final {D,lambda_1} data
    ListaLyap2      final {D,lambda_2} data
    ListaLyap3      final {D,lambda_3} data
    GLyap1          plot of lambda_1(D)
    GLyap2          plot of lambda_2(D)
    GLyap3          plot of lambda_3(D)

Some additional lists are initialized near the beginning of the notebook
(e.g. ListaLyapMax3DTotal, ListaLyapMax3DPositive, and related variables) but
are not required by the principal three-exponent calculation described above
in the supplied version.

How to run the notebook
-----------------------
1. Open

       Variation_3_Lyapunov_Exponents_Figure_3_Supplementary_Material.nb

   in Wolfram Mathematica.

2. Evaluate the main input cell (Shift+Enter), or choose Evaluate Notebook.

3. Mathematica scans all specified D values. The current value of D is printed
   during the calculation.

4. For every D, the nonlinear and variational equations are repeatedly
   integrated, the tangent vectors are orthogonalized and normalized, and the
   three Lyapunov estimates are accumulated.

5. After the scan, the variables

       ListaLyap1, ListaLyap2, ListaLyap3

   contain the numerical Lyapunov-spectrum data.

6. The graphics

       GLyap1, GLyap2, GLyap3

   contain the corresponding plots.

No external input data are required, and the notebook does not automatically
write its results to an external file.

Computational cost
------------------
The supplied settings are computationally demanding.

There are approximately 321 D values. For every D, the notebook performs 1000
NDSolve integrations, each over a time interval of length 100, for a system of
12 coupled equations.

Consequently, a complete evaluation can require substantial computation time.
The exact runtime depends strongly on the computer, Mathematica version, and
the stiffness encountered along the trajectories.

For preliminary tests, K or the number of D values can be reduced. Final
scientific calculations should use sufficiently long integrations to establish
convergence of the exponents.

Numerical considerations
------------------------
1. Convergence
   Lyapunov exponents are asymptotic quantities. The finite-time values should
   be checked for convergence as the total integration time increases.

2. Initial transient
   The supplied implementation starts accumulating logarithmic stretching from
   the initial condition. For studies where transient contamination is a
   concern, one may first integrate the nonlinear system for a transient time
   before beginning the Lyapunov accumulation.

3. Orthogonalization frequency
   The current renormalization interval is 100 time units. Changing this value
   can affect numerical conditioning. Very long intervals can cause tangent
   vectors to become extremely aligned or produce very large/small norms.

4. Machine precision
   The calculation uses MachinePrecision. Near parameter values where an
   exponent is very close to zero, numerical accuracy and convergence should be
   assessed carefully.

5. Ordering of exponents
   ListaLyap1, ListaLyap2, and ListaLyap3 correspond to the tangent directions
   produced by the implemented Gram-Schmidt evolution. Users requiring an
   explicitly ordered spectrum should verify or sort the converged exponents
   before interpreting the labels lambda_1 >= lambda_2 >= lambda_3.

6. Autonomous-flow zero exponent
   Away from equilibria, one exponent of an autonomous continuous-time system
   is theoretically zero along the flow direction. A small nonzero computed
   value is expected in finite-time numerical calculations.

7. Parameter-loop arithmetic
   DD is incremented as a machine number. If exact reproducibility of the
   parameter grid is important, it can be preferable to generate an explicit
   list of D values rather than repeatedly adding 0.000125.

Changing the calculation
------------------------
The main quantities that can be changed are

    sigma                 cooperation parameter
    k                     resource parameter
    C0                    half-saturation parameter
    xc, yc, xp, yp, R0,
    xi                    model parameters
    initial condition     {0.55,0.35,0.8}
    D range               0 <= D <= 0.04
    D increment           0.000125
    K                     number of Lyapunov accumulation blocks
    tmax                  duration of each block.

For a faster exploratory run, reduce K and/or use a coarser D grid. These
shortened calculations should not automatically be regarded as converged
Lyapunov estimates.

Relation to Figure 3
--------------------
This notebook supplies the Lyapunov-exponent component of the dynamical
analysis associated with Figure 3.

Together with the return-map and symbolic-dynamics calculations, it provides a
complementary characterization of the model:

    continuous ecological model
              |
              v
       nonlinear trajectory
              |
              +------------------------------+
              |                              |
              v                              v
       return-map analysis          variational equations
              |                              |
              v                              v
       symbolic dynamics            tangent-vector growth
                                             |
                                             v
                                  Lyapunov spectrum lambda_i(D).

The return map describes the organization of successive oscillations, while the
Lyapunov calculation quantifies the average exponential expansion and
contraction of nearby trajectories.

Summary
-------
For each D between 0 and 0.04, the notebook integrates the three-dimensional
ecological model together with its nine variational equations. Three tangent
vectors are repeatedly evolved, Gram-Schmidt orthogonalized, measured, and
renormalized. The cumulative logarithms of their stretching factors yield
finite-time estimates of the three Lyapunov exponents. The final estimates are
stored as functions of D and plotted separately, providing the Lyapunov-spectrum
information used in the analysis associated with Figure 3.
