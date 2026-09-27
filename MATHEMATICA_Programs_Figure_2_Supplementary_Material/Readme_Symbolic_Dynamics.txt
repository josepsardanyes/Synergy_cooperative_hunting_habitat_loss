README
======

File
----
Symbolic_Dynamics_Supplementary_Material.nb

Purpose
-------
This Wolfram Mathematica notebook computes symbolic sequences associated with the
dynamics of a three-level resource-consumer-predator model. The procedure is
intended as supplementary material for the symbolic-dynamics calculations used
in Figure 3.

For each selected value of the parameter D, the notebook:

1. numerically integrates the continuous-time ecological model;
2. removes an initial transient;
3. extracts successive local minima of the predator variable P(t);
4. constructs a one-dimensional first-return map from these minima;
5. identifies the critical point of the return map; and
6. encodes the forward orbit relative to that critical point using the symbols
   L, A, and R.

The resulting symbolic sequence and the corresponding value of D are printed
to the Mathematica output.

Software
--------
The notebook was created with Wolfram Mathematica 15.0 for Microsoft Windows
(64-bit). It uses only standard built-in Mathematica functions and requires no
external packages or data files.

Model
-----
The state variables are

    RR(t) : resource,
    CC(t) : consumer,
    PP(t) : predator.

The predator-dependent function used in the model is

    FF(P) = P [(1 - sigma) xp + sigma xi P].

The differential equations implemented in the notebook are

    R' = R (1 - D/k - R/k)
         - xc yc C R/(R + R0),

    C' = xc C [yc R/(R + R0) - 1]
         - FF(P) yp C/(C + C0),

    P' = FF(P) yp C/(C + C0)
         - xp P.

The initial conditions are

    R(0) = 0.55,
    C(0) = 0.35,
    P(0) = 0.8.

Parameter values
----------------
The parameter values assigned inside the D loop are

    xc = 0.4
    yc = 2.009
    xp = 0.08
    yp = 2.876
    R0 = 0.16129
    xi = 0.08
    C0 = 0.50

The current notebook settings also use

    sigma = 0.01
    k     = 1.0
    tmax  = 5000

and scan

    D = 0.01, 0.0105, 0.0110, ..., 0.06.

Note that the outer loops are written in a general form, but with the present
upper limits they execute only once for sigma = 0.01 and once for k = 1.0.
Their increments (0.01 and 0.00055, respectively) would become relevant if the
corresponding upper limits were increased.

How the computation works
--------------------------

1. Initialization
   The notebook clears the principal variables and starts the parameter loops.
   The current calculation uses sigma = 0.01 and k = 1.0.

2. Scan over D
   For each D between 0.01 and 0.06 in steps of 0.0005, the model parameters
   and initial conditions are assigned.

3. Numerical integration
   NDSolve integrates the three coupled ODEs on

       0 <= t <= 5000.

   The solver options are

       WorkingPrecision -> MachinePrecision
       Method           -> StiffnessSwitching
       MaxSteps         -> 200000

   The numerical solution is then sampled at integer times t = 0,1,...,5000,
   producing the table

       tab66 = {R(t), C(t), P(t)}.

4. Removal of the transient
   Only predator values corresponding to entries with index i > 3000 are
   retained. Since the sampled table starts at t = 0, this keeps the sampled
   predator dynamics from approximately t = 3000 onward. The resulting series
   is stored in listaP.

5. Detection of local minima
   A sampled value P[i+1] is classified as a local minimum when

       P[i+1] - P[i]   < 0

   and

       P[i+2] - P[i+1] > 0.

   Thus, the code detects changes from decreasing to increasing values in the
   discretely sampled predator time series. The detected minima are stored in

       listmin

   and then copied to

       tabz.

6. Construction of the first-return map
   Consecutive predator minima are paired as

       {P_n, P_(n+1)}.

   These pairs are stored in

       listamap.

   Therefore, listamap represents the numerical first-return map of successive
   local minima. The list

       tabzy

   contains the second coordinate P_(n+1) of every return-map pair.

7. Identification of the critical point
   The code finds

       zymax = Max[tabzy]

   and the position iM at which this maximum occurs. The corresponding
   return-map point is

       pxy1 = listamap[[iM]].

   The first coordinate of this point,

       xc1 = listamap[[iM,1]],

   is used as the critical value c1 that partitions the domain of the return
   map into left and right sides.

   In other words, the critical point is numerically identified through the
   return-map pair having the largest second coordinate.

8. Symbolic encoding
   The symbolic string is initialized as

       S1 = "A",

   where A denotes the critical point/critical symbol.

   The program then follows 36 entries of the return-map orbit, beginning at
   index ic1 = iM. Each value xxx is compared with

       c1 = xc1.

   The symbols are assigned according to

       xxx < c1  ->  "L"
       xxx = c1  ->  "A"
       xxx > c1  ->  "R"

   and concatenated to S1.

   Consequently, the output is an initial "A" followed by 36 symbols generated
   from the orbit relative to the critical partition.

9. Output
   For every value of D, Mathematica prints a line of the form

       <symbolic sequence> D=<value>

   In addition, the notebook prints the current values of sigma and k at the
   beginning of their respective loops.

Important intermediate variables
--------------------------------
    hindequat : system of ODEs and initial conditions
    sol1      : numerical NDSolve solution
    tab66     : sampled {R,C,P} trajectory
    listaP    : post-transient predator time series
    listmin   : detected local minima of P
    tabz      : copy of the local-minimum sequence
    listamap  : return-map pairs {P_n,P_(n+1)}
    tabzy     : second coordinates of the return-map pairs
    zymax     : maximum value in tabzy
    iM        : index associated with zymax
    pxy1      : return-map point at index iM
    xc1       : first coordinate of pxy1; critical partition value
    c1        : symbolic-dynamics partition value (= xc1)
    S1        : symbolic sequence

The notebook also constructs tabzm and tabzym, containing return-map points
whose first coordinate is larger than xc1 and their second coordinates,
respectively. These arrays are retained by the program but are not subsequently
used in the final symbolic-string calculation in the present version.

Running the notebook
--------------------
1. Open Symbolic_Dynamics_Figure_3_Supplementary_Material.nb in Wolfram
   Mathematica.
2. Evaluate the input cell (for example, press Shift+Enter with the cell
   selected, or use "Evaluate Notebook").
3. Wait for the integrations and the complete D scan to finish.
4. Read the printed symbolic sequence associated with each D value in the
   notebook output.

No input files are required and the program does not write external output
files. Results are printed directly in the Mathematica notebook.

Changing the calculation
------------------------
The main quantities that can be modified directly in the notebook are:

    sigma       cooperation parameter used by FF(P)
    k           parameter appearing in the resource equation
    DD          scanned D parameter
    tmax        final integration time
    model parameters xc, yc, xp, yp, R0, xi, C0
    initial conditions R(0), C(0), P(0)
    transient cutoff (currently i > 3000)
    symbolic itinerary length (currently 36 iterations)

To scan wider ranges of sigma or k, increase the corresponding upper bound in
the outer For loops. To alter the D scan, change its initial value, upper bound,
or increment.

Numerical considerations
------------------------
* The trajectory is sampled only at integer times before local minima are
  detected. Therefore, the minima used for the return map are minima of the
  sampled time series, rather than continuously optimized minima of the
  interpolating NDSolve solution.

* A long transient is discarded before constructing the return map. This is
  important because symbolic dynamics should characterize the asymptotic
  attractor rather than the initial approach to it.

* The symbolic encoding requires sufficiently many detected minima after the
  critical-map index to construct all 36 requested symbols. Parameter changes
  that produce too few minima may therefore lead to Part/index errors.

* If the dynamics approaches extinction, an equilibrium, or another regime
  with too few oscillations, listmin/tabzy can become too short or empty. In
  that case Max[tabzy] and subsequent indexing operations will not be defined.

* The test xxx == c1 uses exact equality between stored machine-precision
  numerical values. In the present construction the critical value is taken
  directly from the stored return-map data, but users modifying the algorithm
  should keep machine-precision comparisons in mind.

Summary
-------
The notebook converts the long-term predator oscillations of the continuous
three-species model into a one-dimensional return map based on successive
predator minima. It then locates the critical point of this numerical map and
records a finite symbolic itinerary using the alphabet {L,A,R}. Repeating this
procedure across the specified D interval produces the symbolic-dynamics
information associated with Figure 3.
