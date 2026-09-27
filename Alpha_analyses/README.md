# Possible-coexistence scans

These two C programs map where the resource (`R`), consumer (`C`), and predator (`P`) can coexist in the cooperative-hunting model. They use an adaptive RKF78 integrator, OpenMP parallelization, continuation from the preceding parameter slice, and iterative neighbour seeding to detect coexistence that may depend on initial conditions.

## Programs

- `coexist_RKF78_alphas_possible_space.c` scans one alpha slice at a time over `D` and `sigma`. Grids: `alpha = 0...5` by `0.1`, `D = 0...0.5` by `0.01`, and 201 `sigma` values from `0.01` to `1`. It writes `results_alpha_<alpha>_refined_200.txt`.
- `coexist_RKF78_Ds_possible_space.c` scans one habitat-loss slice at a time over `alpha` and `sigma`. Grids: `D = 0...0.5` by `0.01`, `alpha = 0...20` by `0.2`, and 101 `sigma` values from `0.01` to `1`. It writes `results_D_<D>_refined.txt`.

Output columns are:

- alpha-first program: `alpha D sigma coexistence R C P`
- D-first program: `D alpha sigma coexistence R C P`

`coexistence` is `1` when all three species persist and `0` when at least one becomes extinct. Extinction is declared below `1e-6`. Integrations can continue to `t = 200000`, so complete scans may take a long time. Output files are written to the current working directory and existing files with the same names are replaced.

## Compile and run

```bash
gcc -O3 -std=c11 -Wall -Wextra -pedantic -fopenmp -o coexist_RKF78_alphas_possible_space coexist_RKF78_alphas_possible_space.c -lm
gcc -O3 -std=c11 -Wall -Wextra -pedantic -fopenmp -o coexist_RKF78_Ds_possible_space coexist_RKF78_Ds_possible_space.c -lm

./coexist_RKF78_alphas_possible_space
./coexist_RKF78_Ds_possible_space
```

No command-line arguments are required. Set `OMP_NUM_THREADS` if you want to limit CPU use. The grid limits are named constants near the beginning of each `main` function.

## Numerical specification

The fixed parameters are `K=1`, `xc=0.4`, `yc=2.099`, `yp=2.876`, `xp=0.08`, `xi_base=0.08`, `R0=0.16129`, and `C0=0.5`; the default initial state is `(R,C,P)=(0.55,0.35,0.8)`. The adaptive RKF78 solver uses tolerance `1e-8`, minimum step `1e-5`, and maximum step `0.005` during these scans.

The alpha-first source originally combined `alpha_max=5` with a hard-coded count of 201, unintentionally extending the scan to 20. This release derives the count from `alpha_min`, `alpha_max`, and `alpha_step`, making the documented `0...5` range effective.
