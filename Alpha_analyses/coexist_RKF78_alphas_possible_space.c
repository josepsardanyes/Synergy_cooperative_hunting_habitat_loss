/* Possible-coexistence scan with alpha as the outer continuation parameter.
 * Build and usage instructions are in README.md.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define DIMAX 200
#define MAX(a,b) ((a)<(b) ? (b) : (a))
#define sgn(a) ((a)<0 ? -1 : 1 )

// RKF78 parameters
#define HMAX 0.01
#define HMIN 0.00001

// System dimension
#define DIM 3

static void *xmalloc(size_t bytes)
{
    void *ptr = malloc(bytes);
    if (ptr == NULL) {
        fprintf(stderr, "Error: unable to allocate %zu bytes.\n", bytes);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

// --------------------
// Model parameters
// --------------------
double K  = 1.0;
double xc = 0.4;
double xp = 0.08;
double xi_base = 0.08;   // factor per alpha
double yc = 2.099;
double yp = 2.876;
double R0 = 0.16129;
double C0 = 0.5;

// --------------------
// RKF78 coefficients
// --------------------
static double alfa[13] =
{0., 2./27., 1./9., 1./6., 5./12., 0.5, 5./6., 1./6.,
 2./3., 1./3., 1., 0., 1.};

static double beta[79] =
{0., 2./27., 1./36., 1./12., 1./24., 0., 1./8., 5./12.,
 0., -25./16., 25./16., 0.5e-1, 0., 0., 0.25, 0.2,
 -25./108., 0., 0., 125./108., -65./27., 125./54., 31./300., 0.,
 0., 0., 61./225., -2./9., 13./900., 2., 0., 0.,
 -53./6., 704./45., -107./9., 67./90., 3., -91./108., 0., 0.,
 23./108., -976./135., 311./54., -19./60., 17./6., -1./12., 2383./4100., 0.,
 0., -341./164., 4496./1025., -301./82., 2133./4100., 45./82., 45./164., 18./41.,
 3./205., 0., 0., 0., 0., -6./41., -3./205., -3./41.,
 3./41., 6./41., 0., -1777./4100., 0., 0., -341./164., 4496./1025.,
 -289./82., 2193./4100., 51./82., 33./164., 12./41., 0., 1.};

static double c7[11] =
{41./840., 0., 0., 0., 0., 34./105., 9./35., 9./35.,
 9./280., 9./280., 41./840.};

static double c8[13] =
{0., 0., 0., 0., 0., 34./105., 9./35., 9./35.,
 9./280., 9./280., 0., 41./840., 41./840.};

// --------------------
// RCP model
// --------------------
void model_RCP(double t, double e[], double x[], double dx[])
{
    (void)t;
    double D     = e[0];
    double sigma = e[1];
    double xi    = e[2];

    double R = x[0];
    double C = x[1];
    double P = x[2];

    dx[0] = R*(1 - D/K - R/K) - xc*yc*C*R/(R+R0);
    dx[1] = xc*C*(yc*R/(R+R0) - 1)
            - (yp*C)/(C0+C)*(xp*(1-sigma)*P + xi*sigma*P*P);
    dx[2] = (yp*C)/(C0+C)*(xp*(1-sigma)*P + xi*sigma*P*P) - xp*P;
}

// --------------------
// RKF78 integrator
// --------------------
void rk78(double *at, double e[], double x[], double *ah,
          double tol, double hmin, double hmax, int n,
          void (*deriv)(double, double[], double[], double[]))
{
    double tpon, tol1, err, nor, kh, beth, h1;
    double k[DIMAX][13], x7[DIMAX], x8[DIMAX], xpon[DIMAX], dx[DIMAX];
    int j, l, i, m;

    do {
        m = 0;
        for (i = 0; i < 13; i++) {
            tpon = *at + alfa[i] * (*ah);

            for (j = 0; j < n; j++) xpon[j] = x[j];

            for (l = 0; l < i; l++) {
                m++;
                beth = (*ah) * beta[m];
                for (j = 0; j < n; j++)
                    xpon[j] += beth * k[j][l];
            }

            deriv(tpon, e, xpon, dx);

            for (j = 0; j < n; j++)
                k[j][i] = dx[j];
        }

        err = nor = 0.0;

        for (j = 0; j < n; j++) {
            x7[j] = x8[j] = x[j];

            for (l = 0; l < 11; l++) {
                kh = (*ah) * k[j][l];
                x7[j] += kh * c7[l];
                x8[j] += kh * c8[l];
            }

            x8[j] += (*ah) * (c8[11] * k[j][11] + c8[12] * k[j][12]);

            err += fabs(x8[j] - x7[j]);
            nor += fabs(x8[j]);
        }

        err /= n;

        tol1 = tol * (1 + nor / 100);

        if (err < tol1)
            err = MAX(err, tol1 / 256);

        h1 = *ah;
        *ah *= 0.9 * pow(tol1 / err, 0.125);

        if (fabs(*ah) < hmin) *ah = hmin * sgn(*ah);
        if (fabs(*ah) > hmax) *ah = hmax * sgn(*ah);

    } while ((err >= tol1) && (fabs(*ah) > hmin));

    *at += h1;

    for (j = 0; j < n; j++)
        x[j] = x8[j];
}

// --------------------
// Integrate one point
// --------------------
int integrate_point(double D, double sigma, double xi,
                    double X0[3], double finalX[3])
{
    double t = 0.0;
    double h = 0.001;
    double tol = 1e-8;

    double e[3] = {D, sigma, xi};
    double X[3] = {X0[0], X0[1], X0[2]};

    double Tfinal = 200000.0;

    while (t < Tfinal &&
           X[0] > 1e-6 && X[1] > 1e-6 && X[2] > 1e-6)
    {
        rk78(&t, e, X, &h, tol, HMIN, 0.005, DIM, model_RCP);
    }

    finalX[0] = X[0];
    finalX[1] = X[1];
    finalX[2] = X[2];

    int extinct = (X[0] < 1e-6 || X[1] < 1e-6 || X[2] < 1e-6);

    return !extinct;
}

// --------------------
// Estructura per candidats
// --------------------
typedef struct {
    int idx;
    int a_i;
    int s_i;
    int priority;
} Candidate;

int cmp_priority(const void *a, const void *b)
{
    const Candidate *ca = (const Candidate *)a;
    const Candidate *cb = (const Candidate *)b;
    if (ca->priority < cb->priority) return 1;
    if (ca->priority > cb->priority) return -1;
    return 0;
}

// --------------------
// Mapa ASCII ANSI ULTRA-RÀPID
// --------------------
void print_ascii_map_color_fast(int Nalpha, int Nsigma, int *coex_arr, int iter, double D)
{
    const char *RED   = "\033[31m";
    const char *GREEN = "\033[32m";
    const char *BLUE  = "\033[34m";
    const char *RESET = "\033[0m";

    printf("%s\n=== MAPA ASCII (D=%.3f, iteració %d) ===%s\n", BLUE, D, iter, RESET);

    char *line = xmalloc((Nsigma * 8 + 64) * sizeof(char));

    for (int a_i = 0; a_i < Nalpha; a_i++) {

        int pos = 0;
        int last_state = -1;

        for (int s_i = 0; s_i < Nsigma; s_i++) {
            int idx = a_i * Nsigma + s_i;
            int state = coex_arr[idx];

            if (state != last_state) {
                if (state == 1)
                    pos += sprintf(line + pos, "%s", GREEN);
                else
                    pos += sprintf(line + pos, "%s", RED);
                last_state = state;
            }

            line[pos++] = (state == 1 ? '#' : '.');
        }

        pos += sprintf(line + pos, "%s", RESET);
        line[pos] = '\0';
        printf("%s\n", line);
    }

    printf("%s=== FI MAPA ASCII ===%s\n\n", BLUE, RESET);
    free(line);
}

// --------------------
// MAIN
// --------------------
// --------------------
// MAIN (seeding entre alpha, espai sigma vs D)
// --------------------
int main(void)
{
    const double alpha_min = 0.0;
    const double alpha_max = 5.0;
    const double alpha_step = 0.1;
    const int Nalpha = (int)((alpha_max - alpha_min) / alpha_step + 0.5) + 1;
    const int Nsigma = 201;

    const double Dmin = 0.0;
    const double Dmax = 0.5;
    const double Dstep = 0.01;
    const int ND = (int)((Dmax - Dmin) / Dstep + 0.5) + 1;

    const double Smin = 0.01, Smax = 1.0;
    double X0_default[3] = {0.55, 0.35, 0.8};

    // Arrays per seeding entre alpha
    double *X_prev   = xmalloc((size_t)ND * Nsigma * 3 * sizeof(double));
    int    *coex_prev = xmalloc((size_t)ND * Nsigma * sizeof(int));

    // Inicialitzar seeding
    for (int idx = 0; idx < ND * Nsigma; idx++) {
        coex_prev[idx] = 0;
        X_prev[3*idx + 0] = X0_default[0];
        X_prev[3*idx + 1] = X0_default[1];
        X_prev[3*idx + 2] = X0_default[2];
    }

    // -------------------------
    // Bucle gran: alpha
    // -------------------------
    for (int a_i = 0; a_i < Nalpha; a_i++) {

        double alpha = alpha_min + alpha_step * a_i;
        double xi = alpha * xi_base;

        printf("\n\n=============================\n");
        printf("   ALPHA = %.3f\n", alpha);
        printf("=============================\n\n");

        // Arrays per aquest alpha
        double *sigma_arr = xmalloc((size_t)ND * Nsigma * sizeof(double));
        double *D_arr     = xmalloc((size_t)ND * Nsigma * sizeof(double));
        int    *coex_arr  = xmalloc((size_t)ND * Nsigma * sizeof(int));
        double *X_arr     = xmalloc((size_t)ND * Nsigma * 3 * sizeof(double));

        // -------------------------
        // FASE 1: malla amb seeding entre alpha
        // -------------------------
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int d_i = 0; d_i < ND; d_i++) {
            for (int s_i = 0; s_i < Nsigma; s_i++) {

                int idx = d_i * Nsigma + s_i;

                double D = Dmin + Dstep * d_i;
                double sigma = Smin + (Smax - Smin) * s_i / (Nsigma - 1);

                double X0_use[3];

                if (a_i > 0 && coex_prev[idx] == 1) {
                    X0_use[0] = X_prev[3*idx + 0];
                    X0_use[1] = X_prev[3*idx + 1];
                    X0_use[2] = X_prev[3*idx + 2];
                } else {
                    X0_use[0] = X0_default[0];
                    X0_use[1] = X0_default[1];
                    X0_use[2] = X0_default[2];
                }

                double finalX[3];
                int coex = integrate_point(D, sigma, xi, X0_use, finalX);

                sigma_arr[idx] = sigma;
                D_arr[idx]     = D;
                coex_arr[idx]  = coex;
                X_arr[3*idx + 0] = finalX[0];
                X_arr[3*idx + 1] = finalX[1];
                X_arr[3*idx + 2] = finalX[2];
            }
        }

        print_ascii_map_color_fast(ND, Nsigma, coex_arr, 0, alpha);

        // -------------------------
        // FASE 2: refinament iteratiu
        // -------------------------
        int changed = 1;
        int iter = 0;

        while (changed) {
            changed = 0;
            iter++;

            Candidate *cand = xmalloc((size_t)ND * Nsigma * sizeof(Candidate));
            int ncand = 0;

            for (int d_i = 0; d_i < ND; d_i++) {
                for (int s_i = 0; s_i < Nsigma; s_i++) {

                    int idx = d_i * Nsigma + s_i;
                    if (coex_arr[idx] == 1) continue;

                    int neigh[8][2] = {
                        {d_i-1, s_i}, {d_i+1, s_i},
                        {d_i, s_i-1}, {d_i, s_i+1},
                        {d_i-1, s_i-1}, {d_i-1, s_i+1},
                        {d_i+1, s_i-1}, {d_i+1, s_i+1}
                    };

                    int prio = 0;
                    int n_neigh = 0;
                    int n_neigh_coex = 0;

                    for (int k = 0; k < 4; k++) {
                        int di2 = neigh[k][0];
                        int si2 = neigh[k][1];
                        if (di2 < 0 || di2 >= ND ||
                            si2 < 0 || si2 >= Nsigma)
                            continue;

                        n_neigh++;
                        int idx2 = di2 * Nsigma + si2;
                        if (coex_arr[idx2] == 1) {
                            prio++;
                            n_neigh_coex++;
                        }
                    }

                    if (prio == 0) continue;
                    if (n_neigh_coex == n_neigh) continue;

                    cand[ncand].idx = idx;
                    cand[ncand].a_i = d_i;
                    cand[ncand].s_i = s_i;
                    cand[ncand].priority = prio;
                    ncand++;
                }
            }

            if (ncand == 0) {
                free(cand);
                break;
            }

            qsort(cand, ncand, sizeof(Candidate), cmp_priority);

            #pragma omp parallel for schedule(dynamic)
            for (int c = 0; c < ncand; c++) {

                int idx = cand[c].idx;
                int d_i = cand[c].a_i;
                int s_i = cand[c].s_i;

                if (coex_arr[idx] == 1) continue;

                double D = D_arr[idx];
                double sigma = sigma_arr[idx];


                int neigh[8][2] = {
                    {d_i-1, s_i}, {d_i+1, s_i},
                    {d_i, s_i-1}, {d_i, s_i+1},
                    {d_i-1, s_i-1}, {d_i-1, s_i+1},
                    {d_i+1, s_i-1}, {d_i+1, s_i+1}
                };

                int still_has_coex_neighbor = 0;
                for (int k = 0; k < 4; k++) {
                    int di2 = neigh[k][0];
                    int si2 = neigh[k][1];
                    if (di2 < 0 || di2 >= ND ||
                        si2 < 0 || si2 >= Nsigma)
                        continue;

                    int idx2 = di2 * Nsigma + si2;
                    if (coex_arr[idx2] == 1) {
                        still_has_coex_neighbor = 1;
                        break;
                    }
                }
                if (!still_has_coex_neighbor) continue;

                for (int k = 0; k < 4; k++) {
                    int di2 = neigh[k][0];
                    int si2 = neigh[k][1];
                    if (di2 < 0 || di2 >= ND ||
                        si2 < 0 || si2 >= Nsigma)
                        continue;

                    int idx2 = di2 * Nsigma + si2;
                    if (coex_arr[idx2] == 1) {

                        double X0_try[3] = {
                            X_arr[3*idx2 + 0],
                            X_arr[3*idx2 + 1],
                            X_arr[3*idx2 + 2]
                        };

                        double finalX2[3];
                        int coex2 = integrate_point(D, sigma, xi, X0_try, finalX2);

                        if (coex2 == 1) {
                            coex_arr[idx] = 1;
                            X_arr[3*idx + 0] = finalX2[0];
                            X_arr[3*idx + 1] = finalX2[1];
                            X_arr[3*idx + 2] = finalX2[2];
                            changed = 1;
                            break;
                        }
                    }
                }
            }

            free(cand);

            print_ascii_map_color_fast(ND, Nsigma, coex_arr, iter, alpha);
        }

        // -------------------------
        // Escriure fitxer refinat
        // -------------------------
        char fname[256];
        snprintf(fname, sizeof(fname), "results_alpha_%.3f_refined_200.txt", alpha);

        FILE *fout = fopen(fname, "w");
        if (fout == NULL) {
            perror(fname);
            return EXIT_FAILURE;
        }
        for (int d_i = 0; d_i < ND; d_i++) {
            for (int s_i = 0; s_i < Nsigma; s_i++) {
                int idx = d_i * Nsigma + s_i;

                fprintf(fout, "%.10f %.10f %.10f %d %.10f %.10f %.10f\n",
                        alpha,
                        D_arr[idx],
                        sigma_arr[idx],
                        coex_arr[idx],
                        X_arr[3*idx + 0],
                        X_arr[3*idx + 1],
                        X_arr[3*idx + 2]);
            }
        }
        fclose(fout);

        // guardar per seeding al alpha següent
        for (int idx = 0; idx < ND * Nsigma; idx++) {
            coex_prev[idx] = coex_arr[idx];
            X_prev[3*idx + 0] = X_arr[3*idx + 0];
            X_prev[3*idx + 1] = X_arr[3*idx + 1];
            X_prev[3*idx + 2] = X_arr[3*idx + 2];
        }

        free(sigma_arr);
        free(D_arr);
        free(coex_arr);
        free(X_arr);

        printf("Alpha=%.3f refinat en %d iteracions\n", alpha, iter);
    }

    free(X_prev);
    free(coex_prev);

    return EXIT_SUCCESS;
}
