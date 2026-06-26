#include "pds_core.h"
#include "pds_thermo.h"
#include "pds_stability.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* L8: Advanced Topics in Dissipative Structure Theory */

void stochastic_resonance_simulate(double x0, double A, double omega,
                                    double D, double dt, int n_steps,
                                    double* trajectory, double* signal) {
    double x = x0;
    for (int i = 0; i < n_steps; i++) {
        double t = i * dt;
        double force = x - x*x*x;
        force += A * cos(omega * t);
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        if (u1 < 1e-10) u1 = 1e-10;
        double noise = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        force += sqrt(2.0 * D / dt) * noise;
        x += dt * force;
        if (trajectory) trajectory[i] = x;
        if (signal) signal[i] = A * cos(omega * t);
    }
}

void cgle_rhs_1d(const double* A_real, const double* A_imag, int n,
                  double mu, double alpha, double beta, double dx,
                  double* dAr_dt, double* dAi_dt) {
    for (int i = 0; i < n; i++) {
        double Ar = A_real[i], Ai = A_imag[i];
        double A2 = Ar*Ar + Ai*Ai;
        dAr_dt[i] = mu*Ar - A2*(Ar - alpha*Ai);
        dAi_dt[i] = mu*Ai - A2*(Ai + alpha*Ar);
        double lap_Ar = 0.0, lap_Ai = 0.0;
        if (i > 0 && i < n-1) {
            lap_Ar = (A_real[i+1] - 2.0*Ar + A_real[i-1]) / (dx*dx);
            lap_Ai = (A_imag[i+1] - 2.0*Ai + A_imag[i-1]) / (dx*dx);
        }
        dAr_dt[i] += lap_Ar - beta*lap_Ai;
        dAi_dt[i] += lap_Ai + beta*lap_Ar;
    }
}

double lyapunov_function_brusselator(double x, double y, double xs, double ys) {
    double dx = x - xs, dy = y - ys;
    return 0.5 * (dx*dx + dy*dy);
}

void game_of_life_step(const bool* grid, bool* next, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int neighbors = 0;
            for (int di = -1; di <= 1; di++) {
                for (int dj = -1; dj <= 1; dj++) {
                    if (di == 0 && dj == 0) continue;
                    int ni = (i + di + rows) % rows;
                    int nj = (j + dj + cols) % cols;
                    if (grid[ni * cols + nj]) neighbors++;
                }
            }
            bool alive = grid[i * cols + j];
            next[i * cols + j] = (alive && (neighbors == 2 || neighbors == 3))
                                 || (!alive && neighbors == 3);
        }
    }
}

typedef struct { double lower; double central; double upper; } FuzzyNumber;

FuzzyNumber fuzzy_multiply(FuzzyNumber a, FuzzyNumber b) {
    FuzzyNumber r;
    double p1 = a.lower * b.lower, p2 = a.lower * b.upper;
    double p3 = a.upper * b.lower, p4 = a.upper * b.upper;
    r.lower = fmin(fmin(p1, p2), fmin(p3, p4));
    r.upper = fmax(fmax(p1, p2), fmax(p3, p4));
    r.central = a.central * b.central;
    return r;
}

double fuzzy_entropy_production(FuzzyNumber flux, FuzzyNumber force) {
    FuzzyNumber prod = fuzzy_multiply(flux, force);
    if (prod.lower < 0.0) prod.lower = 0.0;
    return prod.central;
}
