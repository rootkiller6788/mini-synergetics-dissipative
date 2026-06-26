#include "pds_core.h"
#include "pds_thermo.h"
#include "pds_brusselator.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

void oregonator_rhs(const double* state, const double* params, int n, double* deriv) {
    (void)n;
    double X = state[0], Y = state[1], Z = state[2];
    double A = params[0], B = params[1];
    double k1 = 1.28, k2 = 2.4e6, k3 = 33.6, k4 = 3.0e3, kc = 1.0;
    deriv[0] = k1*A*Y - k2*X*Y + k3*A*X - 2.0*k4*X*X;
    deriv[1] = -k1*A*Y - k2*X*Y + 0.5*kc*B*Z;
    deriv[2] = 2.0*k3*A*X - kc*B*Z;
}

void glycolytic_oscillator_rhs(const double* state, const double* params,
                                int n, double* deriv) {
    (void)n;
    double x = state[0], y = state[1];
    double nu = params[0];
    deriv[0] = -x + 0.5*y + x*x*y;
    deriv[1] = nu - 0.5*y - x*x*y;
}

double climate_energy_balance(double T, const double* params, int n) {
    (void)n;
    double S = params[0], alpha = params[1], epsilon = params[2], C = params[3];
    double sigma = 5.67e-8;
    return ((1.0 - alpha) * S / 4.0 - epsilon * sigma * T*T*T*T) / C;
}

void pds_applications_summary(void) {
    printf("=== PDS Applications (L7) ==="); putchar(0x0a);
    printf("1. Chemical Oscillations (BZ Reaction): Oregonator model"); putchar(0x0a);
    printf("   - NASA ISS microgravity experiments (2020)"); putchar(0x0a);
    printf("2. Climate Systems: Earth as dissipative structure"); putchar(0x0a);
    printf("   - Solar-driven, IR-cooled, multiple steady states"); putchar(0x0a);
    printf("3. Biochemical Oscillations: Glycolytic oscillator"); putchar(0x0a);
    printf("   - Yeast cell metabolism, Selkov model (1968)"); putchar(0x0a);
    printf("4. Industrial Reactor: Toyota lean production"); putchar(0x0a);
    printf("5. Nuclear Safety: Fukushima reactor entropy analysis"); putchar(0x0a);
    printf("6. Smart Grid: Power grid frequency regulation"); putchar(0x0a);
    printf("7. Ecological Systems: Predator-prey dynamics"); putchar(0x0a);
}

void pds_application_integrate(
    void (*rhs)(const double*, const double*, int, double*),
    double* state, int n_vars, const double* params, int n_params,
    double t_start, double t_end, double dt,
    void (*callback)(double, const double*, int)) {
    double t = t_start;
    double x[16];
    for (int i = 0; i < n_vars && i < 16; i++) x[i] = state[i];
    int n_steps = (int)((t_end - t_start) / dt);
    for (int step = 0; step < n_steps; step++) {
        double dxdt[16];
        rhs(x, params, n_params, dxdt);
        for (int i = 0; i < n_vars && i < 16; i++) {
            x[i] += dt * dxdt[i];
            if (x[i] < 0.0) x[i] = 0.0;
        }
        t += dt;
        if (callback && step % 100 == 0) callback(t, x, n_vars);
    }
    for (int i = 0; i < n_vars && i < 16; i++) state[i] = x[i];
}
