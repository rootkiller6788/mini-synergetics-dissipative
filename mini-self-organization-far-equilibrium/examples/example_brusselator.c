#include "sofe_core.h"
#include "sofe_dynamics.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Brusselator: Self-Organization via Chemical Oscillations ===\n");
    printf("Model: A -> X, B+X -> Y+D, 2X+Y -> 3X, X -> E\n");
    printf("Parameters: A=1.0, B=3.0 (B > 1+A^2 gives limit cycle)\n\n");

    double y[] = {1.0, 1.0};
    double p[] = {1.0, 3.0};
    double t = 0.0, dt = 0.01;
    int n_steps = 2000;

    printf("t\t\tX\t\tY\n");
    for (int step = 0; step < n_steps; step++) {
        if (step % 100 == 0) {
            printf("%.2f\t\t%.4f\t\t%.4f\n", t, y[0], y[1]);
        }
        sofe_integrate_step(y, 2, t, dt, sofe_brusselator_rates, p, 2, SOFE_INT_RK4);
        t += dt;
    }

    printf("\nFinal state: X=%.4f, Y=%.4f\n", y[0], y[1]);
    printf("The Brusselator demonstrates temporal self-organization\n");
    printf("(limit cycle oscillations) far from equilibrium.\n");
    return 0;
}
