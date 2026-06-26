#include "pds_brusselator.h"
#include "pds_stability.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    puts("Bifurcation Analysis of the Brusselator");
    puts("========================================");
    puts("Hopf bifurcation condition: b_c = 1 + a^2");
    double a = 1.0;
    printf("Fixed a = %.1f", a); putchar(0x0a);
    printf("  b_c = 1 + %.1f^2 = %.2f", a, 1.0 + a*a); putchar(0x0a); putchar(0x0a);
    puts("Parameter scan (varying b):");
    puts("  b      | Stability | Oscillatory | Notes");
    puts("  -------|-----------|-------------|-------");
    for (int i = 0; i < 6; i++) {
        double b = 0.5 + i * 0.5;
        PDSBrusselatorParams params = {a, b, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        double bc = pds_brusselator_bifurcation_point(&params);
        PDSBrusselatorState ss; ss.x = a; ss.y = b / a;
        double J[2][2];
        pds_brusselator_jacobian(&ss, &params, J);
        double trace = J[0][0] + J[1][1];
        double det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        double disc = trace*trace - 4.0*det;
        const char* st = (trace < 0) ? "STABLE" : "UNSTABLE";
        const char* os = (disc < 0) ? "YES" : "no";
        const char* nt = (b < bc) ? "thermo branch" : ((b > bc) ? "** LIMIT CYCLE **" : "** BIFURCATION **");
        printf("  %.1f    | %-9s | %-11s | %s", b, st, os, nt); putchar(0x0a);
    }
    return 0;
}
