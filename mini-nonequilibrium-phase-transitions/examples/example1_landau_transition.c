#include <stdio.h>
#include <math.h>
#include "landau_theory.h"
#include "phase_transition.h"

/* Example 1: Landau Theory of Continuous Phase Transition
 *
 * Demonstrates the Landau mean-field theory for a continuous
 * (second-order) phase transition. Computes the equilibrium
 * order parameter, susceptibility, and specific heat as
 * functions of reduced temperature.
 *
 * Knowledge points: L4 Landau theory, L2 order parameter,
 *   L6 canonical problem: mean-field Ising transition.
 */

int main(void) {
    printf("============================================================\n");
    printf("  Example 1: Landau Theory of Continuous Phase Transition\n");
    printf("============================================================\n\n");

    /* Create Landau coefficients for Ising universality class */
    LandauCoefficients lc = nept_landau_ising_create(1.0, 2.0);

    printf("Landau coefficients: a=%.1f, b=%.1f (phi^4 theory)\n", lc.a, lc.b);
    printf("Transition type: %s\n\n",
           nept_landau_transition_type(lc) == NEPT_CONTINUOUS
           ? "Continuous (second-order)" : "First-order");

    printf("%-12s %-12s %-12s %-12s %-12s\n",
           "t=(T-Tc)/Tc", "eta_eq", "chi", "F(eta_eq)", "C jump");
    printf("--------------------------------------------------------------\n");

    double t_values[] = {0.1, 0.05, 0.02, 0.01, 0.0, -0.01, -0.02, -0.05, -0.1, -0.2};

    for (int i = 0; i < 10; i++) {
        double t = t_values[i];
        MeanFieldResults res;
        double eta = nept_landau_equilibrium(lc, t, &res);
        double chi = nept_landau_susceptibility(lc, t);

        printf("%-12.3f %-12.4f %-12.4f %-12.6f %-12.4f\n",
               t, eta, chi, res.free_energy, res.specific_heat_jump);
    }

    printf("\nSpecific heat jump at Tc: Delta_C = %.4f\n",
           nept_landau_specific_heat_jump(lc, 300.0));

    printf("\nGinzburg criterion analysis:\n");
    GinzburgCriterion gc = nept_ginzburg_criterion(1.0, 1.0, 75.0, 300.0, 3);
    printf("  t_G (Ginzburg reduced temp) = %.6f\n", gc.ginzburg_reduced_temp);
    printf("  Mean-field valid? %s\n", gc.mean_field_valid ? "YES" : "NO");
    printf("  Upper critical dimension d_c = %d\n", gc.upper_critical_dimension);

    printf("\nKey takeaways:\n");
    printf("  1. Order parameter eta ~ |t|^beta with beta=0.5 (mean-field)\n");
    printf("  2. Susceptibility chi ~ |t|^-gamma with gamma=1.0\n");
    printf("  3. Specific heat has finite jump Delta_C at Tc\n");
    printf("  4. Ginzburg criterion: mean-field valid for d >= 4\n");

    return 0;
}
