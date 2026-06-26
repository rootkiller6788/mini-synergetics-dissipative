#include <stdio.h>
#include <math.h>
#include "critical_phenomena.h"

/* Example 4: Critical Exponents and Scaling Relations
 *
 * Demonstrates the computation and verification of critical
 * exponents for various universality classes. Shows all four
 * scaling relations (Rushbrooke, Widom, Fisher, Josephson).
 *
 * Knowledge points: L4 scaling relations, L2 universality,
 *   L6 canonical problem: critical exponents estimation.
 */

int main(void) {
    printf("============================================================\n");
    printf("  Example 4: Critical Exponents & Scaling Relations\n");
    printf("============================================================\n\n");

    /* Compare universality classes */
    CriticalExponents mf3 = nept_critical_exponents_mean_field(3);
    CriticalExponents is2 = nept_critical_exponents_ising_2d();
    CriticalExponents is3 = nept_critical_exponents_ising_3d();
    CriticalExponents xy3 = nept_critical_exponents_xy_3d();
    CriticalExponents heis3 = nept_critical_exponents_heisenberg_3d();

    printf("%-20s %6s %6s %6s %6s %6s %6s %6s\n",
           "Universality Class", "alpha", "beta", "gamma", "delta", "nu", "eta", "d");
    printf("----------------------------------------------------------------\n");

    CriticalExponents classes[] = {mf3, is2, is3, xy3, heis3};
    const char *names[] = {"Mean-Field (3D)", "2D Ising", "3D Ising",
                           "3D XY", "3D Heisenberg"};

    for (int i = 0; i < 5; i++) {
        printf("%-20s %6.3f %6.3f %6.3f %6.1f %6.3f %6.3f %6d\n",
               names[i],
               classes[i].alpha, classes[i].beta,
               classes[i].gamma, classes[i].delta,
               classes[i].nu, classes[i].eta_fisher,
               classes[i].dimension);
    }

    /* Verify scaling relations for each class */
    printf("\nScaling Relation Verification:\n");
    printf("%-20s %12s %12s %12s %12s\n",
           "Class", "Rushbrooke", "Widom", "Fisher", "Josephson");
    printf("----------------------------------------------------------------\n");

    double tol = 0.1;
    for (int i = 0; i < 5; i++) {
        bool results[4];
        nept_check_all_scaling_relations(&classes[i], tol, results);
        printf("%-20s %12s %12s %12s %12s\n",
               names[i],
               results[0] ? "PASS" : "FAIL",
               results[1] ? "PASS" : "FAIL",
               results[2] ? "PASS" : "FAIL",
               results[3] ? "PASS" : "FAIL");
    }

    /* Demonstrate exponent estimation from synthetic data */
    printf("\nExponent Estimation Demo:\n");
    printf("  Synthetic data: eta = B * |t|^beta with beta=0.3265 (3D Ising)\n");

    double eta_data[8], t_data[8];
    for (int i = 0; i < 8; i++) {
        t_data[i] = -(0.001 * pow(2.0, (double)i));
        eta_data[i] = 1.0 * pow(fabs(t_data[i]), 0.3265);
    }

    double amp;
    double beta_est = nept_estimate_beta(eta_data, t_data, 8, &amp);
    printf("  Estimated beta = %.4f (true = 0.3265)\n", beta_est);
    printf("  Estimated amplitude B = %.4f\n", amp);

    /* Universality class matching */
    printf("\nUniversality Class Matching:\n");
    UniversalityClass uc = nept_match_universality_class(&is3, 0.1);
    printf("  3D Ising exponents match: %s\n",
           uc == UC_ISING_3D ? "ISING_3D (correct)" : "UNKNOWN");

    printf("\nKey takeaways:\n");
    printf("  1. Critical exponents are universal (same for all systems in a class)\n");
    printf("  2. Four exact scaling relations connect the exponents\n");
    printf("  3. Hyperscaling (Josephson) holds only for d <= d_c = 4\n");
    printf("  4. Mean-field exponents are correct for d >= 4\n");
    printf("  5. Real exponents (3D) deviate significantly from mean-field\n");

    return 0;
}
