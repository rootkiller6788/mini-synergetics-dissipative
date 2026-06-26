/**
 * example_stability_analysis.c — Turing Instability Analysis Tool
 *
 * Systematically analyzes the Turing instability conditions across
 * multiple models, performing parameter sweeps to map out Turing spaces.
 * This demonstrates the linear stability analysis pipeline:
 *   HSS → Jacobian → Turing conditions → dispersion → Turing space
 *
 * Knowledge Coverage: L6 — Canonical problem: parameter space mapping
 *                     L4 — Verification of Turing instability theorem
 */

#include "turing_common.h"
#include "turing_models.h"
#include "turing_analysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static void analyze_model(ModelType type, const char *label,
                          ModelParams *params, double Du, double Dv) {
    printf("\n========================================\n");
    printf("  %s\n", label);
    printf("========================================\n");

    model_params_print(type, params);
    printf("  Du=%.4g, Dv=%.4g (Dv/Du=%.2f)\n", Du, Dv, Dv/Du);

    /* Find HSS */
    HomogeneousSteadyState hss;
    int hss_ret = model_find_hss(type, params, &hss);
    if (hss_ret != 0 || !hss.converged) {
        printf("  ⚠ No homogeneous steady state found.\n");
        return;
    }
    printf("  HSS: u*=%.6g, v*=%.6g (residual=%.2e)\n",
           hss.u_star, hss.v_star, hss.residual);

    /* Jacobian */
    double J[2][2];
    model_jacobian(type, hss.u_star, hss.v_star, params, J);
    printf("  Jacobian: [[%.4g, %.4g], [%.4g, %.4g]]\n",
           J[0][0], J[0][1], J[1][0], J[1][1]);

    /* Turing conditions */
    TuringConditions tc;
    turing_conditions_compute(J, Du, Dv, &tc);
    turing_conditions_print(&tc);

    /* Dispersion */
    if (tc.turing_unstable) {
        DispersionRelation *dr = turing_dispersion_compute(
            J, Du, Dv, 200, 0.01, 5.0);
        if (dr) {
            turing_dispersion_print(dr);
            dispersion_free(dr);
        }
    }

    /* Hopf check */
    int hopf = is_hopf_unstable(J);
    printf("  Hopf-unstable: %s\n", hopf ? "YES (temporal oscillations)" : "NO");

    /* Turing-Hopf interaction */
    if (tc.turing_unstable && hopf) {
        printf("  ➜ Turing-Hopf codimension-2 point detected!\n");
    }

    /* Critical wavenumber */
    double kc = turing_critical_wavenumber(J, Du, Dv);
    if (kc > 0.0) {
        printf("  Critical wavenumber k_c = %.4g (wavelength λ = %.4g)\n",
               kc, 2.0 * M_PI / kc);
    }
}

static void parameter_sweep_demo(ModelType type,
                                  int param_idx,
                                  const char *param_name,
                                  double p_min, double p_max, int n_steps,
                                  ModelParams *base, double Du, double Dv) {
    printf("\n--- Parameter Sweep: %s (%s) ---\n",
           model_name(type), param_name);
    printf("%-12s %-10s %-10s %-12s %s\n",
           param_name, "tr(J)", "det(J)", "cross", "Turing?");

    for (int i = 0; i < n_steps; i++) {
        double p_val = p_min + (p_max - p_min) * (double)i / (double)(n_steps-1);
        ModelParams p = *base;

        /* Set the sweep parameter */
        switch (param_idx) {
            case 0: p.F = p_val; break;
            case 1: p.k = p_val; break;
            case 2: p.A = p_val; break;
            case 3: p.B = p_val; break;
            case 4: p.a = p_val; break;
            case 5: p.b = p_val; break;
            default: break;
        }

        HomogeneousSteadyState hss;
        if (model_find_hss(type, &p, &hss) != 0 || !hss.converged) {
            printf("%-12.4g %-10s %-10s %-12s %s\n",
                   p_val, "---", "---", "---", "NO HSS");
            continue;
        }

        double J[2][2];
        model_jacobian(type, hss.u_star, hss.v_star, &p, J);

        double tr = J[0][0] + J[1][1];
        double det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        double cross = Dv*J[0][0] + Du*J[1][1];
        double disc = cross*cross - 4.0*Du*Dv*det;

        int turing = (tr < 0 && det > 0 && cross > 0 && disc > 0);
        printf("%-12.4g %-10.4g %-10.4g %-12.4g %s\n",
               p_val, tr, det, cross, turing ? "YES ✓" : "no");
    }
}

int main(void) {
    printf("=== Turing Instability Analysis Tool ===\n");

    /* Model 1: Brusselator */
    {
        ModelParams p;
        model_params_default(MODEL_BRUSSELATOR, &p);
        p.A = 1.0; p.B = 3.0;
        analyze_model(MODEL_BRUSSELATOR, "Brusselator", &p, 1.0, 20.0);
    }

    /* Model 2: Schnakenberg */
    {
        ModelParams p;
        model_params_default(MODEL_SCHNAKENBERG, &p);
        p.a = 0.1; p.b = 0.9;
        analyze_model(MODEL_SCHNAKENBERG, "Schnakenberg", &p, 1.0, 40.0);
    }

    /* Model 3: Gray-Scott */
    {
        ModelParams p;
        model_params_default(MODEL_GRAY_SCOTT, &p);
        p.F = 0.04; p.k = 0.06;
        analyze_model(MODEL_GRAY_SCOTT, "Gray-Scott", &p, 2e-5, 1e-5);
    }

    /* Parameter sweep: Schnakenberg, vary b */
    printf("\n\n=== Schnakenberg: Parameter Sweep (vary b) ===\n");
    {
        ModelParams base;
        model_params_default(MODEL_SCHNAKENBERG, &base);
        base.a = 0.1;
        parameter_sweep_demo(MODEL_SCHNAKENBERG, 5, "b",
                             0.01, 2.0, 21, &base, 1.0, 40.0);
    }

    /* Parameter sweep: Brusselator, vary B */
    printf("\n\n=== Brusselator: Parameter Sweep (vary B) ===\n");
    {
        ModelParams base;
        model_params_default(MODEL_BRUSSELATOR, &base);
        base.A = 1.0;
        parameter_sweep_demo(MODEL_BRUSSELATOR, 3, "B",
                             0.5, 5.0, 21, &base, 1.0, 20.0);
    }

    /* Turing Space for Schnakenberg */
    printf("\n\n=== Schnakenberg Turing Space (a vs b) ===\n");
    {
        ModelParams base;
        model_params_default(MODEL_SCHNAKENBERG, &base);
        TuringSpace *ts = turing_space_alloc(
            MODEL_SCHNAKENBERG, &base, 1.0, 40.0,
            4, 0.05, 0.30, 40,   /* param a */
            5, 0.05, 1.50, 40);  /* param b */
        if (ts) {
            SolverConfig sc = {0};
            turing_space_compute(ts, &sc);
            turing_space_print_summary(ts);
            turing_space_free(ts);
        }
    }

    printf("\nAnalysis complete.\n");
    return 0;
}
