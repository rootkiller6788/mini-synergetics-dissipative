#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Example 3: Turing Pattern Formation Instability Analysis
 *
 * Analyzes the Turing instability condition for a two-component
 * reaction-diffusion system. Computes the dispersion relation
 * and determines whether spatial patterns will form.
 *
 * The Turing instability is a non-equilibrium phase transition
 * from a uniform state to a spatially periodic pattern.
 *
 * Knowledge points: L6 canonical problem (Turing instability),
 *   L7 application (pattern formation, morphogenesis).
 */

/* Pattern formation declarations */
extern void *pattern_turing_analyze(double,double,double,double,
                                     double,double,double,double,int);
extern void pattern_disperse_free(void*);
extern double pattern_get_k_critical(void*);
extern int pattern_is_turing_unstable(void*);
extern double pattern_get_wavelength(void*);
extern double pattern_sh_equilibrium(double,double,double);
extern int pattern_selection(double,double,double,double*,double*,int*);

int main(void) {
    printf("============================================================\n");
    printf("  Example 3: Turing Pattern Formation Analysis\n");
    printf("============================================================\n\n");

    /* Activator-inhibitor system parameters (Gierer-Meinhardt type):
     *   du/dt = a*u - b*v + Du*d^2u/dx^2 + ...
     *   dv/dt = c*u - d*v + Dv*d^2v/dx^2 + ...
     *
     * Jacobian at homogeneous SS:
     *   J = [[a11, a12], [a21, a22]] where:
     *   a11 = activator self-enhancement (positive)
     *   a12 = inhibitor coupling (negative)
     *   a21 = activator-inhibitor cross (positive)
     *   a22 = inhibitor decay (negative)
     */

    printf("Case 1: Activator-Inhibitor with Dv > Du (Turing possible)\n");
    double a11=0.5, a12=-1.0, a21=1.0, a22=-1.5;
    double D1=1.0, D2=20.0;  /* Dv/Du = 20 */
    printf("  J = [[%.1f, %.1f], [%.1f, %.1f]], Du=%.1f, Dv=%.1f\n",
           a11, a12, a21, a22, D1, D2);

    void *dr = pattern_turing_analyze(a11, a12, a21, a22, D1, D2,
                                       0.01, 2.0, 200);
    if (dr) {
        int unstable = pattern_is_turing_unstable(dr);
        double kc = pattern_get_k_critical(dr);
        double lambda = pattern_get_wavelength(dr);
        printf("  Turing unstable: %s\n", unstable ? "YES" : "NO");
        printf("  Critical wave number k_c = %.4f\n", kc);
        printf("  Pattern wavelength = %.3f (spatial units)\n", lambda);
        pattern_disperse_free(dr);
    }

    printf("\nCase 2: Equal diffusion (no Turing instability)\n");
    D1 = 1.0; D2 = 1.0;  /* Du = Dv */
    printf("  J = [[%.1f, %.1f], [%.1f, %.1f]], Du=%.1f, Dv=%.1f\n",
           a11, a12, a21, a22, D1, D2);

    dr = pattern_turing_analyze(a11, a12, a21, a22, D1, D2,
                                 0.01, 2.0, 200);
    if (dr) {
        int unstable = pattern_is_turing_unstable(dr);
        printf("  Turing unstable: %s\n", unstable ? "YES" : "NO");
        pattern_disperse_free(dr);
    }

    /* Swift-Hohenberg pattern selection: stripes vs hexagons */
    printf("\nSwift-Hohenberg Pattern Selection:\n");
    double r = 0.2;   /* above threshold */
    double g2 = 0.1;  /* quadratic -> hexagons possible */
    double g3 = 1.0;  /* cubic saturation */

    double s_amp = 0, h_amp = 0;
    int pref = 0;
    pattern_selection(r, g2, g3, &s_amp, &h_amp, &pref);

    printf("  Control parameter r = %.2f\n", r);
    printf("  Stripe amplitude = %.4f\n", s_amp);
    printf("  Hexagon amplitude = %.4f\n", h_amp);
    printf("  Preferred pattern: %s\n",
           pref == 1 ? "HEXAGONS" : (pref == 0 ? "STRIPES" : "UNIFORM"));

    printf("\nKey takeaways:\n");
    printf("  1. Turing instability requires D_inhibitor >> D_activator\n");
    printf("  2. The critical wave number determines the pattern wavelength\n");
    printf("  3. Swift-Hohenberg is the normal form for pattern formation\n");
    printf("  4. Hexagons vs stripes: competition via quadratic nonlinearity\n");
    printf("  5. Pattern formation is a non-equilibrium phase transition\n");
    printf("     from homogeneous to spatially modulated state\n");

    return 0;
}
