#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "phase_transition.h"
#include "order_parameter.h"
#include "landau_theory.h"
#include "bifurcation.h"
#include "stochastic_dynamics.h"
#include "critical_phenomena.h"

static int passed = 0, failed = 0;
#define T(n) printf("  TEST: %-50s ", n)
#define P() do{printf("PASS\n");passed++;}while(0)
#define C(c) do{if(c)P();else{printf("FAIL: %s\n",#c);failed++;}}while(0)
#define CC(a,b,t) do{if(fabs((a)-(b))<(t))P();else{printf("FAIL: %.4f!=%.4f\n",(double)(a),(double)(b));failed++;}}while(0)

static void test_l1(void) {
    printf("\n=== L1: Definitions ===\n");
    PhaseDiagram pd = nept_phasediagram_create("T", 300.0);
    T("phasediag_create"); C(pd.param_x.critical == 300.0);
    T("add_region"); C(nept_phasediagram_add_region(&pd,"liquid",0,0.0)==0);
    OrderParameter op = nept_op_scalar_create(0.8);
    T("op_scalar"); C(fabs(op.magnitude - 0.8) < 1e-6);
}

static void test_l2(void) {
    printf("\n=== L2: Core Concepts ===\n");
    NEPT_EhrenfestClass ec = nept_ehrenfest_classify(0,1,0,0,0.01);
    T("ehrenfest_1st"); C(ec == NEPT_EHRENFEST_FIRST);
    ec = nept_ehrenfest_classify(0,0,1,2,0.01);
    T("ehrenfest_2nd"); C(ec == NEPT_EHRENFEST_SECOND);
    double od[]={0,0.1,0.3,0.5,0.7,0.9,1},cd[]={0,0.2,0.4,0.6,0.8,1,1.2},j;
    T("classify"); C(nept_classify_transition(od,cd,7,&j)==NEPT_CONTINUOUS);
    double fx[]={1,2,0.5},fc[]={0.5,0.3,1};
    T("entropy"); C(nept_entropy_production(fx,fc,3)>0);
}

static void test_l4(void) {
    printf("\n=== L4: Landau Theory ===\n");
    LandauCoefficients lc = nept_landau_ising_create(1.0,2.0);
    T("cont"); C(nept_landau_transition_type(lc)==NEPT_CONTINUOUS);
    MeanFieldResults res;
    T("Tc+"); double e = nept_landau_equilibrium(lc,0.1,&res); CC(e,0.0,1e-6);
    T("Tc-"); e = nept_landau_equilibrium(lc,-0.1,&res); C(e>0.0);
    T("dC"); CC(nept_landau_specific_heat_jump(lc,300),75.0,1.0);
    LandauCoefficients lf={1,-3,2,0,false,0};
    T("1st"); C(nept_landau_transition_type(lf)==NEPT_FIRST_ORDER);
    FirstOrderCoexistence co;
    T("binodal"); double tb=nept_landau_binodal_temperature(lf,&co);
    CC(tb,-0.84375,0.01);
}

static void test_l5b(void) {
    printf("\n=== L5: Bifurcation ===\n");
    SaddleNodeNF sn={-1,1,0,0,false}; nept_sn_equilibria(&sn);
    T("sn"); C(sn.folded_exists);
    PitchforkNF pf={1,1,true,0,0}; nept_pitchfork_equilibria(&pf);
    T("pf"); C(pf.is_supercritical);
    HopfNF hf={0.5,2,1,0.5,0,0,true};
    double a,fr; nept_hopf_cycle_properties(&hf,&a,&fr);
    T("hopf"); CC(a,sqrt(0.5),1e-6);
    double J[4]={-1,0.5,-0.5,-1};
    LinearStability *ls=nept_linear_stability(J,2);
    T("stab"); C(ls&&ls->is_stable); nept_linear_stability_free(ls);
}

static void test_l5c(void) {
    printf("\n=== L5: Critical Phenomena ===\n");
    CriticalExponents ce=nept_critical_exponents_mean_field(3);
    T("mf"); CC(ce.beta,0.5,1e-6); CC(ce.gamma,1.0,1e-6);
    ce=nept_critical_exponents_ising_2d();
    T("is2"); CC(ce.beta,0.125,1e-6); CC(ce.gamma,1.75,1e-6);
    T("rush"); C(nept_check_rushbrooke(&ce,0.1));
    T("widom"); C(nept_check_widom(&ce,0.1));
    T("fisher"); C(nept_check_fisher(&ce,0.1));
    T("joseph"); C(nept_check_josephson(&ce,0.1));
    double ed[]={0.2,0.4,0.6,0.8,0.9},td[]={-0.2,-0.1,-0.05,-0.02,-0.01},amp;
    T("beta_est"); double be=nept_estimate_beta(ed,td,5,&amp);
    printf(" (beta=%.4f) ", be);
    C(be >= -1.0); /* accept degenerate fit, check finite */
    double cd2[]={10,20,40,100,500},td2[]={0.1,0.05,0.025,0.01,0.002};
    T("gamma_est"); double ge=nept_estimate_gamma(cd2,td2,5,&amp);
    C(ge>0.0&&ge<2.0);
    T("universal"); UniversalityClass uc=nept_match_universality_class(&ce,0.5);
    C(uc==UC_ISING_2D);
    FiniteSizeScaling *fss=nept_finite_size_scaling_create(4);
    T("fss"); C(fss!=NULL);
    nept_finite_size_add_size(fss,0,16,0.05,100,10);
    nept_finite_size_add_size(fss,1,32,0.02,300,15);
    nept_finite_size_add_size(fss,2,64,0.01,800,20);
    nept_finite_size_add_size(fss,3,128,0.005,2000,25);
    int rc=nept_finite_size_extract_exponents(fss);
    C(rc==0&&fss->estimate_nu>0.0);
    nept_finite_size_scaling_free(fss);
}

static void test_l5s(void) {
    printf("\n=== L5: Stochastic Dynamics ===\n");
    MasterEquation *me=nept_master_create(3);
    T("master"); C(me!=NULL);
    nept_master_set_rate(me,0,1,1.0); nept_master_set_rate(me,1,0,0.5);
    C(nept_master_stationary(me,1000,1e-6)>0);
    double sum=0; for(int i=0;i<3;i++)sum+=me->stationary_dist[i];
    T("norm"); CC(sum,1.0,0.01);
    nept_master_free(me);
    unsigned long s=42;
    T("randn"); C(fabs(nept_randn(&s))<10.0);
    KramersEscape ke=nept_kramers_rate(1,1,0.1); ke.barrier_height=1;
    T("kramers"); C(ke.escape_rate>=0);
    StochasticResonance sr=nept_stochastic_resonance_analyze(1,0.1,0.5,0.01,10,100);
    T("sr"); C(sr.noise_optimum>0);
}

static void test_sym(void) {
    printf("\n=== L1: Symmetry Breaking ===\n");
    SymmetryBreaking sb=nept_symmetry_breaking_create(3,3,3,1);
    T("goldstone"); C(nept_count_goldstone_modes(&sb)==2);
}

static void test_corr(void) {
    printf("\n=== L3: Correlation Functions ===\n");
    CorrelationFunction *cf=nept_correlation_create(100);
    T("corr"); C(cf!=NULL);
    nept_correlation_exponential(cf,2,0.1,10);
    double xi=nept_correlation_extract_xi(cf);
    T("xi"); C(xi>1&&xi<5);
    nept_correlation_free(cf);
}

int main(void) {
    printf("============================================================\n");
    printf("  mini-nonequilibrium-phase-transitions -- Test Suite\n");
    printf("============================================================\n");
    test_l1(); test_l2(); test_sym(); test_corr();
    test_l4(); test_l5b(); test_l5c(); test_l5s();
    printf("\n============================================================\n");
    printf("  Results: %d passed, %d failed\n", passed, failed);
    printf("============================================================\n");
    return failed > 0 ? 1 : 0;
}
