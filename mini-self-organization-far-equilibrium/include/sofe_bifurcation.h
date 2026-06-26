#ifndef SOFE_BIFURCATION_H
#define SOFE_BIFURCATION_H

#include "sofe_core.h"
#include "sofe_dynamics.h"

typedef struct {
    SOFEBifurcationType type;
    double parameter_value;
    int codimension;
    double* normal_form_coefficients;
    int n_coeffs;
    double* critical_eigenvector;
    int n_dims;
    bool is_supercritical;
    double unfolding_param1;
    double unfolding_param2;
} SOFEBifurcationPoint;

typedef struct {
    SOFEBifurcationPoint* points;
    int n_points;
    int capacity;
    double parameter_start;
    double parameter_end;
} SOFEBifurcationDiagram;

typedef struct {
    double parameter_value;
    double* steady_state;
    double* eigenvalues_real;
    double* eigenvalues_imag;
    int n_vars;
    int branch_id;
} SOFEBranchPoint;

typedef struct {
    SOFEBranchPoint** branches;
    int* branch_lengths;
    int n_branches;
    int capacity;
    double* parameter_range;
    int n_params;
} SOFEBifurcationTree;

typedef struct {
    double* coefficients;
    int order;
    double truncation_error;
    double validity_radius;
} SOFENormalForm;

SOFEBifurcationPoint* sofe_bifurcation_detect(double* eigenvalues, int n,
                                                double parameter,
                                                double prev_parameter,
                                                double* prev_eigenvalues);
void sofe_bifurcation_point_free(SOFEBifurcationPoint* bp);
void sofe_bifurcation_point_print(SOFEBifurcationPoint* bp);

SOFEBifurcationDiagram* sofe_bifurcation_diagram_create(void);
void sofe_bifurcation_diagram_free(SOFEBifurcationDiagram* bd);
void sofe_bifurcation_diagram_add_point(SOFEBifurcationDiagram* bd,
                                         SOFEBifurcationPoint* bp);
void sofe_bifurcation_diagram_print(SOFEBifurcationDiagram* bd);

SOFEBifurcationTree* sofe_bifurcation_tree_create(int n_params);
void sofe_bifurcation_tree_free(SOFEBifurcationTree* bt);
int sofe_bifurcation_tree_add_branch(SOFEBifurcationTree* bt);
void sofe_bifurcation_tree_add_point(SOFEBifurcationTree* bt,
                                      int branch_id, double param,
                                      double* state, double* eig_real,
                                      double* eig_imag, int n_vars);

SOFENormalForm* sofe_normal_form_compute(SOFEBifurcationType type,
                                           double* jacobian, int n,
                                           double parameter);
void sofe_normal_form_free(SOFENormalForm* nf);
void sofe_normal_form_print(SOFENormalForm* nf);
double sofe_normal_form_evaluate(SOFENormalForm* nf, double x);
double sofe_normal_form_amplitude(SOFENormalForm* nf, double parameter);

int sofe_bifurcation_count_stable_branches(SOFEBifurcationTree* bt,
                                            double parameter);
int sofe_bifurcation_count_unstable_branches(SOFEBifurcationTree* bt,
                                              double parameter);
bool sofe_bifurcation_is_hysteresis(SOFEBifurcationDiagram* bd);

double sofe_pitchfork_bifurcation_amplitude(double reduced_param,
                                              double coeff_cubic,
                                              bool supercritical);
double sofe_hopf_bifurcation_amplitude(double reduced_param,
                                         double coeff_cubic,
                                         double frequency);
double sofe_saddle_node_distance(double param, double critical_param,
                                  double coeff);

#endif
