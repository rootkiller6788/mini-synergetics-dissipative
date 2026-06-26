#ifndef PDS_BIFURCATION_H
#define PDS_BIFURCATION_H

#include "pds_core.h"

typedef enum {
    PDS_BIF_NONE = 0,
    PDS_BIF_SADDLE_NODE = 1,
    PDS_BIF_PITCHFORK_SUPERCRITICAL = 2,
    PDS_BIF_PITCHFORK_SUBCRITICAL = 3,
    PDS_BIF_HOPF_SUPERCRITICAL = 4,
    PDS_BIF_HOPF_SUBCRITICAL = 5,
    PDS_BIF_TRANSCRITICAL = 6,
    PDS_BIF_PERIOD_DOUBLING = 7,
    PDS_BIF_HOMOCLINIC = 8
} PDSBifurcationType;

typedef struct {
    PDSBifurcationType type;
    double parameter_value;
    const char* parameter_name;
    int parameter_index;
    bool is_supercritical;
    double* eigenvalues_at_bifurcation;
    int n_eigenvalues;
    double amplitude_scaling;
} PDSBifurcationPoint;

typedef struct {
    PDSBifurcationPoint* bifurcations;
    int n_bifurcations;
    int capacity;
    double param_min, param_max;
    int param_index;
    double* param_values;
    double** steady_states;
    int* n_steady_states;
    int n_param_steps;
} PDSBifurcationDiagram;

PDSBifurcationDiagram* pds_bif_diagram_create(double p_min, double p_max,
                                                int n_steps, int param_idx);
void pds_bif_diagram_free(PDSBifurcationDiagram* bd);
void pds_bif_diagram_compute(PDSBifurcationDiagram* bd,
                              void (*rhs_fn)(const double*, const double*, int, double*),
                              const double* base_params, int n_params, int n_vars,
                              const double* initial_guess);
void pds_bif_detect(PDSBifurcationDiagram* bd);
void pds_bif_print(const PDSBifurcationDiagram* bd);
double pds_bif_normal_form_saddle_node(double param, double param_c);
double pds_bif_normal_form_hopf(double param, double param_c, double* amplitude);

#endif