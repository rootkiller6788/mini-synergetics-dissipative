#include "sofe_core.h"
#include "sofe_bifurcation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

SOFEBifurcationPoint* sofe_bifurcation_detect(double* eigenvalues, int n,
                                                double parameter,
                                                double prev_parameter,
                                                double* prev_eigenvalues) {
    (void)prev_parameter;
    if (!eigenvalues || n <= 0) return NULL;
    SOFEBifurcationPoint* bp = (SOFEBifurcationPoint*)malloc(
        sizeof(SOFEBifurcationPoint));
    assert(bp != NULL);
    memset(bp, 0, sizeof(SOFEBifurcationPoint));
    bp->parameter_value = parameter;
    bp->n_dims = n;
    bp->critical_eigenvector = (double*)malloc(n * sizeof(double));
    bool found = false;

    for (int i = 0; i < n; i++) {
        double re = eigenvalues[2*i];
        double re_prev = prev_eigenvalues ? prev_eigenvalues[2*i] : 0.0;
        bp->critical_eigenvector[i] = 1.0;
        if (fabs(re) < 1e-8) {
            bp->type = SOFE_BIF_SADDLE_NODE;
            found = true;
        } else if (re * re_prev < 0.0 && fabs(eigenvalues[2*i+1]) > 1e-8) {
            bp->type = (re > 0.0) ? SOFE_BIF_HOPF_SUPER : SOFE_BIF_HOPF_SUB;
            bp->is_supercritical = (re > 0.0);
            found = true;
        }
    }
    if (!found) { free(bp->critical_eigenvector); free(bp); return NULL; }
    return bp;
}

void sofe_bifurcation_point_free(SOFEBifurcationPoint* bp) {
    if (!bp) return;
    free(bp->normal_form_coefficients);
    free(bp->critical_eigenvector);
    free(bp);
}

void sofe_bifurcation_point_print(SOFEBifurcationPoint* bp) {
    if (!bp) { printf("(null)\n"); return; }
    printf("Bifurcation Point:\n");
    printf("  Type: %s\n", sofe_bifurcation_name(bp->type));
    printf("  Parameter: %.6e\n", bp->parameter_value);
    printf("  Codimension: %d\n", bp->codimension);
    printf("  Supercritical: %s\n", bp->is_supercritical ? "Yes" : "No");
}

SOFEBifurcationDiagram* sofe_bifurcation_diagram_create(void) {
    SOFEBifurcationDiagram* bd = (SOFEBifurcationDiagram*)malloc(
        sizeof(SOFEBifurcationDiagram));
    assert(bd != NULL);
    memset(bd, 0, sizeof(SOFEBifurcationDiagram));
    bd->capacity = 16;
    bd->points = (SOFEBifurcationPoint*)malloc(
        bd->capacity * sizeof(SOFEBifurcationPoint));
    return bd;
}

void sofe_bifurcation_diagram_free(SOFEBifurcationDiagram* bd) {
    if (!bd) return;
    for (int i = 0; i < bd->n_points; i++)
        free(bd->points[i].critical_eigenvector);
    free(bd->points); free(bd);
}

void sofe_bifurcation_diagram_add_point(SOFEBifurcationDiagram* bd,
                                         SOFEBifurcationPoint* bp) {
    if (!bd || !bp) return;
    if (bd->n_points >= bd->capacity) {
        bd->capacity *= 2;
        bd->points = (SOFEBifurcationPoint*)realloc(bd->points,
            bd->capacity * sizeof(SOFEBifurcationPoint));
        assert(bd->points != NULL);
    }
    memcpy(&bd->points[bd->n_points], bp, sizeof(SOFEBifurcationPoint));
    bd->points[bd->n_points].critical_eigenvector =
        (double*)malloc(bp->n_dims * sizeof(double));
    memcpy(bd->points[bd->n_points].critical_eigenvector,
           bp->critical_eigenvector, bp->n_dims * sizeof(double));
    bd->n_points++;
}

void sofe_bifurcation_diagram_print(SOFEBifurcationDiagram* bd) {
    if (!bd) { printf("(null)\n"); return; }
    printf("Bifurcation Diagram: %d points\n", bd->n_points);
    for (int i = 0; i < bd->n_points; i++) {
        printf("  [%d] %s at param=%.4f\n", i,
               sofe_bifurcation_name(bd->points[i].type),
               bd->points[i].parameter_value);
    }
}
SOFEBifurcationTree* sofe_bifurcation_tree_create(int n_params) {
    SOFEBifurcationTree* bt = (SOFEBifurcationTree*)malloc(
        sizeof(SOFEBifurcationTree));
    assert(bt != NULL);
    memset(bt, 0, sizeof(SOFEBifurcationTree));
    bt->n_params = n_params;
    bt->parameter_range = (double*)malloc(2 * n_params * sizeof(double));
    bt->capacity = 8;
    bt->branches = (SOFEBranchPoint**)malloc(bt->capacity * sizeof(SOFEBranchPoint*));
    bt->branch_lengths = (int*)malloc(bt->capacity * sizeof(int));
    memset(bt->branch_lengths, 0, bt->capacity * sizeof(int));
    return bt;
}

void sofe_bifurcation_tree_free(SOFEBifurcationTree* bt) {
    if (!bt) return;
    for (int i = 0; i < bt->n_branches; i++) {
        for (int j = 0; j < bt->branch_lengths[i]; j++) {
            free(bt->branches[i][j].steady_state);
            free(bt->branches[i][j].eigenvalues_real);
            free(bt->branches[i][j].eigenvalues_imag);
        }
        free(bt->branches[i]);
    }
    free(bt->branches); free(bt->branch_lengths);
    free(bt->parameter_range); free(bt);
}

int sofe_bifurcation_tree_add_branch(SOFEBifurcationTree* bt) {
    if (!bt) return -1;
    if (bt->n_branches >= bt->capacity) {
        bt->capacity *= 2;
        bt->branches = (SOFEBranchPoint**)realloc(bt->branches,
            bt->capacity * sizeof(SOFEBranchPoint*));
        bt->branch_lengths = (int*)realloc(bt->branch_lengths,
            bt->capacity * sizeof(int));
    }
    int idx = bt->n_branches;
    bt->branches[idx] = (SOFEBranchPoint*)malloc(100 * sizeof(SOFEBranchPoint));
    bt->branch_lengths[idx] = 0;
    bt->n_branches++;
    return idx;
}

void sofe_bifurcation_tree_add_point(SOFEBifurcationTree* bt,
                                      int branch_id, double param,
                                      double* state, double* eig_real,
                                      double* eig_imag, int n_vars) {
    if (!bt || branch_id < 0 || branch_id >= bt->n_branches) return;
    int len = bt->branch_lengths[branch_id];
    SOFEBranchPoint* bp = &bt->branches[branch_id][len];
    bp->parameter_value = param;
    bp->n_vars = n_vars;
    bp->branch_id = branch_id;
    bp->steady_state = (double*)malloc(n_vars * sizeof(double));
    bp->eigenvalues_real = (double*)malloc(n_vars * sizeof(double));
    bp->eigenvalues_imag = (double*)malloc(n_vars * sizeof(double));
    if (state) memcpy(bp->steady_state, state, n_vars * sizeof(double));
    if (eig_real) memcpy(bp->eigenvalues_real, eig_real, n_vars * sizeof(double));
    if (eig_imag) memcpy(bp->eigenvalues_imag, eig_imag, n_vars * sizeof(double));
    bt->branch_lengths[branch_id]++;
}

SOFENormalForm* sofe_normal_form_compute(SOFEBifurcationType type,
                                           double* jacobian, int n,
                                           double parameter) {
    if (!jacobian || n <= 0) return NULL;
    SOFENormalForm* nf = (SOFENormalForm*)malloc(sizeof(SOFENormalForm));
    assert(nf != NULL);
    memset(nf, 0, sizeof(SOFENormalForm));
    nf->coefficients = (double*)malloc(4 * sizeof(double));
    nf->order = 3;

    switch (type) {
    case SOFE_BIF_SADDLE_NODE:
        nf->coefficients[0] = parameter;
        nf->coefficients[1] = 0.0;
        nf->coefficients[2] = jacobian[0];
        break;
    case SOFE_BIF_PITCHFORK_SUPER:
        nf->coefficients[0] = parameter;
        nf->coefficients[1] = 0.0;
        nf->coefficients[2] = -1.0;
        break;
    case SOFE_BIF_HOPF_SUPER:
        nf->coefficients[0] = parameter;
        nf->coefficients[1] = 1.0;
        nf->coefficients[2] = -1.0;
        break;
    default:
        nf->coefficients[0] = parameter;
        nf->coefficients[1] = 0.0;
        nf->coefficients[2] = 1.0;
        break;
    }
    nf->validity_radius = 1.0;
    return nf;
}

void sofe_normal_form_free(SOFENormalForm* nf) {
    if (!nf) return;
    free(nf->coefficients);
    free(nf);
}

void sofe_normal_form_print(SOFENormalForm* nf) {
    if (!nf) { printf("(null)\n"); return; }
    printf("Normal Form (order %d):\n", nf->order);
    printf("  dx/dt = %.4f + %.4f*x + %.4f*x^2 + %.4f*x^3\n",
           nf->coefficients[0], nf->coefficients[1],
           nf->coefficients[2], nf->coefficients[3]);
}

double sofe_normal_form_evaluate(SOFENormalForm* nf, double x) {
    if (!nf) return 0.0;
    double result = nf->coefficients[0];
    double xn = 1.0;
    for (int i = 1; i <= nf->order && i < 4; i++) {
        xn *= x;
        result += nf->coefficients[i] * xn;
    }
    return result;
}

double sofe_normal_form_amplitude(SOFENormalForm* nf, double parameter) {
    if (!nf) return 0.0;
    nf->coefficients[0] = parameter;
    double a = nf->coefficients[3];
    double b = nf->coefficients[2];
    if (fabs(a) < 1e-15) return 0.0;
    double disc = -b * parameter / a;
    return (disc > 0.0) ? sqrt(disc) : 0.0;
}

int sofe_bifurcation_count_stable_branches(SOFEBifurcationTree* bt,
                                            double parameter) {
    if (!bt) return 0;
    int count = 0;
    for (int i = 0; i < bt->n_branches; i++) {
        for (int j = 0; j < bt->branch_lengths[i]; j++) {
            if (fabs(bt->branches[i][j].parameter_value - parameter) < 1e-6) {
                bool stable = true;
                for (int k = 0; k < bt->branches[i][j].n_vars; k++)
                    if (bt->branches[i][j].eigenvalues_real[k] > 0.0)
                        stable = false;
                if (stable) count++;
            }
        }
    }
    return count;
}

int sofe_bifurcation_count_unstable_branches(SOFEBifurcationTree* bt,
                                              double parameter) {
    if (!bt) return 0;
    int count = 0;
    for (int i = 0; i < bt->n_branches; i++) {
        for (int j = 0; j < bt->branch_lengths[i]; j++) {
            if (fabs(bt->branches[i][j].parameter_value - parameter) < 1e-6) {
                bool stable = true;
                for (int k = 0; k < bt->branches[i][j].n_vars; k++)
                    if (bt->branches[i][j].eigenvalues_real[k] > 0.0)
                        stable = false;
                if (!stable) count++;
            }
        }
    }
    return count;
}

bool sofe_bifurcation_is_hysteresis(SOFEBifurcationDiagram* bd) {
    if (!bd || bd->n_points < 2) return false;
    for (int i = 0; i < bd->n_points; i++) {
        if (bd->points[i].type == SOFE_BIF_SADDLE_NODE) return true;
    }
    return false;
}

double sofe_pitchfork_bifurcation_amplitude(double reduced_param,
                                              double coeff_cubic,
                                              bool supercritical) {
    if (fabs(coeff_cubic) < 1e-15) return 0.0;
    if (supercritical && reduced_param > 0.0)
        return sqrt(reduced_param / fabs(coeff_cubic));
    if (!supercritical && reduced_param < 0.0)
        return sqrt(-reduced_param / fabs(coeff_cubic));
    return 0.0;
}

double sofe_hopf_bifurcation_amplitude(double reduced_param,
                                         double coeff_cubic,
                                         double frequency) {
    (void)frequency;
    if (fabs(coeff_cubic) < 1e-15) return 0.0;
    if (reduced_param > 0.0)
        return sqrt(reduced_param / fabs(coeff_cubic));
    return 0.0;
}

double sofe_saddle_node_distance(double param, double critical_param,
                                  double coeff) {
    return sqrt(fabs(param - critical_param) / (fabs(coeff) + 1e-15));
}
