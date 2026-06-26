#include "pds_bifurcation.h"
#include "pds_stability.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

PDSBifurcationDiagram* pds_bif_diagram_create(double p_min, double p_max,
                                                int n_steps, int param_idx) {
    if (p_min > p_max || n_steps <= 0) return NULL;
    PDSBifurcationDiagram* bd = calloc(1, sizeof(PDSBifurcationDiagram));
    if (!bd) return NULL;
    bd->param_min = p_min;
    bd->param_max = p_max;
    bd->n_param_steps = n_steps;
    bd->param_index = param_idx;
    bd->capacity = 16;
    bd->bifurcations = calloc(bd->capacity, sizeof(PDSBifurcationPoint));
    bd->param_values = calloc(n_steps, sizeof(double));
    bd->steady_states = malloc(n_steps * sizeof(double*));
    bd->n_steady_states = calloc(n_steps, sizeof(int));
    for (int i = 0; i < n_steps; i++)
        bd->param_values[i] = p_min + (p_max - p_min) * i / (double)(n_steps - 1);
    return bd;
}

void pds_bif_diagram_free(PDSBifurcationDiagram* bd) {
    if (!bd) return;
    if (bd->bifurcations) {
        for (int i = 0; i < bd->n_bifurcations; i++)
            free(bd->bifurcations[i].eigenvalues_at_bifurcation);
        free(bd->bifurcations);
    }
    if (bd->steady_states) {
        for (int i = 0; i < bd->n_param_steps; i++) free(bd->steady_states[i]);
        free(bd->steady_states);
    }
    free(bd->param_values); free(bd->n_steady_states); free(bd);
}

void pds_bif_diagram_compute(PDSBifurcationDiagram* bd,
                              void (*rhs_fn)(const double*, const double*, int, double*),
                              const double* base_params, int n_params, int n_vars,
                              const double* initial_guess) {
    if (!bd || !rhs_fn || !base_params || n_vars <= 0) return;
    double* params = malloc(n_params * sizeof(double));
    double* x = malloc(n_vars * sizeof(double));
    if (!params || !x) { free(params); free(x); return; }
    if (initial_guess) memcpy(x, initial_guess, n_vars * sizeof(double));
    else memset(x, 0, n_vars * sizeof(double));
    for (int step = 0; step < bd->n_param_steps; step++) {
        memcpy(params, base_params, n_params * sizeof(double));
        params[bd->param_index] = bd->param_values[step];
        double dt = 0.001;
        for (int r = 0; r < 5000; r++) {
            double dxdt[16];
            rhs_fn(x, params, n_params, dxdt);
            double norm = 0.0;
            for (int v = 0; v < n_vars && v < 16; v++) {
                x[v] += dt * dxdt[v];
                if (x[v] < 0.0) x[v] = 0.0;
                norm += dxdt[v] * dxdt[v];
            }
            if (sqrt(norm) < 1e-8) break;
        }
        bd->steady_states[step] = malloc(n_vars * sizeof(double));
        if (bd->steady_states[step]) {
            memcpy(bd->steady_states[step], x, n_vars * sizeof(double));
            bd->n_steady_states[step] = n_vars;
        }
    }
    free(params); free(x);
}

void pds_bif_detect(PDSBifurcationDiagram* bd) {
    if (!bd || bd->n_param_steps < 2) return;
    for (int step = 1; step < bd->n_param_steps; step++) {
        if (!bd->steady_states[step] || !bd->steady_states[step-1]) continue;
        double* x_prev = bd->steady_states[step-1];
        double* x_curr = bd->steady_states[step];
        int nv = bd->n_steady_states[step];
        double max_change = 0.0;
        for (int v = 0; v < nv; v++) {
            double change = fabs(x_curr[v] - x_prev[v]);
            if (change > max_change) max_change = change;
        }
        if (max_change > 0.5 && bd->n_bifurcations < bd->capacity) {
            PDSBifurcationPoint* bp = &bd->bifurcations[bd->n_bifurcations];
            bp->parameter_value = bd->param_values[step];
            bp->parameter_index = bd->param_index;
            bp->type = PDS_BIF_SADDLE_NODE;
            bp->is_supercritical = true;
            bp->n_eigenvalues = nv;
            bp->eigenvalues_at_bifurcation = calloc(nv * 2, sizeof(double));
            bp->amplitude_scaling = 0.5;
            bd->n_bifurcations++;
        }
    }
}

static const char* bif_name(PDSBifurcationType t) {
    switch (t) {
        case PDS_BIF_SADDLE_NODE: return "Saddle-Node";
        case PDS_BIF_PITCHFORK_SUPERCRITICAL: return "Supercritical Pitchfork";
        case PDS_BIF_PITCHFORK_SUBCRITICAL: return "Subcritical Pitchfork";
        case PDS_BIF_HOPF_SUPERCRITICAL: return "Supercritical Hopf";
        case PDS_BIF_HOPF_SUBCRITICAL: return "Subcritical Hopf";
        case PDS_BIF_TRANSCRITICAL: return "Transcritical";
        case PDS_BIF_PERIOD_DOUBLING: return "Period-Doubling";
        case PDS_BIF_HOMOCLINIC: return "Homoclinic";
        default: return "Unknown";
    }
}

void pds_bif_print(const PDSBifurcationDiagram* bd) {
    if (!bd) { puts("BifurcationDiagram: NULL"); return; }
    puts("=== Bifurcation Diagram ===");
    printf("  Parameter range: [%.4f, %.4f], %d steps\n",
           bd->param_min, bd->param_max, bd->n_param_steps);
    printf("  Detected bifurcations: %d\n", bd->n_bifurcations);
    for (int i = 0; i < bd->n_bifurcations; i++) {
        PDSBifurcationPoint* bp = &bd->bifurcations[i];
        printf("  [%d] %s at param = %.4f (%s)\n",
               i, bif_name(bp->type), bp->parameter_value,
               bp->is_supercritical ? "supercritical" : "subcritical");
    }
}

double pds_bif_normal_form_saddle_node(double param, double param_c) {
    return param - param_c;
}

double pds_bif_normal_form_hopf(double param, double param_c, double* amplitude) {
    double mu = param - param_c;
    if (amplitude) *amplitude = (mu > 0.0) ? sqrt(mu) : 0.0;
    return mu;
}