/* ============================================================================
 * sofe_core.c -- Core Implementation of Self-Organization Far-Equilibrium
 *
 * Implements fundamental data structures and operations for:
 *   - Reaction network construction and management
 *   - Spatial grid and field operations
 *   - Order parameter management
 *   - System-level lifecycle and state tracking
 *
 * References:
 *   Nicolis & Prigogine, Self-Organization in Nonequilibrium Systems (1977)
 *   Haken, Synergetics: An Introduction (1977)
 *   Feinberg, Foundations of Chemical Reaction Network Theory (2019)
 *   Kondepudi & Prigogine, Modern Thermodynamics (2015), Ch. 15-19
 * ============================================================================ */

#include "sofe_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ============================================================================
 * SOFE Complex Number Operations
 * ============================================================================ */

SOFEComplex sofe_complex_add(SOFEComplex a, SOFEComplex b) {
    SOFEComplex r;
    r.real = a.real + b.real;
    r.imag = a.imag + b.imag;
    return r;
}

SOFEComplex sofe_complex_mul(SOFEComplex a, SOFEComplex b) {
    SOFEComplex r;
    r.real = a.real * b.real - a.imag * b.imag;
    r.imag = a.real * b.imag + a.imag * b.real;
    return r;
}

double sofe_complex_abs(SOFEComplex z) {
    return sqrt(z.real * z.real + z.imag * z.imag);
}

/* ============================================================================
 * Vector and Matrix Operations
 * ============================================================================ */

double sofe_vector_norm(double* v, int n) {
    if (!v || n <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

double sofe_vector_dot(double* a, double* b, int n) {
    if (!a || !b || n <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
}

void sofe_matrix_vector_mul(double* out, SOFEMatrix* A, double* x) {
    if (!out || !A || !x) return;
    for (int i = 0; i < A->rows; i++) {
        out[i] = 0.0;
        for (int j = 0; j < A->cols; j++) {
            out[i] += A->data[i * A->cols + j] * x[j];
        }
    }
}

void sofe_matrix_print(SOFEMatrix* m) {
    if (!m) { printf("(null matrix)\n"); return; }
    printf("Matrix %dx%d:\n", m->rows, m->cols);
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            printf("% 8.4f ", m->data[i * m->cols + j]);
        }
        printf("\n");
    }
}

void sofe_vector_print(SOFEVector* v) {
    if (!v) { printf("(null vector)\n"); return; }
    printf("Vector[%d]: [", v->dim);
    for (int i = 0; i < v->dim; i++) {
        printf("%.4f%s", v->components[i], i < v->dim - 1 ? ", " : "");
    }
    printf("]\n");
}

/* ============================================================================
 * System Lifecycle
 * ============================================================================ */

SOFESystem* sofe_system_create(const char* name) {
    SOFESystem* sys = (SOFESystem*)malloc(sizeof(SOFESystem));
    assert(sys != NULL);
    memset(sys, 0, sizeof(SOFESystem));
    sys->name = (char*)malloc(strlen(name) + 1);
    assert(sys->name != NULL);
    strcpy(sys->name, name);
    sys->regime = SOFE_REGIME_NEAR_EQUILIBRIUM;
    sys->stability = SOFE_STABLE_NODE;
    sys->temperature = 298.15;
    sys->has_macroscopic_order = false;
    sys->order_param_capacity = 4;
    sys->order_params = (SOFEOrderParameter*)malloc(
        sys->order_param_capacity * sizeof(SOFEOrderParameter));
    assert(sys->order_params != NULL);
    memset(sys->order_params, 0,
           sys->order_param_capacity * sizeof(SOFEOrderParameter));
    sys->n_order_params = 0;
    sys->n_control_params = 0;
    sys->history_capacity = 1000;
    sys->history_len = 0;
    sys->history_entropy_prod = (double*)malloc(
        sys->history_capacity * sizeof(double));
    sys->history_order_param = (double*)malloc(
        sys->history_capacity * sizeof(double));
    assert(sys->history_entropy_prod && sys->history_order_param);
    sys->lyapunov_exponent_max = 0.0;
    return sys;
}

void sofe_system_free(SOFESystem* sys) {
    if (!sys) return;
    free(sys->name);
    if (sys->network) sofe_network_free(sys->network);
    if (sys->concentration_field) sofe_field_free(sys->concentration_field);
    for (int i = 0; i < sys->n_order_params; i++) {
        free(sys->order_params[i].name);
        free(sys->order_params[i].amplitude);
        free(sys->order_params[i].phase);
    }
    free(sys->order_params);
    free(sys->eigenvalues);
    if (sys->control_params) {
        free(sys->control_params);
        for (int i = 0; i < sys->n_control_params; i++) {
            free(sys->control_param_names[i]);
        }
        free(sys->control_param_names);
    }
    free(sys->history_entropy_prod);
    free(sys->history_order_param);
    free(sys);
}

void sofe_system_set_regime(SOFESystem* sys, SOFERegime regime) {
    if (!sys) return;
    sys->regime = regime;
}

void sofe_system_set_temperature(SOFESystem* sys, double T) {
    if (!sys || T <= 0.0) return;
    sys->temperature = T;
    if (sys->network) sys->network->temperature = T;
}

/* ============================================================================
 * Reaction Network Management
 * ============================================================================ */

SOFEReactionNetwork* sofe_network_create(const char* name) {
    SOFEReactionNetwork* net = (SOFEReactionNetwork*)malloc(
        sizeof(SOFEReactionNetwork));
    assert(net != NULL);
    memset(net, 0, sizeof(SOFEReactionNetwork));
    net->name = (char*)malloc(strlen(name) + 1);
    strcpy(net->name, name);
    net->species_capacity = 8;
    net->species = (SOFESpecies*)malloc(
        net->species_capacity * sizeof(SOFESpecies));
    memset(net->species, 0, net->species_capacity * sizeof(SOFESpecies));
    net->reaction_capacity = 4;
    net->reactions = (SOFEReaction*)malloc(
        net->reaction_capacity * sizeof(SOFEReaction));
    memset(net->reactions, 0, net->reaction_capacity * sizeof(SOFEReaction));
    net->temperature = 298.15;
    net->volume = 1.0;
    net->jacobian_computed = false;
    return net;
}

void sofe_network_free(SOFEReactionNetwork* net) {
    if (!net) return;
    free(net->name);
    for (int i = 0; i < net->n_species; i++) free(net->species[i].name);
    free(net->species);
    for (int i = 0; i < net->n_reactions; i++) {
        free(net->reactions[i].name);
        free(net->reactions[i].reactant_indices);
        free(net->reactions[i].reactant_stoich);
        free(net->reactions[i].product_indices);
        free(net->reactions[i].product_stoich);
    }
    free(net->reactions);
    free(net->stoichiometric_matrix);
    free(net->jacobian);
    free(net->flux_vector);
    free(net->affinity_vector);
    free(net);
}

static void ensure_species_capacity(SOFEReactionNetwork* net) {
    if (net->n_species >= net->species_capacity) {
        net->species_capacity *= 2;
        net->species = (SOFESpecies*)realloc(net->species,
            net->species_capacity * sizeof(SOFESpecies));
        assert(net->species != NULL);
    }
}

static void ensure_reaction_capacity(SOFEReactionNetwork* net) {
    if (net->n_reactions >= net->reaction_capacity) {
        net->reaction_capacity *= 2;
        net->reactions = (SOFEReaction*)realloc(net->reactions,
            net->reaction_capacity * sizeof(SOFEReaction));
        assert(net->reactions != NULL);
    }
}

int sofe_network_add_species(SOFEReactionNetwork* net, const char* name,
                              double init_conc, double diffusion_coeff) {
    if (!net || !name) return -1;
    ensure_species_capacity(net);
    int idx = net->n_species;
    net->species[idx].name = (char*)malloc(strlen(name) + 1);
    strcpy(net->species[idx].name, name);
    net->species[idx].concentration = init_conc;
    net->species[idx].diffusion_coefficient = diffusion_coeff;
    net->species[idx].production_rate = 0.0;
    net->species[idx].degradation_rate = 0.0;
    net->species[idx].equilibrium_value = init_conc;
    net->species[idx].is_autocatalytic = false;
    net->species[idx].is_inhibitor = false;
    net->species[idx].is_activator = false;
    net->species[idx].is_order_parameter = false;
    net->n_species++;
    return idx;
}

int sofe_network_add_reaction(SOFEReactionNetwork* net, const char* name,
                               double rate_const, double activation_energy) {
    if (!net || !name) return -1;
    ensure_reaction_capacity(net);
    int idx = net->n_reactions;
    net->reactions[idx].name = (char*)malloc(strlen(name) + 1);
    strcpy(net->reactions[idx].name, name);
    net->reactions[idx].rate_constant = rate_const;
    net->reactions[idx].activation_energy = activation_energy;
    net->reactions[idx].kinetic_type = 0;
    net->reactions[idx].kinetic_order = 1;
    net->reactions[idx].hill_coefficient = 1.0;
    net->reactions[idx].michaelis_constant = 0.0;
    net->reactions[idx].is_autocatalytic_step = false;
    net->reactions[idx].is_rate_limiting = false;
    net->reactions[idx].n_reactants = 0;
    net->reactions[idx].n_products = 0;
    net->n_reactions++;
    return idx;
}

void sofe_network_set_reactants(SOFEReactionNetwork* net, int rxn_idx,
                                 int* species_indices, double* stoich, int n) {
    if (!net || rxn_idx < 0 || rxn_idx >= net->n_reactions || n <= 0) return;
    SOFEReaction* rxn = &net->reactions[rxn_idx];
    free(rxn->reactant_indices);
    free(rxn->reactant_stoich);
    rxn->n_reactants = n;
    rxn->reactant_indices = (int*)malloc(n * sizeof(int));
    rxn->reactant_stoich = (double*)malloc(n * sizeof(double));
    memcpy(rxn->reactant_indices, species_indices, n * sizeof(int));
    memcpy(rxn->reactant_stoich, stoich, n * sizeof(double));
}

void sofe_network_set_products(SOFEReactionNetwork* net, int rxn_idx,
                                int* species_indices, double* stoich, int n) {
    if (!net || rxn_idx < 0 || rxn_idx >= net->n_reactions || n <= 0) return;
    SOFEReaction* rxn = &net->reactions[rxn_idx];
    free(rxn->product_indices);
    free(rxn->product_stoich);
    rxn->n_products = n;
    rxn->product_indices = (int*)malloc(n * sizeof(int));
    rxn->product_stoich = (double*)malloc(n * sizeof(double));
    memcpy(rxn->product_indices, species_indices, n * sizeof(int));
    memcpy(rxn->product_stoich, stoich, n * sizeof(double));
}

void sofe_network_set_kinetics(SOFEReactionNetwork* net, int rxn_idx,
                                int kinetic_type, int order,
                                double hill_n, double km) {
    if (!net || rxn_idx < 0 || rxn_idx >= net->n_reactions) return;
    net->reactions[rxn_idx].kinetic_type = kinetic_type;
    net->reactions[rxn_idx].kinetic_order = order;
    net->reactions[rxn_idx].hill_coefficient = hill_n;
    net->reactions[rxn_idx].michaelis_constant = km;
}

void sofe_network_set_autocatalytic(SOFEReactionNetwork* net, int rxn_idx,
                                     bool is_auto) {
    if (!net || rxn_idx < 0 || rxn_idx >= net->n_reactions) return;
    net->reactions[rxn_idx].is_autocatalytic_step = is_auto;
}

/* ============================================================================
 * Stoichiometric Matrix Computation
 *
 * The stoichiometric matrix N (n_species x n_reactions) encodes the net
 * production of each species in each reaction:
 *   N[i][r] = product_stoich[i] - reactant_stoich[i]
 *
 * d[X]/dt = N * v  where v is the vector of reaction rates.
 *
 * Reference: Feinberg, Foundations of Chemical Reaction Network Theory (2019)
 * ============================================================================ */

void sofe_network_compute_stoichiometric_matrix(SOFEReactionNetwork* net) {
    if (!net) return;
    int ns = net->n_species;
    int nr = net->n_reactions;
    if (ns == 0 || nr == 0) return;

    free(net->stoichiometric_matrix);
    net->stoichiometric_matrix = (double*)malloc(ns * nr * sizeof(double));
    assert(net->stoichiometric_matrix != NULL);
    memset(net->stoichiometric_matrix, 0, ns * nr * sizeof(double));

    for (int r = 0; r < nr; r++) {
        SOFEReaction* rxn = &net->reactions[r];
        for (int j = 0; j < rxn->n_reactants; j++) {
            int sp = rxn->reactant_indices[j];
            net->stoichiometric_matrix[sp * nr + r] -= rxn->reactant_stoich[j];
        }
        for (int j = 0; j < rxn->n_products; j++) {
            int sp = rxn->product_indices[j];
            net->stoichiometric_matrix[sp * nr + r] += rxn->product_stoich[j];
        }
    }
}

/* ============================================================================
 * Reaction Rate Computation (Mass-Action, Michaelis-Menten, Hill Kinetics)
 *
 * Mass action:   v_r = k_r * prod_i [X_i]^(stoich_i)
 * Michaelis:     v_r = k_r * [S] / (K_m + [S])
 * Hill:          v_r = k_r * [S]^n / (K_d^n + [S]^n)
 *
 * With Arrhenius temperature dependence: k(T) = A * exp(-E_a / RT)
 * ============================================================================ */

double sofe_network_reaction_rate(SOFEReactionNetwork* net, int rxn_idx) {
    if (!net || rxn_idx < 0 || rxn_idx >= net->n_reactions) return 0.0;

    SOFEReaction* rxn = &net->reactions[rxn_idx];
    double rate = rxn->rate_constant;
    const double R = 8.314462618; /* J/(mol*K) */
    if (rxn->activation_energy > 0.0 && net->temperature > 0.0) {
        rate *= exp(-rxn->activation_energy / (R * net->temperature));
    }

    switch (rxn->kinetic_type) {
    case 0: /* Mass action */
        for (int j = 0; j < rxn->n_reactants; j++) {
            int sp = rxn->reactant_indices[j];
            double conc = net->species[sp].concentration;
            double order = rxn->reactant_stoich[j];
            if (conc <= 0.0 && order > 0.0) return 0.0;
            rate *= pow(conc, order);
        }
        break;
    case 1: /* Michaelis-Menten */
        {
            int sp = (rxn->n_reactants > 0) ? rxn->reactant_indices[0] : -1;
            if (sp < 0) return 0.0;
            double conc = net->species[sp].concentration;
            double km = rxn->michaelis_constant;
            if (km <= 0.0) km = 1.0;
            rate *= conc / (km + conc);
        }
        break;
    case 2: /* Hill kinetics */
        {
            int sp = (rxn->n_reactants > 0) ? rxn->reactant_indices[0] : -1;
            if (sp < 0) return 0.0;
            double conc = net->species[sp].concentration;
            double n_hill = rxn->hill_coefficient;
            double km = rxn->michaelis_constant;
            if (km <= 0.0) km = 1.0;
            double conc_n = pow(conc, n_hill);
            double km_n = pow(km, n_hill);
            rate *= conc_n / (km_n + conc_n);
        }
        break;
    default: break;
    }

    return rate;
}

/* ============================================================================
 * Flux Vector Computation
 *
 * J_r = v_r  (reaction flux equals reaction rate)
 * Used for entropy production: sigma = (1/T) * sum_r A_r * J_r >= 0
 * ============================================================================ */

void sofe_network_compute_fluxes(SOFEReactionNetwork* net) {
    if (!net) return;
    int nr = net->n_reactions;
    if (net->flux_vector) free(net->flux_vector);
    net->flux_vector = (double*)malloc(nr * sizeof(double));
    assert(net->flux_vector != NULL);
    for (int r = 0; r < nr; r++) {
        net->flux_vector[r] = sofe_network_reaction_rate(net, r);
    }
}

/* ============================================================================
 * Jacobian Matrix Computation
 *
 * J[i][j] = d(d[X_i]/dt) / d[X_j] = sum_r N[i][r] * (dv_r / d[X_j])
 *
 * For mass-action: dv_r/d[X_j] = (stoich_rj / [X_j]) * v_r  when stoich_rj > 0
 * For Michaelis:   dv_r/d[S] = k * K_m / (K_m + [S])^2
 * For Hill:        dv_r/d[S] = k * n * K_d^n * [S]^(n-1) / (K_d^n + [S]^n)^2
 * ============================================================================ */

void sofe_network_compute_jacobian(SOFEReactionNetwork* net) {
    if (!net) return;
    int ns = net->n_species;
    int nr = net->n_reactions;
    if (ns == 0) return;

    if (!net->stoichiometric_matrix) sofe_network_compute_stoichiometric_matrix(net);

    free(net->jacobian);
    net->jacobian = (double*)malloc(ns * ns * sizeof(double));
    assert(net->jacobian != NULL);
    memset(net->jacobian, 0, ns * ns * sizeof(double));
    if (nr == 0) return;

    for (int r = 0; r < nr; r++) {
        SOFEReaction* rxn = &net->reactions[r];
        double v_r = sofe_network_reaction_rate(net, r);

        for (int j = 0; j < rxn->n_reactants; j++) {
            int sp_j = rxn->reactant_indices[j];
            double c_j = net->species[sp_j].concentration;
            double stoich_j = rxn->reactant_stoich[j];
            double dv_dxj = 0.0;

            if (c_j > 1e-15) {
                switch (rxn->kinetic_type) {
                case 0:
                    dv_dxj = stoich_j * v_r / c_j;
                    break;
                case 1: {
                    double km = rxn->michaelis_constant;
                    if (km <= 0.0) km = 1.0;
                    dv_dxj = rxn->rate_constant * km / ((km + c_j) * (km + c_j));
                    break;
                }
                case 2: {
                    double n_h = rxn->hill_coefficient;
                    double km = rxn->michaelis_constant;
                    if (km <= 0.0) km = 1.0;
                    double cn = pow(c_j, n_h);
                    double kn = pow(km, n_h);
                    double denom = (kn + cn) * (kn + cn);
                    dv_dxj = rxn->rate_constant * n_h * kn * pow(c_j, n_h - 1.0) / denom;
                    break;
                }
                default: break;
                }
            }

            for (int i = 0; i < ns; i++) {
                double N_ir = net->stoichiometric_matrix[i * nr + r];
                net->jacobian[i * ns + sp_j] += N_ir * dv_dxj;
            }
        }
    }
    net->jacobian_computed = true;
}

/* ============================================================================
 * Species Configuration
 * ============================================================================ */

void sofe_species_set_autocatalytic(SOFESpecies* sp, bool is_auto) {
    if (!sp) return;
    sp->is_autocatalytic = is_auto;
}

void sofe_species_set_diffusion(SOFESpecies* sp, double D) {
    if (!sp || D < 0.0) return;
    sp->diffusion_coefficient = D;
}

void sofe_species_set_equilibrium(SOFESpecies* sp, double eq_val) {
    if (!sp || eq_val < 0.0) return;
    sp->equilibrium_value = eq_val;
}

/* ============================================================================
 * Grid Operations — Spatial Discretization for Reaction-Diffusion PDEs
 * ============================================================================ */

SOFEGrid* sofe_grid_create_1d(int nx, double dx, bool periodic) {
    if (nx <= 0 || dx <= 0.0) return NULL;
    SOFEGrid* grid = (SOFEGrid*)malloc(sizeof(SOFEGrid));
    assert(grid != NULL);
    memset(grid, 0, sizeof(SOFEGrid));
    grid->type = SOFE_GRID_1D;
    grid->nx = nx; grid->ny = 1; grid->nz = 1;
    grid->dx = dx; grid->dy = 1.0; grid->dz = 1.0;
    grid->n_points = nx;
    grid->is_periodic_x = periodic;
    grid->max_neighbors = 2;
    return grid;
}

SOFEGrid* sofe_grid_create_2d(int nx, int ny, double dx, double dy,
                               bool periodic_x, bool periodic_y) {
    if (nx <= 0 || ny <= 0 || dx <= 0.0 || dy <= 0.0) return NULL;
    SOFEGrid* grid = (SOFEGrid*)malloc(sizeof(SOFEGrid));
    assert(grid != NULL);
    memset(grid, 0, sizeof(SOFEGrid));
    grid->type = SOFE_GRID_2D;
    grid->nx = nx; grid->ny = ny; grid->nz = 1;
    grid->dx = dx; grid->dy = dy; grid->dz = 1.0;
    grid->n_points = nx * ny;
    grid->is_periodic_x = periodic_x;
    grid->is_periodic_y = periodic_y;
    grid->max_neighbors = 4;
    return grid;
}

SOFEGrid* sofe_grid_create_3d(int nx, int ny, int nz,
                               double dx, double dy, double dz,
                               bool px, bool py, bool pz) {
    if (nx <= 0 || ny <= 0 || nz <= 0 || dx <= 0.0 || dy <= 0.0 || dz <= 0.0)
        return NULL;
    SOFEGrid* grid = (SOFEGrid*)malloc(sizeof(SOFEGrid));
    assert(grid != NULL);
    memset(grid, 0, sizeof(SOFEGrid));
    grid->type = SOFE_GRID_3D;
    grid->nx = nx; grid->ny = ny; grid->nz = nz;
    grid->dx = dx; grid->dy = dy; grid->dz = dz;
    grid->n_points = nx * ny * nz;
    grid->is_periodic_x = px;
    grid->is_periodic_y = py;
    grid->is_periodic_z = pz;
    grid->max_neighbors = 6;
    return grid;
}

void sofe_grid_free(SOFEGrid* grid) {
    if (!grid) return;
    free(grid->coordinates);
    free(grid->neighbor_indices);
    free(grid->neighbor_offsets);
    free(grid);
}

int sofe_grid_index(SOFEGrid* grid, int ix, int iy, int iz) {
    if (!grid) return -1;
    switch (grid->type) {
    case SOFE_GRID_1D: return ix;
    case SOFE_GRID_2D: return ix + iy * grid->nx;
    case SOFE_GRID_3D: return ix + iy * grid->nx + iz * grid->nx * grid->ny;
    default: return -1;
    }
}

/* ============================================================================
 * Grid Neighbor Computation — Discrete Laplacian Stencil
 *
 * 5-point stencil in 2D: nabla^2 u = (u_L+u_R-2u)/dx^2 + (u_D+u_U-2u)/dy^2
 * 7-point stencil in 3D.
 * Supports periodic boundary conditions.
 * ============================================================================ */

void sofe_grid_get_neighbors(SOFEGrid* grid, int p, int* nb, int* nf) {
    if (!grid || !nb || !nf) return;
    *nf = 0;
    int nx = grid->nx, ny = grid->ny, nz = grid->nz;
    switch (grid->type) {
    case SOFE_GRID_1D: {
        int ix = p;
        if (ix > 0) nb[(*nf)++] = ix - 1;
        else if (grid->is_periodic_x) nb[(*nf)++] = nx - 1;
        if (ix < nx - 1) nb[(*nf)++] = ix + 1;
        else if (grid->is_periodic_x) nb[(*nf)++] = 0;
        break;
    }
    case SOFE_GRID_2D: {
        int ix = p % nx, iy = p / nx;
        if (ix > 0) nb[(*nf)++] = p - 1;
        else if (grid->is_periodic_x) nb[(*nf)++] = p - 1 + nx;
        if (ix < nx - 1) nb[(*nf)++] = p + 1;
        else if (grid->is_periodic_x) nb[(*nf)++] = p + 1 - nx;
        if (iy > 0) nb[(*nf)++] = p - nx;
        else if (grid->is_periodic_y) nb[(*nf)++] = p - nx + nx * ny;
        if (iy < ny - 1) nb[(*nf)++] = p + nx;
        else if (grid->is_periodic_y) nb[(*nf)++] = p + nx - nx * ny;
        break;
    }
    case SOFE_GRID_3D: {
        int nxy = nx * ny;
        int iz = p / nxy, rem = p % nxy;
        int iy = rem / nx, ix = rem % nx;
        if (ix > 0) nb[(*nf)++] = p - 1;
        else if (grid->is_periodic_x) nb[(*nf)++] = p - 1 + nx;
        if (ix < nx - 1) nb[(*nf)++] = p + 1;
        else if (grid->is_periodic_x) nb[(*nf)++] = p + 1 - nx;
        if (iy > 0) nb[(*nf)++] = p - nx;
        else if (grid->is_periodic_y) nb[(*nf)++] = p - nx + nxy;
        if (iy < ny - 1) nb[(*nf)++] = p + nx;
        else if (grid->is_periodic_y) nb[(*nf)++] = p + nx - nxy;
        if (iz > 0) nb[(*nf)++] = p - nxy;
        else if (grid->is_periodic_z) nb[(*nf)++] = p - nxy + nxy * nz;
        if (iz < nz - 1) nb[(*nf)++] = p + nxy;
        else if (grid->is_periodic_z) nb[(*nf)++] = p + nxy - nxy * nz;
        break;
    }
    default: break;
    }
}

/* ============================================================================
 * Field Operations — Spatial Concentration Fields
 * ============================================================================ */

SOFEField* sofe_field_create(SOFEGrid* grid, int n_species) {
    if (!grid || n_species <= 0) return NULL;
    SOFEField* field = (SOFEField*)malloc(sizeof(SOFEField));
    assert(field != NULL);
    field->grid = grid;
    field->n_species = n_species;
    field->time = 0.0;
    field->field = (double**)malloc(n_species * sizeof(double*));
    assert(field->field != NULL);
    for (int s = 0; s < n_species; s++) {
        field->field[s] = (double*)malloc(grid->n_points * sizeof(double));
        assert(field->field[s] != NULL);
        memset(field->field[s], 0, grid->n_points * sizeof(double));
    }
    return field;
}

void sofe_field_free(SOFEField* field) {
    if (!field) return;
    if (field->field) {
        for (int s = 0; s < field->n_species; s++) free(field->field[s]);
        free(field->field);
    }
    free(field);
}

void sofe_field_set(SOFEField* field, int species, int point, double val) {
    if (!field || species < 0 || species >= field->n_species ||
        point < 0 || point >= field->grid->n_points) return;
    field->field[species][point] = val;
}

double sofe_field_get(SOFEField* field, int species, int point) {
    if (!field || species < 0 || species >= field->n_species ||
        point < 0 || point >= field->grid->n_points) return 0.0;
    return field->field[species][point];
}

void sofe_field_copy(SOFEField* dst, SOFEField* src) {
    if (!dst || !src || dst->n_species != src->n_species ||
        dst->grid->n_points != src->grid->n_points) return;
    for (int s = 0; s < src->n_species; s++)
        memcpy(dst->field[s], src->field[s],
               src->grid->n_points * sizeof(double));
    dst->time = src->time;
}

/* ============================================================================
 * Discrete Laplacian — 5-point (2D) / 7-point (3D) Stencil
 *
 * Laplacian approximates diffusion: du/dt = D * nabla^2 u
 * Central difference: d^2u/dx^2 ~ [u(x+h) + u(x-h) - 2u(x)] / h^2
 * ============================================================================ */

double sofe_field_laplacian(SOFEField* field, int species, int point) {
    if (!field || species < 0 || species >= field->n_species) return 0.0;
    SOFEGrid* g = field->grid;
    int nb[12]; int nf = 0;
    sofe_grid_get_neighbors(g, point, nb, &nf);
    double cv = field->field[species][point];
    double lap = 0.0;
    switch (g->type) {
    case SOFE_GRID_1D:
        for (int i = 0; i < nf; i++)
            lap += (field->field[species][nb[i]] - cv) / (g->dx * g->dx);
        break;
    case SOFE_GRID_2D:
        for (int i = 0; i < 2 && i < nf; i++)
            lap += (field->field[species][nb[i]] - cv) / (g->dx * g->dx);
        for (int i = 2; i < 4 && i < nf; i++)
            lap += (field->field[species][nb[i]] - cv) / (g->dy * g->dy);
        break;
    case SOFE_GRID_3D:
        for (int i = 0; i < 2 && i < nf; i++)
            lap += (field->field[species][nb[i]] - cv) / (g->dx * g->dx);
        for (int i = 2; i < 4 && i < nf; i++)
            lap += (field->field[species][nb[i]] - cv) / (g->dy * g->dy);
        for (int i = 4; i < 6 && i < nf; i++)
            lap += (field->field[species][nb[i]] - cv) / (g->dz * g->dz);
        break;
    default: break;
    }
    return lap;
}

/* ============================================================================
 * Order Parameter Management (Haken's Synergetics)
 *
 * Order parameters are the slow collective variables that enslave the
 * fast-relaxing modes near a bifurcation. They describe the macroscopic
 * state of the self-organizing system.
 *
 * Amplitude of order parameter: |xi| = projection onto critical eigenmode
 * ============================================================================ */

int sofe_system_add_order_parameter(SOFESystem* sys, const char* name,
                                     SOFEOrderType type,
                                     SOFESymmetryType broken_symmetry,
                                     int n_components) {
    if (!sys || !name || n_components <= 0) return -1;
    if (sys->n_order_params >= sys->order_param_capacity) {
        sys->order_param_capacity *= 2;
        sys->order_params = (SOFEOrderParameter*)realloc(sys->order_params,
            sys->order_param_capacity * sizeof(SOFEOrderParameter));
        assert(sys->order_params != NULL);
    }
    int idx = sys->n_order_params;
    SOFEOrderParameter* op = &sys->order_params[idx];
    memset(op, 0, sizeof(SOFEOrderParameter));
    op->name = (char*)malloc(strlen(name) + 1);
    strcpy(op->name, name);
    op->type = type;
    op->broken_symmetry = broken_symmetry;
    op->n_components = n_components;
    op->amplitude = (double*)malloc(n_components * sizeof(double));
    op->phase = (double*)malloc(n_components * sizeof(double));
    memset(op->amplitude, 0, n_components * sizeof(double));
    memset(op->phase, 0, n_components * sizeof(double));
    op->is_slow_variable = true;
    op->relaxation_time = 1.0;
    sys->n_order_params++;
    return idx;
}

void sofe_system_set_order_amplitude(SOFESystem* sys, int op_idx,
                                      int component, double amplitude) {
    if (!sys || op_idx < 0 || op_idx >= sys->n_order_params) return;
    if (component < 0 || component >= sys->order_params[op_idx].n_components) return;
    sys->order_params[op_idx].amplitude[component] = amplitude;
}

double sofe_system_get_order_amplitude(SOFESystem* sys, int op_idx,
                                        int component) {
    if (!sys || op_idx < 0 || op_idx >= sys->n_order_params) return 0.0;
    if (component < 0 || component >= sys->order_params[op_idx].n_components) return 0.0;
    return sys->order_params[op_idx].amplitude[component];
}

/* ============================================================================
 * Compute Order Parameters from Spatial Field
 *
 * xi(t) = integral{ u*(x) * delta_c(x,t) dx }
 * Implementation: spatial standard deviation from mean as amplitude proxy.
 * ============================================================================ */

void sofe_system_compute_order_parameters(SOFESystem* sys) {
    if (!sys || !sys->concentration_field || sys->n_order_params == 0) return;
    SOFEField* f = sys->concentration_field;
    for (int op = 0; op < sys->n_order_params; op++) {
        SOFEOrderParameter* o = &sys->order_params[op];
        for (int comp = 0; comp < o->n_components; comp++) {
            int sp_idx = comp % f->n_species;
            double sum = 0.0, sum_sq = 0.0;
            int np = f->grid->n_points;
            for (int p = 0; p < np; p++) {
                double val = f->field[sp_idx][p];
                sum += val; sum_sq += val * val;
            }
            double mean = sum / np;
            o->amplitude[comp] = sqrt(sum_sq / np - mean * mean);
            o->susceptibility = o->amplitude[comp] * o->amplitude[comp];
        }
    }
}

/* ============================================================================
 * Control Parameter Management
 * ============================================================================ */

int sofe_system_add_control_param(SOFESystem* sys, const char* name,
                                   double initial_value) {
    if (!sys || !name) return -1;
    int idx = sys->n_control_params++;
    sys->control_params = (double*)realloc(sys->control_params,
        sys->n_control_params * sizeof(double));
    sys->control_param_names = (char**)realloc(sys->control_param_names,
        sys->n_control_params * sizeof(char*));
    sys->control_params[idx] = initial_value;
    sys->control_param_names[idx] = (char*)malloc(strlen(name) + 1);
    strcpy(sys->control_param_names[idx], name);
    if (idx == 0) sys->bifurcation_parameter = initial_value;
    return idx;
}

void sofe_system_set_control_param(SOFESystem* sys, int idx, double value) {
    if (!sys || idx < 0 || idx >= sys->n_control_params) return;
    sys->control_params[idx] = value;
    if (idx == 0) sys->bifurcation_parameter = value;
}

double sofe_system_get_control_param(SOFESystem* sys, int idx) {
    if (!sys || idx < 0 || idx >= sys->n_control_params) return 0.0;
    return sys->control_params[idx];
}

/* ============================================================================
 * History Tracking — Time Series Recording
 * ============================================================================ */

void sofe_system_record_state(SOFESystem* sys) {
    if (!sys) return;
    if (sys->history_len >= sys->history_capacity) {
        sys->history_capacity *= 2;
        sys->history_entropy_prod = (double*)realloc(sys->history_entropy_prod,
            sys->history_capacity * sizeof(double));
        sys->history_order_param = (double*)realloc(sys->history_order_param,
            sys->history_capacity * sizeof(double));
        assert(sys->history_entropy_prod && sys->history_order_param);
    }
    sys->history_entropy_prod[sys->history_len] = sys->entropy_production_rate;
    double amp_sum = 0.0;
    for (int i = 0; i < sys->n_order_params; i++)
        amp_sum += sys->order_params[i].amplitude[0];
    sys->history_order_param[sys->history_len] = amp_sum;
    sys->history_len++;
}

void sofe_system_clear_history(SOFESystem* sys) {
    if (!sys) return;
    sys->history_len = 0;
}

/* ============================================================================
 * Stability Classification from Eigenvalues
 * 
 * Given eigenvalues lambda_i = alpha_i + i*omega_i:
 *   All alpha_i < 0      -> Stable
 *   Any alpha_i > 0      -> Unstable
 *   max(Re) = 0          -> Bifurcation point
 *
 * Classification subtypes: Node (real), Focus (complex), Saddle (mixed signs)
 * ============================================================================ */

SOFEStabilityClass sofe_classify_stability(double* eigenvalues, int n) {
    if (!eigenvalues || n <= 0) return SOFE_STABLE_NODE;
    double max_real = -1e300, min_real = 1e300;
    bool has_imag = false, has_pos = false, has_neg = false, has_zero = false;
    for (int i = 0; i < n; i++) {
        double re = eigenvalues[2 * i];
        double im = eigenvalues[2 * i + 1];
        if (fabs(im) > 1e-12) has_imag = true;
        if (re > max_real) max_real = re;
        if (re < min_real) min_real = re;
        if (re > 1e-10) has_pos = true;
        if (re < -1e-10) has_neg = true;
        if (fabs(re) < 1e-10) has_zero = true;
    }
    if (max_real < -1e-10) return has_imag ? SOFE_STABLE_FOCUS : SOFE_STABLE_NODE;
    if (min_real > 1e-10) return has_imag ? SOFE_UNSTABLE_FOCUS : SOFE_UNSTABLE_NODE;
    if (has_pos && has_neg) return SOFE_SADDLE;
    if (has_zero) return SOFE_CENTER;
    return SOFE_STABLE_NODE;
}

/* ============================================================================
 * String Conversions for Enum Types (L1: Definitions)
 * ============================================================================ */

const char* sofe_stability_name(SOFEStabilityClass sc) {
    static const char* names[] = {
        "Stable Node","Stable Focus","Unstable Node","Unstable Focus",
        "Saddle","Center","Stable Limit Cycle","Unstable Limit Cycle",
        "Strange Attractor"
    };
    return (sc <= 8) ? names[sc] : "Unknown";
}

const char* sofe_bifurcation_name(SOFEBifurcationType bt) {
    static const char* names[] = {
        "None","Saddle-Node","Supercritical Pitchfork","Subcritical Pitchfork",
        "Supercritical Hopf","Subcritical Hopf","Transcritical",
        "Turing","Takens-Bogdanov","Codimension-2"
    };
    return (bt <= 9) ? names[bt] : "Unknown";
}

const char* sofe_regime_name(SOFERegime reg) {
    static const char* names[] = {
        "Equilibrium","Near-Equilibrium","Nonlinear Non-Equilibrium",
        "Far-Equilibrium","Critical","Chaotic Non-Equilibrium"
    };
    return (reg <= 5) ? names[reg] : "Unknown";
}

const char* sofe_symmetry_name(SOFESymmetryType sym) {
    static const char* names[] = {
        "None","Z2","O(2)","SO(3)","Translational","Gauge","Chiral","Scale"
    };
    return (sym <= 7) ? names[sym] : "Unknown";
}

/* ============================================================================
 * System Summary — Human-Readable State Report
 * ============================================================================ */

void sofe_system_print_summary(SOFESystem* sys) {
    if (!sys) { printf("(null system)\n"); return; }
    printf("========================================\n");
    printf("SOFE System: %s\n", sys->name);
    printf("  Regime:            %s\n", sofe_regime_name(sys->regime));
    printf("  Temperature:       %.2f K\n", sys->temperature);
    printf("  Stability:         %s\n", sofe_stability_name(sys->stability));
    if (sys->network) {
        printf("  Species:           %d\n", sys->network->n_species);
        printf("  Reactions:         %d\n", sys->network->n_reactions);
    }
    printf("  Order Parameters:  %d\n", sys->n_order_params);
    printf("  Control Params:    %d\n", sys->n_control_params);
    printf("  Entropy Prod Rate: %.6e\n", sys->entropy_production_rate);
    printf("  Free Energy:       %.6e\n", sys->free_energy);
    printf("  Fluctuation Var:   %.6e\n", sys->fluctuation_variance);
    printf("  History Length:    %d\n", sys->history_len);
    printf("  Macroscopic Order: %s\n",
           sys->has_macroscopic_order ? "Yes" : "No");
    printf("========================================\n");
}
