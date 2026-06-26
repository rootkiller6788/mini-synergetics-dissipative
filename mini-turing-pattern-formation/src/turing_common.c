/**
 * turing_common.c — Core Utilities for Turing Pattern Formation Framework
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Phil. Trans. R. Soc. Lond. B 237(641):37-72
 *
 * Implements:
 *   - ScalarField2D allocation, deallocation, and basic operations
 *   - ReactionDiffusion2D allocation and deallocation
 *   - Boundary condition application (Dirichlet, Neumann, Periodic, Robin)
 *   - Laplacian computation with multiple stencils
 *   - Field initialization patterns (random, Gaussian, sinusoidal, dots, gradient)
 *   - Safe math primitives
 *   - Model parameter defaults
 *
 * Knowledge Coverage:
 *   L1: Field types, boundary conditions, Laplacian discretization
 *   L2: Grid representation, ghost cells
 *   L3: 5-point/7-point/9-point stencils, finite difference operators
 */

#include "turing_common.h"
#include "turing_models.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * Pseudo-Random Number Generator: xoshiro256** (Blackman & Vigna 2018)
 *
 * Fast, statistically robust PRNG for initial condition generation.
 * Period: 2^256 - 1. State: 4 × 64-bit integers.
 *
 * Reference: Blackman, D. & Vigna, S. (2018) "Scrambled Linear Pseudorandom
 *            Number Generators", arXiv:1805.01407
 * ==========================================================================*/

typedef struct {
    uint64_t s[4];
} Xoshiro256SS;

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t xoshiro256ss_next(Xoshiro256SS *state) {
    uint64_t result = rotl64(state->s[1] * 5, 7) * 9;
    uint64_t t = state->s[1] << 17;
    state->s[2] ^= state->s[0];
    state->s[3] ^= state->s[1];
    state->s[1] ^= state->s[2];
    state->s[0] ^= state->s[3];
    state->s[2] ^= t;
    state->s[3] = rotl64(state->s[3], 45);
    return result;
}

static void xoshiro256ss_seed(Xoshiro256SS *state, uint64_t seed) {
    /* SplitMix64 seeding to avoid zero states */
    uint64_t z = seed + 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 4; i++) {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        state->s[i] = z;
    }
    /* Ensure state is non-zero */
    int all_zero = 1;
    for (int i = 0; i < 4; i++) {
        if (state->s[i] != 0) { all_zero = 0; break; }
    }
    if (all_zero) state->s[0] = 1;
}

/* Return uniform double in [0, 1) */
static double xoshiro_uniform(Xoshiro256SS *state) {
    uint64_t x = xoshiro256ss_next(state);
    return (x >> 11) * 0x1.0p-53;  /* 53-bit mantissa precision */
}

/* Box-Muller transform: return standard normal N(0,1) */
static double xoshiro_normal(Xoshiro256SS *state) {
    double u1 = xoshiro_uniform(state);
    double u2 = xoshiro_uniform(state);
    /* Avoid log(0) */
    if (u1 < 1e-16) u1 = 1e-16;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ==========================================================================
 * Safe Math Primitives
 * ==========================================================================*/

double safe_exp(double x) {
    /* Clamp to avoid overflow: exp(709) ≈ 8.2e307, near double max */
    if (x > 709.0) return 1.0e300;
    if (x < -745.0) return 0.0;  /* exp(-745) ≈ 5e-324, near double min */
    return exp(x);
}

double safe_div(double num, double den) {
    if (fabs(den) < 1e-30) return 0.0;
    return num / den;
}

int is_finite_val(double x) {
    return isfinite(x);
}

/* ==========================================================================
 * ScalarField2D Lifecycle
 * ==========================================================================*/

ScalarField2D *field_alloc(int Nx, int Ny, double Lx, double Ly, int ghost) {
    if (Nx <= 0 || Ny <= 0 || Lx <= 0.0 || Ly <= 0.0 || ghost < 0) {
        return NULL;
    }

    ScalarField2D *f = (ScalarField2D *)calloc(1, sizeof(ScalarField2D));
    if (!f) return NULL;

    int total_Nx = Nx + 2 * ghost;
    int total_Ny = Ny + 2 * ghost;
    size_t size = (size_t)total_Nx * (size_t)total_Ny;

    f->data = (double *)calloc(size, sizeof(double));
    if (!f->data) {
        free(f);
        return NULL;
    }

    f->Nx = Nx;
    f->Ny = Ny;
    f->Lx = Lx;
    f->Ly = Ly;
    f->dx = Lx / (double)(Nx - 1);
    f->dy = Ly / (double)(Ny - 1);
    f->ghost = ghost;

    return f;
}

void field_free(ScalarField2D *f) {
    if (f) {
        free(f->data);
        free(f);
    }
}

double *field_interior(ScalarField2D *f) {
    if (!f) return NULL;
    int total_Ny = f->Ny + 2 * f->ghost;
    return f->data + f->ghost * total_Ny + f->ghost;
}

/* ==========================================================================
 * ReactionDiffusion2D Lifecycle
 * ==========================================================================*/

ReactionDiffusion2D *rd_alloc(int Nx, int Ny, double Lx, double Ly,
                               double Du, double Dv, double dt) {
    if (Nx <= 0 || Ny <= 0 || Du < 0.0 || Dv < 0.0 || dt <= 0.0) {
        return NULL;
    }

    ReactionDiffusion2D *rd = (ReactionDiffusion2D *)calloc(1,
                                 sizeof(ReactionDiffusion2D));
    if (!rd) return NULL;

    /* Allocate fields with 1 ghost layer for BC handling */
    ScalarField2D *u = field_alloc(Nx, Ny, Lx, Ly, 1);
    ScalarField2D *v = field_alloc(Nx, Ny, Lx, Ly, 1);

    if (!u || !v) {
        field_free(u);
        field_free(v);
        free(rd);
        return NULL;
    }

    rd->u = *u;
    rd->v = *v;
    free(u);
    free(v);

    rd->Du = Du;
    rd->Dv = Dv;
    rd->dt = dt;
    rd->t = 0.0;
    rd->step = 0;

    return rd;
}

void rd_free(ReactionDiffusion2D *rd) {
    if (rd) {
        free(rd->u.data);
        free(rd->v.data);
        free(rd);
    }
}

/* ==========================================================================
 * Field Basic Operations
 * ==========================================================================*/

void field_fill(ScalarField2D *f, double value) {
    if (!f || !f->data) return;
    int total_Nx = f->Nx + 2 * f->ghost;
    int total_Ny = f->Ny + 2 * f->ghost;
    size_t size = (size_t)total_Nx * (size_t)total_Ny;
    for (size_t i = 0; i < size; i++) {
        f->data[i] = value;
    }
}

void field_fill_interior(ScalarField2D *f, double value) {
    if (!f) return;
    for (int i = 0; i < f->Nx; i++) {
        for (int j = 0; j < f->Ny; j++) {
            field_set(f, i, j, value);
        }
    }
}

int field_copy(const ScalarField2D *src, ScalarField2D *dst) {
    if (!src || !dst || !src->data || !dst->data) return -1;
    if (src->Nx != dst->Nx || src->Ny != dst->Ny) return -1;
    if (src->ghost != dst->ghost) return -1;

    int total_Nx = src->Nx + 2 * src->ghost;
    int total_Ny = src->Ny + 2 * src->ghost;
    size_t size = (size_t)total_Nx * (size_t)total_Ny;
    memcpy(dst->data, src->data, size * sizeof(double));
    return 0;
}

double field_diff_norm(const ScalarField2D *a, const ScalarField2D *b) {
    if (!a || !b) return -1.0;
    if (a->Nx != b->Nx || a->Ny != b->Ny) return -1.0;

    double sum_sq = 0.0;
    for (int i = 0; i < a->Nx; i++) {
        for (int j = 0; j < a->Ny; j++) {
            double diff = field_get(a, i, j) - field_get(b, i, j);
            sum_sq += diff * diff;
        }
    }
    return sqrt(sum_sq / (double)(a->Nx * a->Ny));
}

/* ==========================================================================
 * Boundary Condition Application
 *
 * Populates ghost cells according to boundary condition type.
 *
 * Dirichlet:  u_ghost = fixed_value (mirror: u_ghost = -u_interior for zero)
 * Neumann:    u_ghost = u_interior (zero gradient / no-flux)
 * Periodic:   u_ghost = u_opposite_interior (wrap-around)
 * Robin:      α*u_ghost + β*(u_ghost - u_interior)/dx = γ
 *             → u_ghost = (γ*dx + β*u_interior) / (α*dx + β)
 * ==========================================================================*/

void field_apply_bc(ScalarField2D *f, BoundaryCondition bc_x,
                    BoundaryCondition bc_y) {
    if (!f || f->ghost < 1) return;

    int Nx = f->Nx;
    int Ny = f->Ny;
    int g = f->ghost; /* number of ghost layers */
    int row_stride = Ny + 2 * g;

    for (int layer = 1; layer <= g; layer++) {
        /* === X-boundaries (left and right ghost columns) === */
        for (int j = -g; j < Ny + g; j++) {
            int ghost_left_idx   = (-layer + g) * row_stride + (j + g);
            int ghost_right_idx  = (Nx + layer - 1 + g) * row_stride + (j + g);
            int interior_left_idx  = (layer - 1 + g) * row_stride + (j + g);
            int interior_right_idx = (Nx - layer + g) * row_stride + (j + g);

            /* Periodic BC: need f(-dx)=f(L-dx) and f(L+dx)=f(dx)
             * Left ghost at i=-layer: map to interior at Nx-1-layer
             * Right ghost at i=Nx+layer-1: map to interior at layer */
            int periodic_left_idx  = (g + Nx - 1 - layer) * row_stride + (j + g);
            int periodic_right_idx = (g + layer) * row_stride + (j + g);

            switch (bc_x) {
                case BC_DIRICHLET:
                    /* Zero Dirichlet: u_ghost = -u_interior (odd extension) */
                    f->data[ghost_left_idx]  = -f->data[interior_left_idx];
                    f->data[ghost_right_idx] = -f->data[interior_right_idx];
                    break;
                case BC_NEUMANN:
                    /* Zero Neumann: u_ghost = u_interior (even extension) */
                    f->data[ghost_left_idx]  = f->data[interior_left_idx];
                    f->data[ghost_right_idx] = f->data[interior_right_idx];
                    break;
                case BC_PERIODIC:
                    /* Wrap-around with proper shift: f(-dx)=f(L-dx), f(L+dx)=f(dx) */
                    f->data[ghost_left_idx]  = f->data[periodic_left_idx];
                    f->data[ghost_right_idx] = f->data[periodic_right_idx];
                    break;
                case BC_ROBIN:
                    /* Mixed: αu + β∂u/∂n = 0 → u_ghost = u_interior * (β-α*dx)/(β+α*dx) */
                    /* Using α=1, β=1 for simplicity (can generalize) */
                    {
                        double alpha = 1.0, beta = 1.0;
                        double coeff = (beta - alpha * f->dx) / (beta + alpha * f->dx);
                        f->data[ghost_left_idx]  = coeff * f->data[interior_left_idx];
                        f->data[ghost_right_idx] = coeff * f->data[interior_right_idx];
                    }
                    break;
            }
        }

        /* === Y-boundaries (bottom and top ghost rows) === */
        for (int i = 0; i < Nx; i++) {
            int gi = i + g;
            int ghost_bottom_idx = gi * row_stride + (-layer + g);
            int ghost_top_idx    = gi * row_stride + (Ny + layer - 1 + g);
            int interior_bottom_idx = gi * row_stride + (layer - 1 + g);
            int interior_top_idx    = gi * row_stride + (Ny - layer + g);

            /* Periodic Y: f(-dy)=f(Ly-dy), f(Ly+dy)=f(dy) */
            int periodic_bottom_idx = gi * row_stride + (g + Ny - 1 - layer);
            int periodic_top_idx    = gi * row_stride + (g + layer);

            switch (bc_y) {
                case BC_DIRICHLET:
                    f->data[ghost_bottom_idx] = -f->data[interior_bottom_idx];
                    f->data[ghost_top_idx]    = -f->data[interior_top_idx];
                    break;
                case BC_NEUMANN:
                    f->data[ghost_bottom_idx] = f->data[interior_bottom_idx];
                    f->data[ghost_top_idx]    = f->data[interior_top_idx];
                    break;
                case BC_PERIODIC:
                    /* Wrap-around with proper shift */
                    f->data[ghost_bottom_idx] = f->data[periodic_bottom_idx];
                    f->data[ghost_top_idx]    = f->data[periodic_top_idx];
                    break;
                case BC_ROBIN:
                    {
                        double alpha = 1.0, beta = 1.0;
                        double coeff = (beta - alpha * f->dy) / (beta + alpha * f->dy);
                        f->data[ghost_bottom_idx] = coeff * f->data[interior_bottom_idx];
                        f->data[ghost_top_idx]    = coeff * f->data[interior_top_idx];
                    }
                    break;
            }
        }
    }
}

/* ==========================================================================
 * Laplacian Computation
 *
 * Finite difference approximation of ∇²u = ∂²u/∂x² + ∂²u/∂y².
 *
 * 5-point stencil (2nd order):
 *   ∇²u[i][j] = (u[i+1][j] + u[i-1][j] - 2u[i][j])/dx²
 *             + (u[i][j+1] + u[i][j-1] - 2u[i][j])/dy²
 *
 * 7-point stencil (adds diagonal neighbors for isotropy):
 *   Adds: (u[i+1][j+1] + u[i-1][j-1] + u[i+1][j-1] + u[i-1][j+1] - 4u[i][j])/(4h²)
 *
 * 9-point Mehrstellen stencil (4th order on uniform grid):
 *   ∇²u ≈ (1/(6h²)) * [-20*u[i][j] + 4*sum(neighbors) + sum(diagonals)]
 *
 * Reference: Thomas, J.W. (1995) Numerical PDEs: Finite Difference Methods
 * Reference: Strikwerda, J.C. (2004) Finite Difference Schemes and PDEs
 * ==========================================================================*/

int laplacian_compute(const ScalarField2D *f, ScalarField2D *lap,
                      LaplacianStencil stencil) {
    if (!f || !lap || !f->data || !lap->data) return -1;
    if (f->Nx != lap->Nx || f->Ny != lap->Ny) return -1;
    if (f->ghost != lap->ghost) return -1;

    int Nx = f->Nx;
    int Ny = f->Ny;
    double dx = f->dx;
    double dy = f->dy;
    double idx2 = 1.0 / (dx * dx);
    double idy2 = 1.0 / (dy * dy);

    switch (stencil) {
        case LAPLACIAN_5PT: {
            /* Standard 5-point stencil — O(h²) accuracy */
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    double uc = field_get(f, i, j);
                    double ux = (field_get(f, i+1, j) - 2.0*uc + field_get(f, i-1, j)) * idx2;
                    double uy = (field_get(f, i, j+1) - 2.0*uc + field_get(f, i, j-1)) * idy2;
                    field_set(lap, i, j, ux + uy);
                }
            }
            break;
        }
        case LAPLACIAN_7PT: {
            /* 7-point stencil with diagonal neighbors for improved isotropy */
            double h2 = dx * dx;  /* assume dx ≈ dy */
            double isoh2 = 1.0 / (4.0 * h2);  /* diagonal contribution */
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    double uc = field_get(f, i, j);
                    /* Axis-aligned terms */
                    double ax = (field_get(f, i+1, j) + field_get(f, i-1, j)
                                 - 2.0*uc) * idx2;
                    double ay = (field_get(f, i, j+1) + field_get(f, i, j-1)
                                 - 2.0*uc) * idy2;
                    /* Diagonal terms (weighted 1/4) */
                    double diag = (field_get(f, i+1, j+1) + field_get(f, i-1, j-1)
                                 + field_get(f, i+1, j-1) + field_get(f, i-1, j+1)
                                 - 4.0*uc) * isoh2;
                    field_set(lap, i, j, ax + ay + diag);
                }
            }
            break;
        }
        case LAPLACIAN_9PT: {
            /* 9-point Mehrstellen compact scheme — O(h⁴) accuracy
             * Requires uniform grid (dx = dy = h).
             * ∇²u[i][j] = (1/(6h²)) * [-20*u[i][j]
             *   + 4*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1])
             *   + (u[i+1][j+1] + u[i-1][j-1] + u[i+1][j-1] + u[i-1][j+1])]
             */
            double h = dx;  /* uniform grid assumed */
            double coeff = 1.0 / (6.0 * h * h);
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    double uc = field_get(f, i, j);
                    double sum_neighbors = field_get(f, i+1, j) + field_get(f, i-1, j)
                                         + field_get(f, i, j+1) + field_get(f, i, j-1);
                    double sum_corners   = field_get(f, i+1, j+1) + field_get(f, i-1, j-1)
                                         + field_get(f, i+1, j-1) + field_get(f, i-1, j+1);
                    double lap_val = coeff * (-20.0 * uc + 4.0 * sum_neighbors + sum_corners);
                    field_set(lap, i, j, lap_val);
                }
            }
            break;
        }
    }
    return 0;
}

/* ==========================================================================
 * Spectral Laplacian via Discrete Fourier Transform
 *
 * For periodic boundary conditions, the Laplacian is diagonalized by
 * the Fourier basis: F(∇²u) = -k² · F(u)
 *
 * This implementation uses a simple O(N²logN) row-column FFT.
 * For production use, link with FFTW3 for O(N log N) performance.
 *
 * Reference: Trefethen, L.N. (2000) Spectral Methods in MATLAB
 * ==========================================================================*/

/**
 * Compute 1D FFT of real data, length must be power of 2.
 * in-place Cooley-Tukey radix-2 decimation-in-time.
 *
 * @param real  real part array, length n
 * @param imag  imaginary part array, length n (zeroed on input)
 * @param n     length (must be power of 2)
 * @param inv   0=forward, 1=inverse
 */
static void fft_1d(double *real, double *imag, int n, int inv) {
    /* Bit-reversal permutation */
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double tr = real[i], ti = imag[i];
            real[i] = real[j]; imag[i] = imag[j];
            real[j] = tr; imag[j] = ti;
        }
    }

    /* Butterfly steps */
    for (int len = 2; len <= n; len <<= 1) {
        double angle = 2.0 * M_PI / len * (inv ? 1.0 : -1.0);
        double w_re = cos(angle);
        double w_im = sin(angle);
        for (int i = 0; i < n; i += len) {
            double cur_re = 1.0, cur_im = 0.0;
            for (int j = 0; j < len/2; j++) {
                double u_re = real[i + j];
                double u_im = imag[i + j];
                double v_re = real[i + j + len/2];
                double v_im = imag[i + j + len/2];
                double t_re = cur_re * v_re - cur_im * v_im;
                double t_im = cur_re * v_im + cur_im * v_re;
                real[i + j] = u_re + t_re;
                imag[i + j] = u_im + t_im;
                real[i + j + len/2] = u_re - t_re;
                imag[i + j + len/2] = u_im - t_im;
                double nxt_re = cur_re * w_re - cur_im * w_im;
                double nxt_im = cur_re * w_im + cur_im * w_re;
                cur_re = nxt_re; cur_im = nxt_im;
            }
        }
    }
    if (inv) {
        double ninv = 1.0 / n;
        for (int i = 0; i < n; i++) { real[i] *= ninv; imag[i] *= ninv; }
    }
}

static int next_pow2(int n) {
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

int laplacian_spectral(const ScalarField2D *f, ScalarField2D *lap) {
    if (!f || !lap) return -1;
    if (f->Nx != lap->Nx || f->Ny != lap->Ny) return -1;

    int Nx = f->Nx;
    int Ny = f->Ny;
    int Nx_pad = next_pow2(Nx);
    int Ny_pad = next_pow2(Ny);
    int n_total = Nx_pad * Ny_pad;

    double *real = (double *)calloc((size_t)n_total, sizeof(double));
    double *imag = (double *)calloc((size_t)n_total, sizeof(double));
    if (!real || !imag) { free(real); free(imag); return -1; }

    /* Copy data, zero-pad */
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            real[i * Ny_pad + j] = field_get(f, i, j);

    /* 1D FFTs along rows */
    for (int i = 0; i < Nx_pad; i++) {
        double *row_r = real + i * Ny_pad;
        double *row_i = imag + i * Ny_pad;
        fft_1d(row_r, row_i, Ny_pad, 0);
    }

    /* Transpose + 1D FFTs along columns */
    double *t_real = (double *)calloc((size_t)n_total, sizeof(double));
    double *t_imag = (double *)calloc((size_t)n_total, sizeof(double));
    if (!t_real || !t_imag) {
        free(real); free(imag); free(t_real); free(t_imag); return -1;
    }
    for (int i = 0; i < Nx_pad; i++)
        for (int j = 0; j < Ny_pad; j++) {
            t_real[j * Nx_pad + i] = real[i * Ny_pad + j];
            t_imag[j * Nx_pad + i] = imag[i * Ny_pad + j];
        }
    for (int j = 0; j < Ny_pad; j++) {
        fft_1d(t_real + j * Nx_pad, t_imag + j * Nx_pad, Nx_pad, 0);
    }
    /* Transpose back */
    for (int i = 0; i < Nx_pad; i++)
        for (int j = 0; j < Ny_pad; j++) {
            real[i * Ny_pad + j] = t_real[j * Nx_pad + i];
            imag[i * Ny_pad + j] = t_imag[j * Nx_pad + i];
        }
    free(t_real); free(t_imag);

    /* Apply spectral Laplacian: multiply by -k² */
    double kx_factor = 2.0 * M_PI / f->Lx;
    double ky_factor = 2.0 * M_PI / f->Ly;

    for (int i = 0; i < Nx_pad; i++) {
        int ki = (i <= Nx_pad/2) ? i : (i - Nx_pad);
        double kx = ki * kx_factor;
        for (int j = 0; j < Ny_pad; j++) {
            int kj = (j <= Ny_pad/2) ? j : (j - Ny_pad);
            double ky = kj * ky_factor;
            double k2 = kx * kx + ky * ky;
            double neg_k2 = -k2;
            int idx = i * Ny_pad + j;
            real[idx] *= neg_k2;
            imag[idx] *= neg_k2;
        }
    }

    /* Inverse FFT: rows then columns */
    for (int i = 0; i < Nx_pad; i++) {
        fft_1d(real + i * Ny_pad, imag + i * Ny_pad, Ny_pad, 1);
    }
    t_real = (double *)calloc((size_t)n_total, sizeof(double));
    t_imag = (double *)calloc((size_t)n_total, sizeof(double));
    for (int i = 0; i < Nx_pad; i++)
        for (int j = 0; j < Ny_pad; j++) {
            t_real[j * Nx_pad + i] = real[i * Ny_pad + j];
            t_imag[j * Nx_pad + i] = imag[i * Ny_pad + j];
        }
    for (int j = 0; j < Ny_pad; j++) {
        fft_1d(t_real + j * Nx_pad, t_imag + j * Nx_pad, Nx_pad, 1);
    }
    for (int i = 0; i < Nx_pad; i++)
        for (int j = 0; j < Ny_pad; j++) {
            real[i * Ny_pad + j] = t_real[j * Nx_pad + i];
        }
    free(t_real); free(t_imag);

    /* Copy result back (discard padding) */
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            field_set(lap, i, j, real[i * Ny_pad + j]);

    free(real); free(imag);
    return 0;
}

/* ==========================================================================
 * Initial Conditions
 *
 * All initial condition generators place values in the field interior.
 * Ghost cells are NOT set — call field_apply_bc() after initialization.
 * ==========================================================================*/

void field_init_random(ScalarField2D *f, double mean, double amplitude,
                       uint64_t seed) {
    if (!f) return;
    Xoshiro256SS rng;
    xoshiro256ss_seed(&rng, seed);

    for (int i = 0; i < f->Nx; i++) {
        for (int j = 0; j < f->Ny; j++) {
            double noise = xoshiro_normal(&rng);
            field_set(f, i, j, mean + amplitude * noise);
        }
    }
}

void field_init_gaussian(ScalarField2D *f, double cx, double cy,
                         double sigma, double amplitude) {
    if (!f || sigma <= 0.0) return;
    double two_sigma2 = 2.0 * sigma * sigma;

    for (int i = 0; i < f->Nx; i++) {
        double x = (double)i * f->dx;
        double dx2 = (x - cx) * (x - cx);
        for (int j = 0; j < f->Ny; j++) {
            double y = (double)j * f->dy;
            double dy2 = (y - cy) * (y - cy);
            double r2 = dx2 + dy2;
            field_set(f, i, j, amplitude * exp(-r2 / two_sigma2));
        }
    }
}

void field_init_sinusoidal(ScalarField2D *f, double kx, double ky,
                           double amplitude, double offset) {
    if (!f) return;

    for (int i = 0; i < f->Nx; i++) {
        double x = (double)i * f->dx;
        double sx = cos(kx * x);
        for (int j = 0; j < f->Ny; j++) {
            double y = (double)j * f->dy;
            double sy = cos(ky * y);
            field_set(f, i, j, offset + amplitude * sx * sy);
        }
    }
}

void field_init_random_dots(ScalarField2D *f, int num_dots, double radius,
                            double amplitude, uint64_t seed) {
    if (!f || num_dots <= 0) return;

    /* Clear field first */
    field_fill_interior(f, 0.0);

    Xoshiro256SS rng;
    xoshiro256ss_seed(&rng, seed);

    for (int dot = 0; dot < num_dots; dot++) {
        /* Random center within domain */
        double cx = xoshiro_uniform(&rng) * f->Lx;
        double cy = xoshiro_uniform(&rng) * f->Ly;

        /* Gaussian blob at (cx, cy) */
        double two_sigma2 = 2.0 * radius * radius;
        for (int i = 0; i < f->Nx; i++) {
            double x = (double)i * f->dx;
            double dx2 = (x - cx) * (x - cx);
            if (dx2 > 9.0 * two_sigma2) continue; /* skip far points */
            for (int j = 0; j < f->Ny; j++) {
                double y = (double)j * f->dy;
                double dy2 = (y - cy) * (y - cy);
                double r2 = dx2 + dy2;
                if (r2 < 9.0 * two_sigma2) {
                    double val = field_get(f, i, j);
                    field_set(f, i, j, val + amplitude * exp(-r2 / two_sigma2));
                }
            }
        }
    }
}

void field_init_gradient_x(ScalarField2D *f, double left_val, double right_val) {
    if (!f || f->Nx < 2) return;

    for (int i = 0; i < f->Nx; i++) {
        double frac = (double)i / (double)(f->Nx - 1);
        double val = left_val + frac * (right_val - left_val);
        for (int j = 0; j < f->Ny; j++) {
            field_set(f, i, j, val);
        }
    }
}

/* ==========================================================================
 * Model Parameter Defaults
 *
 * Default parameters are chosen to produce Turing patterns based on
 * published literature values.
 *
 * References:
 *   Gray-Scott: Pearson (1993) Science 261:189-192
 *   Gierer-Meinhardt: Gierer & Meinhardt (1972) Kybernetik 12(1):30-39
 *   Brusselator: Prigogine & Lefever (1968) J Chem Phys 48(4):1695-1700
 *   Schnakenberg: Murray (2003) Mathematical Biology II, §2.4
 *   FitzHugh-Nagumo: FitzHugh (1961) Biophys J 1(6):445-466
 *   Lengyel-Epstein: Lengyel & Epstein (1991) Science 251:650-652
 * ==========================================================================*/

void model_params_default(ModelType type, ModelParams *p) {
    if (!p) return;
    memset(p, 0, sizeof(ModelParams));

    switch (type) {
        case MODEL_GRAY_SCOTT:
            p->F = 0.04;         /* feed rate: 0.02-0.06 for spots/stripes */
            p->k = 0.06;         /* kill rate: 0.055-0.065 for patterns */
            break;

        case MODEL_GIERER_MEINHARDT:
            p->a     = 0.1;      /* activator production */
            p->b     = 1.0;      /* inhibitor production */
            p->c     = 0.001;    /* basal production */
            p->rho_a = 1.0;      /* activator saturation */
            p->rho_h = 1.0;      /* inhibitor saturation */
            p->mu_a  = 0.01;     /* activator decay */
            p->mu_h  = 0.01;     /* inhibitor decay */
            break;

        case MODEL_BRUSSELATOR:
            p->A = 1.0;          /* reservoir concentration of A */
            p->B = 3.0;          /* B > 1 + A² for Turing (here: 3 > 2 ✓) */
            break;

        case MODEL_SCHNAKENBERG:
            p->a = 0.1;          /* dimensionless parameter */
            p->b = 0.9;          /* dimensionless parameter */
            break;

        case MODEL_FITZHUGH_NAGUMO:
            p->epsilon = 0.08;   /* timescale separation */
            p->beta     = 0.7;   /* threshold */
            p->gamma    = 0.5;   /* recovery rate */
            p->I_ext    = 0.0;   /* external stimulus */
            break;

        case MODEL_LENGYEL_EPSTEIN:
            p->a = 12.0;         /* iodide concentration parameter */
            p->b = 1.5;          /* chlorite parameter */
            p->c = 1.0;          /* starch indicator complexation */
            break;

        case MODEL_THOMAS:
            p->a   = 100.0;      /* substrate supply */
            p->b   = 0.5;        /* oxygen supply */
            p->rho = 10.0;       /* reaction rate */
            p->K   = 0.1;        /* inhibition constant */
            break;

        default:
            break;
    }
}

int model_params_validate(ModelType type, const ModelParams *p) {
    if (!p) return 0;

    switch (type) {
        case MODEL_GRAY_SCOTT:
            return (p->F > 0.0 && p->k > 0.0);
        case MODEL_GIERER_MEINHARDT:
            return (p->a > 0.0 && p->b > 0.0 && p->mu_a > 0.0 && p->mu_h > 0.0);
        case MODEL_BRUSSELATOR:
            return (p->A > 0.0 && p->B > 0.0);
        case MODEL_SCHNAKENBERG:
            return (p->a >= 0.0 && p->b >= 0.0);
        case MODEL_FITZHUGH_NAGUMO:
            return (p->epsilon > 0.0 && p->beta >= 0.0 && p->gamma >= 0.0);
        case MODEL_LENGYEL_EPSTEIN:
            return (p->a > 0.0 && p->b > 0.0);
        case MODEL_THOMAS:
            return (p->a > 0.0 && p->b > 0.0 && p->rho > 0.0);
        default:
            return 0;
    }
}

void model_params_print(ModelType type, const ModelParams *p) {
    if (!p) return;

    printf("Model: %s\n", model_name(type));
    printf("Parameters:\n");
    printf("  F=%.6g  k=%.6g  a=%.6g  b=%.6g  c=%.6g\n",
           p->F, p->k, p->a, p->b, p->c);
    printf("  rho_a=%.6g  rho_h=%.6g  mu_a=%.6g  mu_h=%.6g\n",
           p->rho_a, p->rho_h, p->mu_a, p->mu_h);
    printf("  A=%.6g  B=%.6g  epsilon=%.6g  beta=%.6g  gamma=%.6g\n",
           p->A, p->B, p->epsilon, p->beta, p->gamma);
}

/* ==========================================================================
 * Dispersion Relation Lifecycle
 * ==========================================================================*/

DispersionRelation *dispersion_alloc(int n_modes) {
    if (n_modes <= 0) return NULL;

    DispersionRelation *dr = (DispersionRelation *)calloc(1,
                                 sizeof(DispersionRelation));
    if (!dr) return NULL;

    dr->n_modes = n_modes;
    dr->wavenumbers = (double *)calloc((size_t)n_modes, sizeof(double));
    dr->growth_rates = (double *)calloc((size_t)n_modes, sizeof(double));

    if (!dr->wavenumbers || !dr->growth_rates) {
        free(dr->wavenumbers);
        free(dr->growth_rates);
        free(dr);
        return NULL;
    }

    return dr;
}

void dispersion_free(DispersionRelation *dr) {
    if (dr) {
        free(dr->wavenumbers);
        free(dr->growth_rates);
        free(dr);
    }
}
