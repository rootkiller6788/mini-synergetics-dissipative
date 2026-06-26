#include "pds_reaction_diffusion.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static int _idx2d(int i, int j, int ny) { return i * ny + j; }

PDSReactionDiffusion* pds_rd_create(int n_species, int spatial_dims,
                                     const int* grid_size,
                                     const double* domain_length,
                                     const double* diffusion_coeffs) {
    if (n_species <= 0 || spatial_dims < 1 || spatial_dims > 3) return NULL;
    if (!grid_size || !domain_length || !diffusion_coeffs) return NULL;
    PDSReactionDiffusion* rd = calloc(1, sizeof(PDSReactionDiffusion));
    if (!rd) return NULL;
    rd->n_species = n_species;
    rd->spatial_dims = spatial_dims;
    rd->diffusion_coeffs = malloc(n_species * sizeof(double));
    memcpy(rd->diffusion_coeffs, diffusion_coeffs, n_species * sizeof(double));
    rd->grid_size = malloc(spatial_dims * sizeof(int));
    rd->domain_length = malloc(spatial_dims * sizeof(double));
    memcpy(rd->grid_size, grid_size, spatial_dims * sizeof(int));
    memcpy(rd->domain_length, domain_length, spatial_dims * sizeof(double));
    rd->n_points = 1;
    for (int d = 0; d < spatial_dims; d++) rd->n_points *= grid_size[d];
    rd->dx = domain_length[0] / (double)(grid_size[0] - 1);
    rd->u = malloc(n_species * sizeof(double*));
    rd->du_dt = malloc(n_species * sizeof(double*));
    rd->laplacian_workspace = malloc(n_species * sizeof(double*));
    for (int s = 0; s < n_species; s++) {
        rd->u[s] = calloc(rd->n_points, sizeof(double));
        rd->du_dt[s] = calloc(rd->n_points, sizeof(double));
        rd->laplacian_workspace[s] = calloc(rd->n_points, sizeof(double));
    }
    return rd;
}

void pds_rd_free(PDSReactionDiffusion* rd) {
    if (!rd) return;
    free(rd->diffusion_coeffs); free(rd->grid_size); free(rd->domain_length);
    if (rd->u) {
        for (int s = 0; s < rd->n_species; s++) free(rd->u[s]);
        free(rd->u);
    }
    if (rd->du_dt) {
        for (int s = 0; s < rd->n_species; s++) free(rd->du_dt[s]);
        free(rd->du_dt);
    }
    if (rd->laplacian_workspace) {
        for (int s = 0; s < rd->n_species; s++) free(rd->laplacian_workspace[s]);
        free(rd->laplacian_workspace);
    }
    free(rd);
}

void pds_rd_set_homogeneous(PDSReactionDiffusion* rd, const double* u_ss) {
    if (!rd || !u_ss) return;
    for (int s = 0; s < rd->n_species; s++)
        for (int p = 0; p < rd->n_points; p++)
            rd->u[s][p] = u_ss[s];
}

void pds_rd_compute_laplacian(const PDSReactionDiffusion* rd,
                               int species, double* lap) {
    if (!rd || !lap || species < 0 || species >= rd->n_species) return;
    double* u = rd->u[species];
    double dx2 = rd->dx * rd->dx;
    if (rd->spatial_dims == 1) {
        int nx = rd->grid_size[0];
        for (int i = 0; i < nx; i++) {
            double um = (i > 0) ? u[i-1] : u[i];
            double up = (i < nx-1) ? u[i+1] : u[i];
            lap[i] = (up - 2.0 * u[i] + um) / dx2;
        }
    } else if (rd->spatial_dims == 2) {
        int nx = rd->grid_size[0], ny = rd->grid_size[1];
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int idx = _idx2d(i, j, ny);
                double uxm = (i > 0) ? u[_idx2d(i-1, j, ny)] : u[idx];
                double uxp = (i < nx-1) ? u[_idx2d(i+1, j, ny)] : u[idx];
                double uym = (j > 0) ? u[_idx2d(i, j-1, ny)] : u[idx];
                double uyp = (j < ny-1) ? u[_idx2d(i, j+1, ny)] : u[idx];
                lap[idx] = (uxp + uxm + uyp + uym - 4.0 * u[idx]) / dx2;
            }
        }
    }
}

void pds_rd_rhs(PDSReactionDiffusion* rd, PDSReactionFn reaction,
                 const double* params, int n_params) {
    if (!rd || !reaction) return;
    for (int p = 0; p < rd->n_points; p++) {
        double u_local[16], f_local[16];
        for (int s = 0; s < rd->n_species && s < 16; s++)
            u_local[s] = rd->u[s][p];
        reaction(u_local, params, n_params, f_local);
        for (int s = 0; s < rd->n_species && s < 16; s++)
            rd->du_dt[s][p] = f_local[s];
    }
    for (int s = 0; s < rd->n_species; s++) {
        pds_rd_compute_laplacian(rd, s, rd->laplacian_workspace[s]);
        for (int p = 0; p < rd->n_points; p++)
            rd->du_dt[s][p] += rd->diffusion_coeffs[s] * rd->laplacian_workspace[s][p];
    }
}

void pds_rd_euler_step(PDSReactionDiffusion* rd, PDSReactionFn reaction,
                        const double* params, int n_params, double dt) {
    if (!rd || dt <= 0.0) return;
    pds_rd_rhs(rd, reaction, params, n_params);
    for (int s = 0; s < rd->n_species; s++)
        for (int p = 0; p < rd->n_points; p++) {
            rd->u[s][p] += dt * rd->du_dt[s][p];
            if (rd->u[s][p] < 0.0) rd->u[s][p] = 0.0;
        }
}

void pds_rd_add_noise(PDSReactionDiffusion* rd, int species, double amplitude) {
    if (!rd || species < 0 || species >= rd->n_species) return;
    for (int p = 0; p < rd->n_points; p++) {
        double u1 = (double)rand() / (double)RAND_MAX;
        double u2 = (double)rand() / (double)RAND_MAX;
        if (u1 < 1e-10) u1 = 1e-10;
        double gauss = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        rd->u[species][p] += amplitude * gauss;
        if (rd->u[species][p] < 0.0) rd->u[species][p] = 0.0;
    }
}

double pds_rd_dominant_wavenumber(const PDSReactionDiffusion* rd, int species) {
    if (!rd || species < 0 || species >= rd->n_species) return 0.0;
    if (rd->spatial_dims != 1 || rd->n_points < 4) return 0.0;
    double* u = rd->u[species]; int n = rd->n_points;
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += u[i];
    mean /= n;
    int n_crossings = 0;
    for (int i = 1; i < n; i++)
        if ((u[i-1] - mean) * (u[i] - mean) < 0.0) n_crossings++;
    double L = rd->domain_length[0];
    double wavelength = 2.0 * L / (double)(n_crossings > 0 ? n_crossings : 1);
    return 2.0 * M_PI / wavelength;
}

void pds_rd_save_field(const PDSReactionDiffusion* rd, int species,
                        const char* filename) {
    if (!rd || !filename || species < 0 || species >= rd->n_species) return;
    FILE* fp = fopen(filename, "w");
    if (!fp) return;
    if (rd->spatial_dims == 1) {
        fputs("x,u", fp); putc(0x0a, fp);
        for (int i = 0; i < rd->grid_size[0]; i++) {
            double x = (double)i * rd->domain_length[0] / (rd->grid_size[0] - 1);
            fprintf(fp, "%.6f,%.6f", x, rd->u[species][i]);
            putc(0x0a, fp);
        }
    } else if (rd->spatial_dims == 2) {
        fputs("x,y,u", fp); putc(0x0a, fp);
        int nx = rd->grid_size[0], ny = rd->grid_size[1];
        for (int i = 0; i < nx; i++) {
            double x = (double)i * rd->domain_length[0] / (nx - 1);
            for (int j = 0; j < ny; j++) {
                double y = (double)j * rd->domain_length[1] / (ny - 1);
                fprintf(fp, "%.6f,%.6f,%.6f", x, y, rd->u[species][_idx2d(i,j,ny)]);
                putc(0x0a, fp);
            }
        }
    }
    fclose(fp);
}

void pds_rd_print(const PDSReactionDiffusion* rd) {
    if (!rd) { puts("ReactionDiffusion: NULL"); return; }
    puts("=== Reaction-Diffusion System ===");
    printf("  Species: %d, Spatial dims: %d", rd->n_species, rd->spatial_dims);
    putchar(0x0a);
    printf("  Grid: ");
    for (int d = 0; d < rd->spatial_dims; d++) printf("%d ", rd->grid_size[d]);
    printf("(total %d points)", rd->n_points); putchar(0x0a);
    printf("  dx: %.6f", rd->dx); putchar(0x0a);
}
