#ifndef PDS_REACTION_DIFFUSION_H
#define PDS_REACTION_DIFFUSION_H

#include "pds_core.h"

typedef struct {
    int n_species;
    double* diffusion_coeffs;
    int spatial_dims;
    int* grid_size;
    double* domain_length;
    double** u;
    double** du_dt;
    double** laplacian_workspace;
    int n_points;
    double dx;
} PDSReactionDiffusion;

typedef void (*PDSReactionFn)(const double* u, const double* params,
                               int n_params, double* f);

PDSReactionDiffusion* pds_rd_create(int n_species, int spatial_dims,
                                     const int* grid_size,
                                     const double* domain_length,
                                     const double* diffusion_coeffs);
void pds_rd_free(PDSReactionDiffusion* rd);
void pds_rd_set_homogeneous(PDSReactionDiffusion* rd, const double* u_ss);
void pds_rd_compute_laplacian(const PDSReactionDiffusion* rd,
                               int species, double* lap);
void pds_rd_rhs(PDSReactionDiffusion* rd, PDSReactionFn reaction,
                 const double* params, int n_params);
void pds_rd_euler_step(PDSReactionDiffusion* rd, PDSReactionFn reaction,
                        const double* params, int n_params, double dt);
void pds_rd_rk4_step(PDSReactionDiffusion* rd, PDSReactionFn reaction,
                      const double* params, int n_params, double dt);
void pds_rd_add_noise(PDSReactionDiffusion* rd, int species, double amplitude);
double pds_rd_dominant_wavenumber(const PDSReactionDiffusion* rd, int species);
void pds_rd_save_field(const PDSReactionDiffusion* rd, int species,
                        const char* filename);
void pds_rd_print(const PDSReactionDiffusion* rd);

#endif