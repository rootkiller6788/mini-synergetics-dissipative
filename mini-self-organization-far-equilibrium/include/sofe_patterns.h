#ifndef SOFE_PATTERNS_H
#define SOFE_PATTERNS_H

#include "sofe_core.h"
#include "sofe_dynamics.h"

typedef enum {
    SOFE_PATTERN_NONE          = 0,
    SOFE_PATTERN_STRIPES       = 1,
    SOFE_PATTERN_SPOTS         = 2,
    SOFE_PATTERN_HEXAGONS      = 3,
    SOFE_PATTERN_SPIRAL        = 4,
    SOFE_PATTERN_TARGET        = 5,
    SOFE_PATTERN_LABYRINTH     = 6,
    SOFE_PATTERN_TRAVELING_WAVE = 7,
    SOFE_PATTERN_STANDING_WAVE = 8,
    SOFE_PATTERN_DISSIPATIVE_SOLITON = 9
} SOFEPatternType;

typedef struct {
    SOFEPatternType type;
    double dominant_wavelength;
    double dominant_frequency;
    double amplitude;
    double orientation;
    double spatial_coherence;
    double temporal_coherence;
    double* fourier_amplitudes;
    int n_modes;
} SOFEPattern;

typedef struct {
    double* wavelengths;
    double* growth_rates;
    int n_modes;
    double most_unstable_wavelength;
    double critical_parameter;
    bool is_turing_unstable;
    double* activator_diffusivity;
    double* inhibitor_diffusivity;
    int n_species;
} SOFETuringAnalysis;

typedef struct {
    double** correlation;
    int max_distance;
    int n_species;
    double correlation_length;
    bool has_long_range_order;
} SOFESpatialCorrelation;

typedef struct {
    double* spectrum;
    double* wavenumbers;
    int n_points;
    double peak_wavenumber;
    double spectral_width;
    double integrated_power;
} SOFESpatialSpectrum;

SOFEPattern* sofe_pattern_detect(SOFEField* field);
void sofe_pattern_free(SOFEPattern* pat);
const char* sofe_pattern_name(SOFEPatternType pt);
void sofe_pattern_print(SOFEPattern* pat);

SOFETuringAnalysis* sofe_turing_analyze(SOFEReactionDiffusionSystem* rd,
                                          double* diffusivities);
void sofe_turing_free(SOFETuringAnalysis* ta);
void sofe_turing_print(SOFETuringAnalysis* ta);
double sofe_turing_threshold(SOFEReactionDiffusionSystem* rd,
                               int activator_idx, int inhibitor_idx);

SOFESpatialCorrelation* sofe_spatial_correlation_compute(SOFEField* field,
                                                           int max_dist);
void sofe_spatial_correlation_free(SOFESpatialCorrelation* sc);
double sofe_correlation_length_estimate(SOFESpatialCorrelation* sc);

SOFESpatialSpectrum* sofe_spectrum_compute(SOFEField* field, int species);
void sofe_spectrum_free(SOFESpatialSpectrum* ss);
void sofe_spectrum_peaks(SOFESpatialSpectrum* ss,
                          double** peak_wavenumbers, int* n_peaks);

void sofe_field_add_noise(SOFEField* field, double amplitude, int seed);
void sofe_field_smooth(SOFEField* field, int species, double sigma);
void sofe_field_gradient(SOFEField* field, int species,
                          double* grad_x, double* grad_y);
double sofe_field_total_mass(SOFEField* field, int species);
double sofe_field_spatial_entropy(SOFEField* field, int species, int bins);

void sofe_field_apply_perturbation(SOFEField* field,
                                    double (*perturb_func)(double x, double y,
                                                           double z, void* ctx),
                                    void* ctx);

int sofe_field_count_extrema(SOFEField* field, int species,
                              double threshold, int* maxima, int* minima);

#endif
