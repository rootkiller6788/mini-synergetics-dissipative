#include "sofe_core.h"
#include "sofe_patterns.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SOFEPattern* sofe_pattern_detect(SOFEField* field) {
    if (!field) return NULL;
    SOFEPattern* pat = (SOFEPattern*)malloc(sizeof(SOFEPattern));
    assert(pat != NULL);
    memset(pat, 0, sizeof(SOFEPattern));
    int np = field->grid->n_points;
    int ns = field->n_species;
    if (np <= 0 || ns <= 0) { pat->type = SOFE_PATTERN_NONE; return pat; }

    double max_var = 0.0; int max_sp = 0;
    for (int s = 0; s < ns; s++) {
        double sum = 0.0, sum_sq = 0.0;
        for (int p = 0; p < np; p++) {
            double v = field->field[s][p];
            sum += v; sum_sq += v * v;
        }
        double var = sum_sq/np - (sum/np)*(sum/np);
        if (var > max_var) { max_var = var; max_sp = s; }
    }

    pat->amplitude = sqrt(max_var);
    pat->spatial_coherence = (max_var > 0.1) ? 0.8 : 0.1;

    if (max_var < 1e-6) {
        pat->type = SOFE_PATTERN_NONE;
    } else if (field->grid->type == SOFE_GRID_1D) {
        pat->type = SOFE_PATTERN_STRIPES;
        pat->dominant_wavelength = field->grid->dx * np / 4.0;
    } else if (field->grid->type == SOFE_GRID_2D) {
        double cx = 0.0, cy = 0.0, total = 0.0;
        int nx = field->grid->nx, ny = field->grid->ny;
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                int p = ix + iy * nx;
                double v = field->field[max_sp][p];
                cx += ix * v; cy += iy * v; total += v;
            }
        }
        cx /= (total + 1e-30); cy /= (total + 1e-30);
        double radial_var = 0.0;
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                double dx = ix - cx, dy = iy - cy;
                double r = sqrt(dx*dx + dy*dy);
                int p = ix + iy * nx;
                radial_var += field->field[max_sp][p] * r * r;
            }
        }
        radial_var /= (total + 1e-30);
        pat->dominant_wavelength = sqrt(radial_var) * field->grid->dx;
        pat->type = SOFE_PATTERN_SPOTS;
    } else {
        pat->type = SOFE_PATTERN_LABYRINTH;
    }

    pat->temporal_coherence = 0.5;
    return pat;
}
void sofe_pattern_free(SOFEPattern* pat) {
    if (!pat) return;
    free(pat->fourier_amplitudes);
    free(pat);
}

const char* sofe_pattern_name(SOFEPatternType pt) {
    static const char* names[] = {
        "None","Stripes","Spots","Hexagons","Spiral","Target",
        "Labyrinth","Traveling Wave","Standing Wave","Dissipative Soliton"
    };
    return (pt <= 9) ? names[pt] : "Unknown";
}

void sofe_pattern_print(SOFEPattern* pat) {
    if (!pat) { printf("(null pattern)\n"); return; }
    printf("Pattern: %s\n", sofe_pattern_name(pat->type));
    printf("  Amplitude: %.4e  Wavelength: %.4e\n",
           pat->amplitude, pat->dominant_wavelength);
    printf("  Freq: %.4e  Coherence: %.3f/%.3f\n",
           pat->dominant_frequency, pat->spatial_coherence,
           pat->temporal_coherence);
}

SOFETuringAnalysis* sofe_turing_analyze(SOFEReactionDiffusionSystem* rd,
                                          double* diffusivities) {
    if (!rd || !rd->network) return NULL;
    int ns = rd->network->n_species;
    SOFETuringAnalysis* ta = (SOFETuringAnalysis*)malloc(
        sizeof(SOFETuringAnalysis));
    assert(ta != NULL);
    memset(ta, 0, sizeof(SOFETuringAnalysis));
    ta->n_species = ns;
    ta->wavelengths = (double*)malloc(100 * sizeof(double));
    ta->growth_rates = (double*)malloc(100 * sizeof(double));

    if (!rd->network->jacobian) sofe_network_compute_jacobian(rd->network);
    double* J = rd->network->jacobian;

    ta->most_unstable_wavelength = 0.0;
    double max_growth = -1e300;
    int nx = rd->grid->nx;
    double dx = rd->grid->dx;
    double q_min = 2.0 * M_PI / (nx * dx);
    double q_max = M_PI / dx;

    for (int i = 0; i < 100; i++) {
        double q = q_min + (q_max - q_min) * i / 99.0;
        double q2 = q * q;
        double tr = J[0] + J[3] - q2 * (diffusivities[0] + diffusivities[1]);
        double det = (J[0] - q2*diffusivities[0]) * (J[3] - q2*diffusivities[1])
                     - J[1] * J[2];
        double disc = tr*tr/4.0 - det;
        double growth = tr/2.0 + sqrt(disc > 0.0 ? disc : 0.0);
        ta->wavelengths[i] = 2.0 * M_PI / (q + 1e-15);
        ta->growth_rates[i] = growth;
        if (growth > max_growth) {
            max_growth = growth;
            ta->most_unstable_wavelength = ta->wavelengths[i];
        }
        ta->n_modes = i + 1;
    }

    ta->is_turing_unstable = (max_growth > 0.0);
    ta->critical_parameter = ta->most_unstable_wavelength;
    return ta;
}

void sofe_turing_free(SOFETuringAnalysis* ta) {
    if (!ta) return;
    free(ta->wavelengths); free(ta->growth_rates);
    free(ta->activator_diffusivity); free(ta->inhibitor_diffusivity);
    free(ta);
}

void sofe_turing_print(SOFETuringAnalysis* ta) {
    if (!ta) { printf("(null)\n"); return; }
    printf("Turing Analysis:\n");
    printf("  Turing unstable: %s\n", ta->is_turing_unstable ? "Yes" : "No");
    printf("  Unstable wavelength: %.4f\n", ta->most_unstable_wavelength);
    printf("  Max growth rate: %.4e\n",
           ta->growth_rates[ta->n_modes > 0 ? ta->n_modes/2 : 0]);
}

double sofe_turing_threshold(SOFEReactionDiffusionSystem* rd,
                               int activator_idx, int inhibitor_idx) {
    if (!rd || !rd->network) return 0.0;
    double Da = rd->diffusion_coefficients[activator_idx];
    double Di = rd->diffusion_coefficients[inhibitor_idx];
    if (Da >= Di) return -1.0;
    if (!rd->network->jacobian) sofe_network_compute_jacobian(rd->network);
    double* J = rd->network->jacobian;
    int ns = rd->network->n_species;
    double Jaa = J[activator_idx * ns + activator_idx];
    double Jai = J[activator_idx * ns + inhibitor_idx];
    double Jia = J[inhibitor_idx * ns + activator_idx];
    double Jii = J[inhibitor_idx * ns + inhibitor_idx];
    double ratio = Di / Da;
    return ratio * (Jaa * Di + Jii * Da) * (Jaa * Di + Jii * Da)
           / (4.0 * Da * Di * (Jaa * Jii - Jai * Jia) + 1e-30);
}
SOFESpatialCorrelation* sofe_spatial_correlation_compute(SOFEField* field,
                                                           int max_dist) {
    if (!field || max_dist <= 0) return NULL;
    SOFESpatialCorrelation* sc = (SOFESpatialCorrelation*)malloc(
        sizeof(SOFESpatialCorrelation));
    assert(sc != NULL);
    memset(sc, 0, sizeof(SOFESpatialCorrelation));
    sc->max_distance = max_dist;
    sc->n_species = field->n_species;
    sc->correlation = (double**)malloc(field->n_species * sizeof(double*));
    for (int s = 0; s < field->n_species; s++) {
        sc->correlation[s] = (double*)malloc(max_dist * sizeof(double));
        memset(sc->correlation[s], 0, max_dist * sizeof(double));
        int np = field->grid->n_points;
        for (int d = 0; d < max_dist && d < np; d++) {
            double sum = 0.0; int count = 0;
            for (int p = 0; p < np - d; p++) {
                sum += field->field[s][p] * field->field[s][p + d];
                count++;
            }
            sc->correlation[s][d] = (count > 0) ? sum / count : 0.0;
        }
    }
    sc->correlation_length = 0.0;
    double* c0 = sc->correlation[0];
    double norm = c0[0];
    if (norm > 1e-15) {
        for (int d = 1; d < max_dist; d++) {
            if (c0[d] / norm < 0.367879) {
                sc->correlation_length = (double)d; break;
            }
        }
    }
    sc->has_long_range_order = (sc->correlation_length > max_dist * 0.5);
    return sc;
}

void sofe_spatial_correlation_free(SOFESpatialCorrelation* sc) {
    if (!sc) return;
    for (int s = 0; s < sc->n_species; s++) free(sc->correlation[s]);
    free(sc->correlation); free(sc);
}

double sofe_correlation_length_estimate(SOFESpatialCorrelation* sc) {
    return sc ? sc->correlation_length : 0.0;
}

SOFESpatialSpectrum* sofe_spectrum_compute(SOFEField* field, int species) {
    if (!field || species < 0 || species >= field->n_species) return NULL;
    SOFESpatialSpectrum* ss = (SOFESpatialSpectrum*)malloc(
        sizeof(SOFESpatialSpectrum));
    assert(ss != NULL);
    memset(ss, 0, sizeof(SOFESpatialSpectrum));
    int np = field->grid->n_points;
    ss->n_points = np;
    ss->spectrum = (double*)malloc(np * sizeof(double));
    ss->wavenumbers = (double*)malloc(np * sizeof(double));
    for (int k = 0; k < np; k++) {
        double real_part = 0.0, imag_part = 0.0;
        double k_val = 2.0 * M_PI * k / np;
        for (int x = 0; x < np; x++) {
            double theta = k_val * x;
            real_part += field->field[species][x] * cos(theta);
            imag_part -= field->field[species][x] * sin(theta);
        }
        ss->wavenumbers[k] = k_val;
        ss->spectrum[k] = (real_part*real_part + imag_part*imag_part) / np;
        if (ss->spectrum[k] > ss->spectrum[0]) {
            ss->peak_wavenumber = k_val;
        }
        ss->integrated_power += ss->spectrum[k];
    }
    return ss;
}

void sofe_spectrum_free(SOFESpatialSpectrum* ss) {
    if (!ss) return;
    free(ss->spectrum);
    free(ss->wavenumbers);
    free(ss);
}

void sofe_spectrum_peaks(SOFESpatialSpectrum* ss,
                          double** peak_wavenumbers, int* n_peaks) {
    if (!ss || !peak_wavenumbers || !n_peaks) return;
    *n_peaks = 0; *peak_wavenumbers = NULL;
    double threshold = ss->integrated_power / ss->n_points * 3.0;
    for (int i = 1; i < ss->n_points - 1; i++) {
        if (ss->spectrum[i] > threshold &&
            ss->spectrum[i] > ss->spectrum[i-1] &&
            ss->spectrum[i] > ss->spectrum[i+1]) (*n_peaks)++;
    }
    if (*n_peaks == 0) return;
    *peak_wavenumbers = (double*)malloc(*n_peaks * sizeof(double));
    int idx = 0;
    for (int i = 1; i < ss->n_points - 1 && idx < *n_peaks; i++) {
        if (ss->spectrum[i] > threshold &&
            ss->spectrum[i] > ss->spectrum[i-1] &&
            ss->spectrum[i] > ss->spectrum[i+1]) {
            (*peak_wavenumbers)[idx++] = ss->wavenumbers[i];
        }
    }
}

void sofe_field_add_noise(SOFEField* field, double amplitude, int seed) {
    if (!field) return;
    srand((unsigned int)seed);
    for (int s = 0; s < field->n_species; s++) {
        for (int p = 0; p < field->grid->n_points; p++) {
            double r = (double)rand() / RAND_MAX - 0.5;
            field->field[s][p] += amplitude * r;
            if (field->field[s][p] < 0.0) field->field[s][p] = 0.0;
        }
    }
}

void sofe_field_smooth(SOFEField* field, int species, double sigma) {
    if (!field || sigma <= 0.0) return;
    SOFEField* tmp = sofe_field_create(field->grid, 1);
    int np = field->grid->n_points;
    int kernel_size = (int)(3.0 * sigma) + 1;
    for (int p = 0; p < np; p++) {
        double sum = 0.0, weight_sum = 0.0;
        for (int dp = -kernel_size; dp <= kernel_size; dp++) {
            int pp = p + dp;
            if (pp < 0 || pp >= np) continue;
            double w = exp(-dp*dp / (2.0*sigma*sigma));
            sum += field->field[species][pp] * w;
            weight_sum += w;
        }
        tmp->field[0][p] = sum / weight_sum;
    }
    memcpy(field->field[species], tmp->field[0], np * sizeof(double));
    sofe_field_free(tmp);
}

double sofe_field_total_mass(SOFEField* field, int species) {
    if (!field || species < 0 || species >= field->n_species) return 0.0;
    double total = 0.0;
    for (int p = 0; p < field->grid->n_points; p++)
        total += field->field[species][p];
    return total;
}

double sofe_field_spatial_entropy(SOFEField* field, int species, int bins) {
    if (!field || species < 0 || species >= field->n_species || bins <= 0)
        return 0.0;
    int np = field->grid->n_points;
    if (np == 0) return 0.0;
    double min_val = field->field[species][0];
    double max_val = field->field[species][0];
    for (int p = 1; p < np; p++) {
        double v = field->field[species][p];
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
    }
    if (max_val <= min_val) return 0.0;
    double* hist = (double*)calloc(bins, sizeof(double));
    double range = max_val - min_val + 1e-15;
    for (int p = 0; p < np; p++) {
        int bin = (int)((field->field[species][p] - min_val) * bins / range);
        if (bin < 0) bin = 0;
        if (bin >= bins) bin = bins - 1;
        hist[bin] += 1.0;
    }
    double entropy = 0.0;
    for (int b = 0; b < bins; b++) {
        double prob = hist[b] / np;
        if (prob > 1e-15) entropy -= prob * log(prob);
    }
    free(hist);
    return entropy;
}
