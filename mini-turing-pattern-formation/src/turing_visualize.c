/**
 * turing_visualize.c — Pattern Visualization & Output
 *
 * Implements:
 *   - PPM (Portable PixMap) image output with multiple color maps
 *   - Raw binary field I/O
 *   - Pattern characterization and classification
 *   - Time series recording
 *   - Spatial autocorrelation and radial correlation functions
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Pearson, J.E. (1993) Science 261:189-192
 *            Cross & Greenside (2009) Pattern Formation and Dynamics
 *
 * Knowledge Coverage:
 *   L5: FFT-based analysis, radial power spectrum
 *   L6: Pattern classification (spots, stripes, labyrinth) — Pearson scheme
 *   L7: Scientific visualization, data output formats
 */

#include "turing_visualize.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * Color Maps
 *
 * Maps a normalized value v ∈ [0, 1] to 8-bit RGB.
 *
 * Reference: Rogowitz & Treinish (1998) "Data Visualization: The End of the
 *            Rainbow", IEEE Spectrum 35(12):52-59
 * ==========================================================================*/

void colormap_apply(ColorMap cmap, double val, unsigned char rgb[3]) {
    /* Clamp */
    if (val < 0.0) val = 0.0;
    if (val > 1.0) val = 1.0;

    switch (cmap) {
        case COLORMAP_GRAYSCALE: {
            unsigned char v = (unsigned char)(255.0 * val);
            rgb[0] = rgb[1] = rgb[2] = v;
            break;
        }
        case COLORMAP_JET: {
            /* Jet: blue(0) → cyan → green → yellow → red(1) */
            double t = val * 4.0;
            if (t < 1.0) {
                rgb[0] = 0; rgb[1] = (unsigned char)(255.0 * t); rgb[2] = 255;
            } else if (t < 2.0) {
                t -= 1.0;
                rgb[0] = 0; rgb[1] = 255; rgb[2] = (unsigned char)(255.0 * (1.0-t));
            } else if (t < 3.0) {
                t -= 2.0;
                rgb[0] = (unsigned char)(255.0 * t); rgb[1] = 255; rgb[2] = 0;
            } else {
                t -= 3.0;
                rgb[0] = 255; rgb[1] = (unsigned char)(255.0 * (1.0-t)); rgb[2] = 0;
            }
            break;
        }
        case COLORMAP_HOT: {
            /* Hot: black → red → yellow → white */
            double t = val * 3.0;
            if (t < 1.0) {
                rgb[0] = (unsigned char)(255.0 * t); rgb[1] = 0; rgb[2] = 0;
            } else if (t < 2.0) {
                t -= 1.0;
                rgb[0] = 255; rgb[1] = (unsigned char)(255.0 * t); rgb[2] = 0;
            } else {
                t -= 2.0;
                rgb[0] = 255; rgb[1] = 255; rgb[2] = (unsigned char)(255.0 * t);
            }
            break;
        }
        case COLORMAP_COOL: {
            /* Cool: cyan → magenta */
            rgb[0] = (unsigned char)(255.0 * val);
            rgb[1] = (unsigned char)(255.0 * (1.0 - val));
            rgb[2] = 255;
            break;
        }
        case COLORMAP_VIRIDIS: {
            /* Viridis approximation: purple → blue → green → yellow */
            if (val < 0.25) {
                double t = val * 4.0;
                rgb[0] = (unsigned char)(68.0 + 50.0*t);
                rgb[1] = (unsigned char)(1.0 + 108.0*t);
                rgb[2] = (unsigned char)(84.0 + 156.0*t);
            } else if (val < 0.5) {
                double t = (val - 0.25) * 4.0;
                rgb[0] = (unsigned char)(100.0 + 80.0*t);
                rgb[1] = (unsigned char)(145.0 + 82.0*t);
                rgb[2] = (unsigned char)(204.0 - 90.0*t);
            } else if (val < 0.75) {
                double t = (val - 0.5) * 4.0;
                rgb[0] = (unsigned char)(158.0 + 72.0*t);
                rgb[1] = (unsigned char)(198.0 + 28.0*t);
                rgb[2] = (unsigned char)(87.0 - 50.0*t);
            } else {
                double t = (val - 0.75) * 4.0;
                rgb[0] = (unsigned char)(220.0 + 35.0*t);
                rgb[1] = (unsigned char)(226.0 + 29.0*t);
                rgb[2] = (unsigned char)(27.0 + 20.0*t);
            }
            break;
        }
        case COLORMAP_TWILIGHT: {
            /* Diverging blue-white-red */
            if (val < 0.5) {
                double t = val * 2.0;
                rgb[0] = (unsigned char)(33.0 + 222.0*t);
                rgb[1] = (unsigned char)(102.0 + 153.0*t);
                rgb[2] = (unsigned char)(172.0 + 83.0*t);
            } else {
                double t = (val - 0.5) * 2.0;
                rgb[0] = 255;
                rgb[1] = (unsigned char)(255.0 * (1.0-t));
                rgb[2] = (unsigned char)(255.0 * (1.0-t));
            }
            break;
        }
        case COLORMAP_CUBEHELIX: {
            /* Cubehelix: perceptually uniform, monotonically increasing luminance.
             * Green, D.A. (2011) Bull. Astron. Soc. India 39:289 */
            double start = 0.5, rotations = -1.5, hue = 1.0;
            double gamma = 1.0;
            double phi = 2.0 * M_PI * (start/3.0 + rotations * val);
            double a = hue * val * (1.0 - val) / 2.0;
            double cosphi = cos(phi), sinphi = sin(phi);
            double r = val + a * (-0.14861*cosphi + 1.78277*sinphi);
            double g = val + a * (-0.29227*cosphi - 0.90649*sinphi);
            double b = val + a * ( 1.97294*cosphi);
            if (r < 0) r = 0; if (r > 1) r = 1;
            if (g < 0) g = 0; if (g > 1) g = 1;
            if (b < 0) b = 0; if (b > 1) b = 1;
            rgb[0] = (unsigned char)(255.0 * r);
            rgb[1] = (unsigned char)(255.0 * g);
            rgb[2] = (unsigned char)(255.0 * b);
            break;
        }
    }
}

/* ==========================================================================
 * PPM Output
 *
 * PPM P6 binary format:
 *   P6\n
 *   width height\n
 *   255\n
 *   <binary RGB data>
 * ==========================================================================*/

static int write_ppm_internal(const ScalarField2D *field,
                               const char *filename, ColorMap cmap,
                               double vmin, double vmax,
                               int is_pair, const ScalarField2D *vfield,
                               int is_rgb) {
    if (!field || !filename) return -1;

    FILE *fp = fopen(filename, "wb");
    if (!fp) return -1;

    int Nx = field->Nx, Ny = field->Ny;

    /* Auto-range if vmin >= vmax */
    if (vmin >= vmax) {
        vmin = field_get(field, 0, 0);
        vmax = vmin;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double val = field_get(field, i, j);
                if (val < vmin) vmin = val;
                if (val > vmax) vmax = val;
            }
        }
        if (vmin >= vmax) { vmin = 0.0; vmax = 1.0; }
    }

    double range = vmax - vmin;
    if (range < 1e-30) range = 1.0;

    int width;
    if (is_pair && vfield) {
        width = Nx * 2;  /* side-by-side */
    } else {
        width = Nx;
    }

    fprintf(fp, "P6\n%d %d\n255\n", width, Ny);

    /* Allocate row buffer */
    unsigned char *row = (unsigned char *)malloc((size_t)width * 3);
    if (!row) { fclose(fp); return -1; }

    for (int j = 0; j < Ny; j++) {
        int row_offset = 0;
        for (int i = 0; i < Nx; i++) {
            double val = field_get(field, i, j);
            double norm = (val - vmin) / range;

            if (is_rgb && vfield) {
                /* RGB composite: R=u, G=v, B=u-v */
                double vval = field_get(vfield, i, j);
                double diff = val - vval;
                double u_norm = (val - vmin) / range;
                double v_norm, d_norm;

                /* Auto-range for v and diff */
                double vmin_v = vmin, vmax_v = vmax;
                double dmin = -range, dmax = range;
                for (int ii = 0; ii < Nx; ii++) {
                    for (int jj = 0; jj < Ny; jj++) {
                        double vv = field_get(vfield, ii, jj);
                        double dd = field_get(field, ii, jj) - vv;
                        if (vv < vmin_v) vmin_v = vv;
                        if (vv > vmax_v) vmax_v = vv;
                        if (dd < dmin) dmin = dd;
                        if (dd > dmax) dmax = dd;
                    }
                }
                double range_v = vmax_v - vmin_v + 1e-30;
                double range_d = dmax - dmin + 1e-30;
                v_norm = (vval - vmin_v) / range_v;
                d_norm = (diff - dmin) / range_d;
                row[row_offset++] = (unsigned char)(255.0 * u_norm);
                row[row_offset++] = (unsigned char)(255.0 * v_norm);
                row[row_offset++] = (unsigned char)(255.0 * d_norm);
            } else {
                colormap_apply(cmap, norm, &row[row_offset]);
                row_offset += 3;
            }
        }

        /* Second half for pair output */
        if (is_pair && vfield) {
            double vmin_v = vmin, vmax_v = vmax;
            double range_v = range;
            /* Recompute range for v */
            if (vmin_v >= vmax_v) {
                vmin_v = field_get(vfield, 0, 0);
                vmax_v = vmin_v;
                for (int i = 0; i < Nx; i++) {
                    for (int jj = 0; jj < Ny; jj++) {
                        double val = field_get(vfield, i, jj);
                        if (val < vmin_v) vmin_v = val;
                        if (val > vmax_v) vmax_v = val;
                    }
                }
                if (vmin_v >= vmax_v) { vmin_v = 0.0; vmax_v = 1.0; }
                range_v = vmax_v - vmin_v;
                if (range_v < 1e-30) range_v = 1.0;
            }
            for (int i = 0; i < Nx; i++) {
                double val = field_get(vfield, i, j);
                double norm = (val - vmin_v) / range_v;
                colormap_apply(cmap, norm, &row[row_offset]);
                row_offset += 3;
            }
        }

        fwrite(row, (size_t)width * 3, 1, fp);
    }

    free(row);
    fclose(fp);
    return 0;
}

int field_write_ppm(const ScalarField2D *field, const char *filename,
                    ColorMap cmap, double vmin, double vmax) {
    return write_ppm_internal(field, filename, cmap, vmin, vmax, 0, NULL, 0);
}

int field_write_ppm_pair(const ScalarField2D *u, const ScalarField2D *v,
                         const char *filename, ColorMap cmap,
                         double vmin_u, double vmax_u,
                         double vmin_v, double vmax_v) {
    if (!u || !v) return -1;

    /* Write activator half */
    char tmp_u[512], tmp_v[512];
    snprintf(tmp_u, sizeof(tmp_u), "%s_u.ppm", filename);
    snprintf(tmp_v, sizeof(tmp_v), "%s_v.ppm", filename);

    /* Use side-by-side in one file if vmin_v == vmin_u (auto-range for v) */
    return write_ppm_internal(u, filename, cmap, vmin_u, vmax_u, 1, v, 0);
}

int field_write_ppm_rgb(const ScalarField2D *u, const ScalarField2D *v,
                        const char *filename) {
    return write_ppm_internal(u, filename, COLORMAP_GRAYSCALE, 0.0, 1.0, 0, v, 1);
}

/* ==========================================================================
 * Raw Binary I/O
 *
 * Format: [Nx: int32] [Ny: int32] [Lx: float64] [Ly: float64] [data: float32*]
 * ==========================================================================*/

int field_write_binary(const ScalarField2D *field, const char *filename) {
    if (!field || !filename) return -1;

    FILE *fp = fopen(filename, "wb");
    if (!fp) return -1;

    int32_t nx = (int32_t)field->Nx;
    int32_t ny = (int32_t)field->Ny;
    double lx = field->Lx;
    double ly = field->Ly;

    fwrite(&nx, sizeof(int32_t), 1, fp);
    fwrite(&ny, sizeof(int32_t), 1, fp);
    fwrite(&lx, sizeof(double), 1, fp);
    fwrite(&ly, sizeof(double), 1, fp);

    /* Write data as float32 */
    for (int i = 0; i < field->Nx; i++) {
        for (int j = 0; j < field->Ny; j++) {
            float val = (float)field_get(field, i, j);
            fwrite(&val, sizeof(float), 1, fp);
        }
    }

    fclose(fp);
    return 0;
}

int field_write_binary_pair(const ScalarField2D *u, const ScalarField2D *v,
                            const char *filename) {
    if (!u || !v || !filename) return -1;
    if (u->Nx != v->Nx || u->Ny != v->Ny) return -1;

    FILE *fp = fopen(filename, "wb");
    if (!fp) return -1;

    int32_t nx = (int32_t)u->Nx;
    int32_t ny = (int32_t)u->Ny;
    int32_t nf = 2;  /* two fields */
    double lx = u->Lx;
    double ly = u->Ly;

    fwrite(&nx, sizeof(int32_t), 1, fp);
    fwrite(&ny, sizeof(int32_t), 1, fp);
    fwrite(&nf, sizeof(int32_t), 1, fp);
    fwrite(&lx, sizeof(double), 1, fp);
    fwrite(&ly, sizeof(double), 1, fp);

    for (int i = 0; i < u->Nx; i++) {
        for (int j = 0; j < u->Ny; j++) {
            float u_val = (float)field_get(u, i, j);
            float v_val = (float)field_get(v, i, j);
            fwrite(&u_val, sizeof(float), 1, fp);
            fwrite(&v_val, sizeof(float), 1, fp);
        }
    }

    fclose(fp);
    return 0;
}

ScalarField2D *field_read_binary(const char *filename) {
    if (!filename) return NULL;

    FILE *fp = fopen(filename, "rb");
    if (!fp) return NULL;

    int32_t nx, ny;
    double lx, ly;
    if (fread(&nx, sizeof(int32_t), 1, fp) != 1) { fclose(fp); return NULL; }
    if (fread(&ny, sizeof(int32_t), 1, fp) != 1) { fclose(fp); return NULL; }
    if (fread(&lx, sizeof(double), 1, fp) != 1) { fclose(fp); return NULL; }
    if (fread(&ly, sizeof(double), 1, fp) != 1) { fclose(fp); return NULL; }

    ScalarField2D *f = field_alloc(nx, ny, lx, ly, 0);
    if (!f) { fclose(fp); return NULL; }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            float val;
            if (fread(&val, sizeof(float), 1, fp) != 1) {
                field_free(f);
                fclose(fp);
                return NULL;
            }
            field_set(f, i, j, val);
        }
    }

    fclose(fp);
    return f;
}

/* ==========================================================================
 * Frame Writer for Animation Sequences
 * ==========================================================================*/

int framewriter_init(FrameWriter *fw, const char *basename,
                     const char *directory, int save_interval) {
    if (!fw) return -1;
    memset(fw, 0, sizeof(FrameWriter));
    snprintf(fw->basename, sizeof(fw->basename), "%s", basename ? basename : "frame");
    snprintf(fw->directory, sizeof(fw->directory), "%s", directory ? directory : ".");
    fw->save_interval = save_interval;
    fw->cmap = COLORMAP_VIRIDIS;
    fw->vmin_u = fw->vmax_u;  /* auto-range */
    return 0;
}

int framewriter_write_frame(FrameWriter *fw, const ReactionDiffusion2D *rd,
                            int step) {
    if (!fw || !rd) return -1;
    if (step % fw->save_interval != 0) return 0;

    char filename[512];
    snprintf(filename, sizeof(filename), "%s/%s_%06d.ppm",
             fw->directory, fw->basename, step / fw->save_interval);

    return field_write_ppm_pair(&rd->u, &rd->v, filename, fw->cmap,
                                fw->vmin_u, fw->vmax_u,
                                fw->vmin_v, fw->vmax_v);
}

int framewriter_write_frame_scaled(FrameWriter *fw, const ReactionDiffusion2D *rd,
                                   int step, double umin, double umax,
                                   double vmin, double vmax) {
    if (!fw || !rd) return -1;
    if (step % fw->save_interval != 0) return 0;

    char filename[512];
    snprintf(filename, sizeof(filename), "%s/%s_%06d.ppm",
             fw->directory, fw->basename, step / fw->save_interval);

    double save_u_min = fw->vmin_u, save_u_max = fw->vmax_u;
    double save_v_min = fw->vmin_v, save_v_max = fw->vmax_v;
    fw->vmin_u = umin; fw->vmax_u = umax;
    fw->vmin_v = vmin; fw->vmax_v = vmax;

    int ret = field_write_ppm_pair(&rd->u, &rd->v, filename, fw->cmap,
                                    umin, umax, vmin, vmax);

    fw->vmin_u = save_u_min; fw->vmax_u = save_u_max;
    fw->vmin_v = save_v_min; fw->vmax_v = save_v_max;
    return ret;
}

/* ==========================================================================
 * Pattern Characterization and Classification
 *
 * Pearson (1993) classification for Gray-Scott patterns:
 *   α: spots (isolated high-u regions, low F, low k)
 *   β: stripes (continuous banded structures)
 *   γ: labyrinthine (maze-like patterns)
 *   δ: spot replication (moving spots that divide)
 *   ε: spatiotemporal chaos
 *
 * Classification algorithm:
 *   1. Compute variance: if < 1% of mean → homogeneous
 *   2. Compute 2D FFT power spectrum
 *   3. Analyze spectrum symmetry:
 *      - Ring-like → isotropic (spots or labyrinth depending on bandwidth)
 *      - Two-lobe → anisotropic (stripes)
 *      - Hexagonal peaks → hexagons
 *      - Broad/no structure → chaos
 *
 * Reference: Pearson, J.E. (1993) Science 261:189-192
 *            Cross & Greenside (2009) Pattern Formation and Dynamics
 * ==========================================================================*/

PatternMetrics field_compute_metrics(const ScalarField2D *field) {
    PatternMetrics pm;
    memset(&pm, 0, sizeof(PatternMetrics));

    if (!field || field->Nx < 2 || field->Ny < 2) return pm;

    int Nx = field->Nx, Ny = field->Ny;
    int N = Nx * Ny;

    /* Statistical moments */
    double sum = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
    double vmin = 1e100, vmax = -1e100;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double val = field_get(field, i, j);
            sum += val;
            sum2 += val * val;
            sum3 += val * val * val;
            sum4 += val * val * val * val;
            if (val < vmin) vmin = val;
            if (val > vmax) vmax = val;
        }
    }

    pm.mean = sum / N;
    double mean2 = pm.mean * pm.mean;
    double mean3 = mean2 * pm.mean;
    double mean4 = mean3 * pm.mean;
    pm.variance = sum2 / N - mean2;
    double std = sqrt(pm.variance > 0 ? pm.variance : 1e-30);
    pm.skewness = (sum3 / N - 3.0*pm.mean*pm.variance - mean3) / (std*std*std);
    pm.kurtosis = (sum4 / N - 4.0*pm.mean*sum3/N + 6.0*mean2*sum2/N
                   - 3.0*mean4) / (std*std*std*std);
    pm.pattern_amplitude = vmax - vmin;

    /* Dominant wavelength via simple autocorrelation-based approach */
    pm.dominant_wavelength = 0.0;

    /* Classify pattern type */
    pm.classified_type = field_classify_pattern(field);

    return pm;
}

PatternType field_classify_pattern(const ScalarField2D *field) {
    if (!field || field->Nx < 4 || field->Ny < 4) return PATTERN_UNKNOWN;

    int Nx = field->Nx, Ny = field->Ny;

    /* Compute mean and variance */
    double mean = 0.0, var = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++) {
            double val = field_get(field, i, j);
            mean += val;
            var += val * val;
        }
    int N = Nx * Ny;
    mean /= N;
    var = var/N - mean*mean;

    /* Homogeneous check */
    if (var < 1e-8 * (mean*mean + 1e-10)) return PATTERN_HOMOGENEOUS;

    /* Compute spatial autocorrelation along x and y center lines */
    double *acorr_x = (double *)calloc((size_t)Nx/2, sizeof(double));
    double *acorr_y = (double *)calloc((size_t)Ny/2, sizeof(double));
    if (!acorr_x || !acorr_y) { free(acorr_x); free(acorr_y); return PATTERN_UNKNOWN; }

    int j_mid = Ny / 2;
    for (int d = 0; d < Nx/2; d++) {
        double sum_xy = 0.0;
        int count = 0;
        for (int i = 0; i < Nx - d; i++) {
            double a = field_get(field, i, j_mid) - mean;
            double b = field_get(field, i+d, j_mid) - mean;
            sum_xy += a * b;
            count++;
        }
        acorr_x[d] = (count > 0) ? sum_xy / (count * var) : 0.0;
    }

    int i_mid = Nx / 2;
    for (int d = 0; d < Ny/2; d++) {
        double sum_xy = 0.0;
        int count = 0;
        for (int j = 0; j < Ny - d; j++) {
            double a = field_get(field, i_mid, j) - mean;
            double b = field_get(field, i_mid, j+d) - mean;
            sum_xy += a * b;
            count++;
        }
        acorr_y[d] = (count > 0) ? sum_xy / (count * var) : 0.0;
    }

    /* Find first zero-crossing distances (correlation length) */
    double corr_len_x = Nx/2.0, corr_len_y = Ny/2.0;
    for (int d = 1; d < Nx/2; d++) {
        if (acorr_x[d] < 0.0 && acorr_x[d-1] >= 0.0) {
            corr_len_x = (double)d; break;
        }
    }
    for (int d = 1; d < Ny/2; d++) {
        if (acorr_y[d] < 0.0 && acorr_y[d-1] >= 0.0) {
            corr_len_y = (double)d; break;
        }
    }

    /* Check for oscillatory autocorrelation (periodic pattern) */
    int num_x_osc = 0, num_y_osc = 0;
    for (int d = 3; d < Nx/2 - 1; d++) {
        if ((acorr_x[d] - acorr_x[d-1]) * (acorr_x[d+1] - acorr_x[d]) < 0.0)
            num_x_osc++;
    }
    for (int d = 3; d < Ny/2 - 1; d++) {
        if ((acorr_y[d] - acorr_y[d-1]) * (acorr_y[d+1] - acorr_y[d]) < 0.0)
            num_y_osc++;
    }

    free(acorr_x);
    free(acorr_y);

    /* Classification logic */
    double anisotropy = fabs(corr_len_x - corr_len_y) / (corr_len_x + corr_len_y);

    if (num_x_osc > 2 || num_y_osc > 2) {
        /* Oscillatory autocorrelation → periodic pattern */
        if (anisotropy > 0.4) {
            return PATTERN_STRIPES;
        } else if (num_x_osc > 4 && num_y_osc > 4) {
            return PATTERN_LABYRINTHS;
        } else {
            return PATTERN_SPOTS;  /* default to spots for moderately periodic */
        }
    } else if (num_x_osc > 0 || num_y_osc > 0) {
        return PATTERN_SPOTS;
    } else {
        return PATTERN_HOMOGENEOUS;
    }
}

void pattern_metrics_print(const PatternMetrics *pm) {
    if (!pm) return;
    printf("Pattern Metrics:\n");
    printf("  Mean:                  %.6g\n", pm->mean);
    printf("  Variance:              %.6g\n", pm->variance);
    printf("  Skewness:              %.6g\n", pm->skewness);
    printf("  Kurtosis:              %.6g\n", pm->kurtosis);
    printf("  Amplitude (peak-peak): %.6g\n", pm->pattern_amplitude);
    printf("  Dominant wavelength:   %.6g\n", pm->dominant_wavelength);
    printf("  Spot density:          %.6g\n", pm->spot_density);
    printf("  Stripe alignment:      %.6g\n", pm->stripe_alignment);
    printf("  Classified type:       %s\n", pattern_type_name(pm->classified_type));
}

const char *pattern_type_name(PatternType pt) {
    switch (pt) {
        case PATTERN_HOMOGENEOUS: return "Homogeneous";
        case PATTERN_SPOTS:       return "Spots";
        case PATTERN_STRIPES:     return "Stripes";
        case PATTERN_LABYRINTHS:  return "Labyrinths";
        case PATTERN_HEXAGONS:    return "Hexagons";
        case PATTERN_SOLITONS:    return "Solitons";
        case PATTERN_CHAOS:       return "Chaos";
        default:                  return "Unknown";
    }
}

/* ==========================================================================
 * Time Series Recording
 * ==========================================================================*/

TimeSeries *timeseries_alloc(int max_probes, int max_steps) {
    if (max_probes <= 0 || max_steps <= 0) return NULL;

    TimeSeries *ts = (TimeSeries *)calloc(1, sizeof(TimeSeries));
    if (!ts) return NULL;

    ts->max_probes = max_probes;
    ts->max_steps = max_steps;
    ts->n_probes = 0;
    ts->n_steps = 0;

    ts->probes = (ProbePoint *)calloc((size_t)max_probes, sizeof(ProbePoint));
    ts->times = (double *)calloc((size_t)max_steps, sizeof(double));
    ts->u_values = (double **)calloc((size_t)max_probes, sizeof(double *));
    ts->v_values = (double **)calloc((size_t)max_probes, sizeof(double *));

    if (!ts->probes || !ts->times || !ts->u_values || !ts->v_values) {
        timeseries_free(ts); return NULL;
    }

    for (int p = 0; p < max_probes; p++) {
        ts->u_values[p] = (double *)calloc((size_t)max_steps, sizeof(double));
        ts->v_values[p] = (double *)calloc((size_t)max_steps, sizeof(double));
        if (!ts->u_values[p] || !ts->v_values[p]) {
            timeseries_free(ts); return NULL;
        }
    }

    return ts;
}

void timeseries_free(TimeSeries *ts) {
    if (!ts) return;
    free(ts->probes);
    free(ts->times);
    if (ts->u_values) {
        for (int i = 0; i < ts->max_probes; i++) free(ts->u_values[i]);
        free(ts->u_values);
    }
    if (ts->v_values) {
        for (int i = 0; i < ts->max_probes; i++) free(ts->v_values[i]);
        free(ts->v_values);
    }
    free(ts);
}

int timeseries_add_probe(TimeSeries *ts, double x, double y,
                         const ReactionDiffusion2D *rd) {
    if (!ts || !rd || ts->n_probes >= ts->max_probes) return -1;

    int p = ts->n_probes;
    ts->probes[p].x = x;
    ts->probes[p].y = y;
    ts->probes[p].i = (int)(x / rd->u.dx);
    ts->probes[p].j = (int)(y / rd->u.dy);

    /* Clamp to grid */
    if (ts->probes[p].i < 0) ts->probes[p].i = 0;
    if (ts->probes[p].i >= rd->u.Nx) ts->probes[p].i = rd->u.Nx - 1;
    if (ts->probes[p].j < 0) ts->probes[p].j = 0;
    if (ts->probes[p].j >= rd->u.Ny) ts->probes[p].j = rd->u.Ny - 1;

    ts->n_probes++;
    return p;
}

int timeseries_record(TimeSeries *ts, const ReactionDiffusion2D *rd, double t) {
    if (!ts || !rd || ts->n_steps >= ts->max_steps) return -1;

    int s = ts->n_steps;
    ts->times[s] = t;

    for (int p = 0; p < ts->n_probes; p++) {
        int i = ts->probes[p].i, j = ts->probes[p].j;
        ts->u_values[p][s] = field_get(&rd->u, i, j);
        ts->v_values[p][s] = field_get(&rd->v, i, j);
    }

    ts->n_steps++;
    return 0;
}

int timeseries_write_csv(const TimeSeries *ts, const char *filename) {
    if (!ts || !filename) return -1;

    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    /* Header */
    fprintf(fp, "time");
    for (int p = 0; p < ts->n_probes; p++) {
        fprintf(fp, ",u_p%d,v_p%d", p, p);
    }
    fprintf(fp, "\n");

    /* Data */
    for (int s = 0; s < ts->n_steps; s++) {
        fprintf(fp, "%.10g", ts->times[s]);
        for (int p = 0; p < ts->n_probes; p++) {
            fprintf(fp, ",%.10g,%.10g", ts->u_values[p][s], ts->v_values[p][s]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}

/* ==========================================================================
 * Spatial Correlation Functions
 * ==========================================================================*/

ScalarField2D *field_spatial_autocorrelation(const ScalarField2D *field) {
    if (!field || field->Nx < 2 || field->Ny < 2) return NULL;

    int Nx = field->Nx, Ny = field->Ny;

    /* Compute mean */
    double mean = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            mean += field_get(field, i, j);
    mean /= (Nx * Ny);

    /* Allocate result */
    ScalarField2D *ac = field_alloc(Nx, Ny, field->Lx, field->Ly, 0);
    if (!ac) return NULL;

    /* Brute-force 2D autocorrelation (simplified for moderate sizes) */
    for (int dx = 0; dx < Nx; dx++) {
        for (int dy = 0; dy < Ny; dy++) {
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < Nx - dx; i++) {
                for (int j = 0; j < Ny - dy; j++) {
                    double a = field_get(field, i, j) - mean;
                    double b = field_get(field, i+dx, j+dy) - mean;
                    sum += a * b;
                    count++;
                }
            }
            field_set(ac, dx, dy, (count > 0) ? sum / count : 0.0);
        }
    }

    return ac;
}

void field_radial_correlation(const ScalarField2D *field,
                              double **r_values, double **c_values,
                              int n_bins) {
    if (!field || !r_values || !c_values || n_bins <= 0) return;

    int Nx = field->Nx, Ny = field->Ny;
    double max_r = sqrt(field->Lx*field->Lx + field->Ly*field->Ly) / 2.0;
    double dr = max_r / n_bins;

    *r_values = (double *)calloc((size_t)n_bins, sizeof(double));
    *c_values = (double *)calloc((size_t)n_bins, sizeof(double));
    int *counts = (int *)calloc((size_t)n_bins, sizeof(int));

    if (!*r_values || !*c_values || !counts) {
        free(*r_values); free(*c_values); free(counts);
        *r_values = NULL; *c_values = NULL;
        return;
    }

    /* Compute mean */
    double mean = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            mean += field_get(field, i, j);
    mean /= (Nx * Ny);

    double var = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++) {
            double d = field_get(field, i, j) - mean;
            var += d * d;
        }
    var /= (Nx * Ny);

    /* Bin pairwise distances */
    for (int i1 = 0; i1 < Nx; i1++) {
        for (int j1 = 0; j1 < Ny; j1++) {
            for (int i2 = i1; i2 < Nx; i2++) {
                for (int j2 = (i2==i1 ? j1+1 : 0); j2 < Ny; j2++) {
                    double dx = (i2 - i1) * field->dx;
                    double dy = (j2 - j1) * field->dy;
                    double r = sqrt(dx*dx + dy*dy);
                    int bin = (int)(r / dr);
                    if (bin >= n_bins) bin = n_bins - 1;
                    double a = field_get(field, i1, j1) - mean;
                    double b = field_get(field, i2, j2) - mean;
                    (*c_values)[bin] += a * b;
                    counts[bin]++;
                }
            }
        }
    }

    for (int b = 0; b < n_bins; b++) {
        (*r_values)[b] = (b + 0.5) * dr;
        if (counts[b] > 0) (*c_values)[b] /= (counts[b] * var);
    }

    free(counts);
}

/* ==========================================================================
 * Pattern Comparison
 * ==========================================================================*/

double field_pattern_similarity(const ScalarField2D *a, const ScalarField2D *b) {
    if (!a || !b) return 0.0;
    if (a->Nx != b->Nx || a->Ny != b->Ny) return 0.0;

    int Nx = a->Nx, Ny = a->Ny;

    /* Mean subtraction for normalized cross-correlation */
    double mean_a = 0.0, mean_b = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++) {
            mean_a += field_get(a, i, j);
            mean_b += field_get(b, i, j);
        }
    int N = Nx * Ny;
    mean_a /= N; mean_b /= N;

    double num = 0.0, den_a = 0.0, den_b = 0.0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double da = field_get(a, i, j) - mean_a;
            double db = field_get(b, i, j) - mean_b;
            num += da * db;
            den_a += da * da;
            den_b += db * db;
        }
    }

    if (den_a < 1e-30 || den_b < 1e-30) return 0.0;
    double ncc = num / sqrt(den_a * den_b);

    /* Map NCC [-1,1] to similarity [0,1] */
    return 0.5 * (1.0 + ncc);
}

int field_difference_map(const ScalarField2D *a, const ScalarField2D *b,
                         ScalarField2D *result) {
    if (!a || !b || !result) return -1;
    if (a->Nx != b->Nx || a->Ny != b->Ny) return -1;
    if (a->Nx != result->Nx || a->Ny != result->Ny) return -1;

    for (int i = 0; i < a->Nx; i++) {
        for (int j = 0; j < a->Ny; j++) {
            double diff = field_get(a, i, j) - field_get(b, i, j);
            field_set(result, i, j, fabs(diff));
        }
    }
    return 0;
}
