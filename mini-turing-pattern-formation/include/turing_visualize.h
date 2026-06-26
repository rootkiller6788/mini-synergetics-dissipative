/**
 * turing_visualize.h — Pattern Visualization & Output
 *
 * Reference: Turing, A.M. (1952) "The Chemical Basis of Morphogenesis"
 *            Pearson, J.E. (1993) Science 261:189-192 (Gray-Scott classification)
 *
 * Provides output functions for Turing pattern simulation results:
 *   - PPM (Portable PixMap) image output for visualization
 *   - Raw binary output for external analysis
 *   - Pattern characterization and classification
 *   - Timeseries extraction from simulation runs
 *
 * Knowledge Coverage:
 *   L5: FFT-based pattern analysis, morphological metrics
 *   L6: Pattern classification (spots, stripes, labyrinths)
 *   L7: Scientific visualization output
 */

#ifndef TURING_VISUALIZE_H
#define TURING_VISUALIZE_H

#include "turing_common.h"
#include <stdio.h>

/* ---------------------------------------------------------------------------
 * Color Map Types
 * -------------------------------------------------------------------------*/
typedef enum {
    COLORMAP_GRAYSCALE  = 0,
    COLORMAP_JET        = 1,
    COLORMAP_HOT        = 2,
    COLORMAP_COOL       = 3,
    COLORMAP_VIRIDIS    = 4,
    COLORMAP_TWILIGHT   = 5,   /* diverging: blue-white-red */
    COLORMAP_CUBEHELIX  = 6    /* perceptually uniform */
} ColorMap;

/* ---------------------------------------------------------------------------
 * PPM Image Output
 *
 * PPM (P6 binary format) is a simple, widely-supported image format.
 * Can be converted to PNG with: convert output.ppm output.png
 * -------------------------------------------------------------------------*/

/**
 * Write a scalar field to a PPM image file.
 *
 * Parameters:
 *   field    — scalar field to visualize
 *   filename — output file path (convention: .ppm extension)
 *   cmap     — color map for value-to-RGB mapping
 *   vmin     — value mapped to minimum color (auto if vmin >= vmax)
 *   vmax     — value mapped to maximum color (auto if vmin >= vmax)
 *
 * Returns:
 *   0 on success, -1 on file error
 *
 * Complexity: O(Nx·Ny)
 */
int field_write_ppm(const ScalarField2D *field, const char *filename,
                    ColorMap cmap, double vmin, double vmax);

/**
 * Write both activator and inhibitor fields side-by-side in one image.
 * Left half: u field, right half: v field.
 */
int field_write_ppm_pair(const ScalarField2D *u, const ScalarField2D *v,
                         const char *filename, ColorMap cmap,
                         double vmin_u, double vmax_u,
                         double vmin_v, double vmax_v);

/**
 * Write RGB composite: R = u, G = v, B = u-v (difference field).
 * This highlights activator-inhibitor relationships.
 */
int field_write_ppm_rgb(const ScalarField2D *u, const ScalarField2D *v,
                        const char *filename);

/* ---------------------------------------------------------------------------
 * Raw Binary Output
 *
 * Single-precision binary files for efficient storage and external analysis.
 * Format: [Nx(int32)][Ny(int32)][Lx(double)][Ly(double)][data(float32)...]
 * -------------------------------------------------------------------------*/

/** Write field to raw binary file. Returns 0 on success. */
int field_write_binary(const ScalarField2D *field, const char *filename);

/** Write both fields to a single binary file. */
int field_write_binary_pair(const ScalarField2D *u, const ScalarField2D *v,
                            const char *filename);

/** Read field from raw binary file. Returns newly allocated field. */
ScalarField2D *field_read_binary(const char *filename);

/* ---------------------------------------------------------------------------
 * Pattern Frame Sequence (for animation)
 *
 * Writes consecutive frames during simulation.
 * Frame filenames: basename_000001.ppm, basename_000002.ppm, etc.
 * -------------------------------------------------------------------------*/

/** Context for multi-frame output. */
typedef struct {
    char     basename[256];     /* base filename */
    char     directory[256];    /* output directory */
    int      frame_count;       /* total frames written */
    int      save_interval;     /* steps between saves */
    ColorMap cmap;
    double   vmin_u, vmax_u;
    double   vmin_v, vmax_v;
} FrameWriter;

/** Initialize a frame writer for sequential output. */
int framewriter_init(FrameWriter *fw, const char *basename,
                     const char *directory, int save_interval);

/** Write one frame (both fields). Returns 0 on success. */
int framewriter_write_frame(FrameWriter *fw, const ReactionDiffusion2D *rd,
                            int step);

/** Write a frame with fixed value range. */
int framewriter_write_frame_scaled(FrameWriter *fw, const ReactionDiffusion2D *rd,
                                   int step, double umin, double umax,
                                   double vmin, double vmax);

/* ---------------------------------------------------------------------------
 * Pattern Characterization
 * -------------------------------------------------------------------------*/

/**
 * Compute pattern metrics from a scalar field.
 *
 * Metrics include:
 *   - Statistical moments (mean, variance, skewness, kurtosis)
 *   - Dominant wavelength from radial FFT power spectrum
 *   - Pattern amplitude (peak-to-peak)
 *   - Pattern type classification via morphological analysis
 *
 * Reference: Cross, M.C. & Greenside, H. (2009) Pattern Formation and Dynamics
 *
 * Complexity: O(Nx·Ny) for moments, O(N log N) for FFT
 */
PatternMetrics field_compute_metrics(const ScalarField2D *field);

/**
 * Classify pattern type based on field structure.
 * Algorithm:
 *   1. Compute variance: if < threshold → HOMOGENEOUS
 *   2. Compute 2D FFT power spectrum
 *   3. If spectrum has sharp ring → SPOTS (isotropic) or STRIPES (anisotropic)
 *   4. If spectrum is broad → LABYRINTHS or CHAOS
 *   5. If hexagonal peaks detected → HEXAGONS
 *
 * Reference: Pearson (1993) classification scheme
 */
PatternType field_classify_pattern(const ScalarField2D *field);

/** Print pattern metrics in human-readable format. */
void pattern_metrics_print(const PatternMetrics *pm);

/** Get pattern type name as string. */
const char *pattern_type_name(PatternType pt);

/* ---------------------------------------------------------------------------
 * Time Series Extraction
 *
 * Extract concentration at specific probe points over simulation time.
 * -------------------------------------------------------------------------*/

/** Probe point for time series. */
typedef struct {
    double x, y;            /* physical coordinates */
    int    i, j;            /* grid indices (computed from coordinates) */
} ProbePoint;

/** Time series context for tracking multiple probe points. */
typedef struct {
    int          max_probes;
    int          n_probes;
    ProbePoint  *probes;
    int          max_steps;
    int          n_steps;
    double      *times;             /* time values, size max_steps */
    double     **u_values;          /* u_values[i_probe][i_step] */
    double     **v_values;          /* v_values[i_probe][i_step] */
} TimeSeries;

/** Allocate time series structure. */
TimeSeries *timeseries_alloc(int max_probes, int max_steps);

/** Free time series structure. */
void timeseries_free(TimeSeries *ts);

/** Add a probe point at physical coordinates (x, y). Returns probe index. */
int timeseries_add_probe(TimeSeries *ts, double x, double y,
                         const ReactionDiffusion2D *rd);

/** Record current state at all probe points. */
int timeseries_record(TimeSeries *ts, const ReactionDiffusion2D *rd, double t);

/** Write time series to CSV file. */
int timeseries_write_csv(const TimeSeries *ts, const char *filename);

/* ---------------------------------------------------------------------------
 * Spatial Correlation Functions
 *
 * Compute 2D spatial autocorrelation to characterize pattern periodicity.
 * -------------------------------------------------------------------------*/

/**
 * Compute 2D spatial autocorrelation of a scalar field.
 *
 * C(rx, ry) = <u(x,y)·u(x+rx, y+ry)> - <u>²
 *
 * Returns the correlation function as a new ScalarField2D.
 * The correlation field has same dimensions as the input.
 *
 * Complexity: O(Nx·Ny) using FFT (Wiener-Khinchin theorem)
 */
ScalarField2D *field_spatial_autocorrelation(const ScalarField2D *field);

/**
 * Compute the radial correlation function C(r) from the 2D autocorrelation.
 * Useful for quantifying pattern periodicity and coherence length.
 *
 * Returns newly allocated arrays r_values and c_values (caller frees).
 * n_bins: number of radial bins.
 */
void field_radial_correlation(const ScalarField2D *field,
                              double **r_values, double **c_values,
                              int n_bins);

/* ---------------------------------------------------------------------------
 * Pattern Comparison
 * -------------------------------------------------------------------------*/

/**
 * Compute the structural similarity between two pattern fields.
 * Returns a value in [0, 1], where 1 = identical.
 *
 * Uses normalized cross-correlation.
 */
double field_pattern_similarity(const ScalarField2D *a, const ScalarField2D *b);

/**
 * Compute the difference map between two fields.
 * result = |a - b| on output field.
 */
int field_difference_map(const ScalarField2D *a, const ScalarField2D *b,
                         ScalarField2D *result);

/* ---------------------------------------------------------------------------
 * Color Mapping Utilities
 * -------------------------------------------------------------------------*/

/**
 * Map a normalized value [0, 1] to RGB using the specified color map.
 * rgb[0]=R, rgb[1]=G, rgb[2]=B, each in [0, 255].
 */
void colormap_apply(ColorMap cmap, double val, unsigned char rgb[3]);

#endif /* TURING_VISUALIZE_H */
