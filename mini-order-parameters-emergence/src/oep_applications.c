#include "oep_order_parameter.h"
#include "oep_slaving_principle.h"
#include "oep_phase_transition.h"
#include "oep_emergence_metrics.h"
#include "oep_collective_dynamics.h"
#include "oep_landau_ginzburg.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Real-World Applications of Order Parameter Theory
 *
 * This file demonstrates applications of the order parameter
 * framework to concrete systems across multiple domains.
 *
 * L7 Application Areas:
 *   1. Climate: ENSO order parameter from sea surface temperature
 *   2. Neuroscience: Brain coherence as consciousness order parameter
 *   3. Social dynamics: Opinion polarization order parameter
 *   4. Ecology: Ecosystem resilience threshold detection
 *   5. Economics: Market regime transition detection
 *   6. Traffic: Traffic flow order parameter
 *
 * Each function implements a complete, self-contained analysis.
 * ============================================================== */

/* ── Application 1: Climate ENSO Detection ──────────────────── */

typedef struct {
    double* sst_anomalies;    /* Sea surface temperature anomalies (Nino 3.4) */
    int     n_obs;
    double  order_param;      /* ENSO phase: positive = El Nino, negative = La Nina */
    double  threshold;        /* El Nino threshold (typically +0.5°C for 5 months) */
    double  transition_param; /* Control parameter: Walker circulation strength */
    double  critical_value;   /* Bifurcation point estimated from data */
    bool    is_el_nino;
    bool    is_la_nina;
} OEPENSOAnalysis;

OEPENSOAnalysis* oep_app_enso_analyze(const double* sst_data, int n_obs) {
    /* Analyze ENSO (El Niño Southern Oscillation) as an order parameter
     * phenomenon. The equatorial Pacific SST is a macroscopic order
     * parameter that enslaves global atmospheric circulation patterns.
     *
     * Method:
     *   1. Detrend the SST data
     *   2. Identify the order parameter as the leading PC of tropical SST
     *   3. Estimate the bifurcation parameter from the asymmetry
     *      between El Niño and La Niña amplitudes
     *
     * References:
     *   Dijkstra (2005) — Nonlinear Climate Dynamics
     *   Timmermann et al. (2018) — Nature 559:535 */
    if (!sst_data || n_obs < 24) return NULL;

    OEPENSOAnalysis* enso = calloc(1, sizeof(OEPENSOAnalysis));
    if (!enso) return NULL;
    enso->n_obs = n_obs;
    enso->sst_anomalies = malloc((size_t)n_obs * sizeof(double));
    if (!enso->sst_anomalies) { free(enso); return NULL; }

    /* Center data (anomalies) */
    double mean = oep_ts_mean(sst_data, n_obs);
    for (int i = 0; i < n_obs; i++)
        enso->sst_anomalies[i] = sst_data[i] - mean;

    /* Detrend */
    oep_ts_detrend_linear(enso->sst_anomalies, n_obs);

    /* Compute order parameter as smoothed SST anomaly */
    double* smoothed = malloc((size_t)n_obs * sizeof(double));
    if (smoothed) {
        oep_ts_smooth_moving_average(enso->sst_anomalies, smoothed, n_obs, 5);
        enso->order_param = smoothed[n_obs - 1];
        free(smoothed);
    } else {
        enso->order_param = enso->sst_anomalies[n_obs - 1];
    }

    /* Thresholds (NOAA definition) */
    enso->threshold = 0.5;
    enso->is_el_nino = (enso->order_param > enso->threshold);
    enso->is_la_nina = (enso->order_param < -enso->threshold);

    /* Estimate bifurcation parameter from skewness
     * Positive skewness → more strong El Niño events
     * (subcritical character of the ENSO transition) */
    double skew = 0.0, var = oep_ts_variance(enso->sst_anomalies, n_obs);
    if (var > OEP_EPSILON) {
        double std = sqrt(var);
        for (int i = 0; i < n_obs; i++) {
            double z = enso->sst_anomalies[i] / std;
            skew += z * z * z;
        }
        skew /= (double)n_obs;
    }
    enso->transition_param = skew;
    enso->critical_value = 0.0;

    return enso;
}

void oep_app_enso_free(OEPENSOAnalysis* enso) {
    if (!enso) return;
    free(enso->sst_anomalies);
    free(enso);
}

/* ── Application 2: Brain Coherence (Consciousness) ─────────── */

typedef struct {
    double* eeg_channels;    /* Multi-channel EEG data (row-major: time × channels) */
    int     n_channels;
    int     n_samples;
    double  kuramoto_order;  /* Kuramoto coherence r across channels */
    double  phi_integrated;  /* Integrated information Φ */
    double  complexity;      /* Lempel-Ziv complexity proxy */
    bool    is_conscious;    /* Binary classification */
} OEPBrainCoherence;

OEPBrainCoherence* oep_app_brain_coherence_analyze(const double* eeg_data,
                                                    int n_channels, int n_samples) {
    /* Analyze brain coherence as an order parameter for consciousness.
     *
     * The Kuramoto order parameter r across EEG channels measures
     * the degree of phase synchronization. High r in specific frequency
     * bands is associated with conscious processing.
     *
     * References:
     *   Tononi et al. (2016) — Nature Reviews Neuroscience 17:450
     *   Casarotto et al. (2016) — Annals of Neurology 80:718 */
    if (!eeg_data || n_channels < 2 || n_samples < 100) return NULL;

    OEPBrainCoherence* bc = calloc(1, sizeof(OEPBrainCoherence));
    if (!bc) return NULL;

    bc->n_channels = n_channels;
    bc->n_samples = n_samples;
    bc->eeg_channels = malloc((size_t)(n_channels * n_samples) * sizeof(double));
    if (!bc->eeg_channels) { free(bc); return NULL; }
    memcpy(bc->eeg_channels, eeg_data,
           (size_t)(n_channels * n_samples) * sizeof(double));

    /* Estimate phase coherence using Kuramoto order parameter */
    OEPKuramotoSystem* kura = oep_kura_create(n_channels);
    if (kura) {
        oep_kura_set_frequencies_distribution(kura, 0.0, 1.0);
        oep_kura_set_coupling(kura, 2.0);

        /* Use EEG channel means as proxy for oscillator phases */
        for (int ch = 0; ch < n_channels; ch++) {
            double ch_mean = 0.0;
            for (int t = 0; t < n_samples; t++)
                ch_mean += eeg_data[t * n_channels + ch];
            ch_mean /= (double)n_samples;
            kura->oscillators[ch].phase = ch_mean * OEP_PI;
            while (kura->oscillators[ch].phase > 2.0 * OEP_PI)
                kura->oscillators[ch].phase -= 2.0 * OEP_PI;
            while (kura->oscillators[ch].phase < 0.0)
                kura->oscillators[ch].phase += 2.0 * OEP_PI;
        }

        oep_kura_compute_order_parameter(kura);
        bc->kuramoto_order = kura->mean_field_coherence;
        oep_kura_free(kura);
    }

    /* Estimate integrated information (simplified) */
    bc->phi_integrated = bc->kuramoto_order;  /* Surrogate */

    /* Complexity: entropy of the channel correlation matrix */
    double* corr_matrix = malloc((size_t)(n_channels * n_channels) * sizeof(double));
    if (corr_matrix) {
        for (int i = 0; i < n_channels; i++) {
            for (int j = 0; j < n_channels; j++) {
                double* ch_i = malloc((size_t)n_samples * sizeof(double));
                double* ch_j = malloc((size_t)n_samples * sizeof(double));
                if (ch_i && ch_j) {
                    for (int t = 0; t < n_samples; t++) {
                        ch_i[t] = eeg_data[t * n_channels + i];
                        ch_j[t] = eeg_data[t * n_channels + j];
                    }
                    corr_matrix[i * n_channels + j] =
                        0.5 * (1.0 + oep_ts_cross_correlation(ch_i, ch_j, n_samples, 0));
                }
                free(ch_i); free(ch_j);
            }
        }
        /* Eigenvalue entropy of correlation matrix as complexity proxy */
        double sum_eig = 0.0;
        for (int i = 0; i < n_channels; i++) sum_eig += corr_matrix[i * n_channels + i];
        double complexity = 0.0;
        for (int i = 0; i < n_channels; i++) {
            double p = (sum_eig > OEP_EPSILON) ? corr_matrix[i * n_channels + i] / sum_eig : 0.0;
            if (p > 0.0) complexity -= p * log(p);
        }
        bc->complexity = complexity / log((double)n_channels);
        free(corr_matrix);
    }

    /* Consciousness threshold: high coherence + high complexity */
    bc->is_conscious = (bc->kuramoto_order > 0.3) && (bc->complexity > 0.4);

    return bc;
}

void oep_app_brain_free(OEPBrainCoherence* bc) {
    if (!bc) return;
    free(bc->eeg_channels);
    free(bc);
}

/* ── Application 3: Opinion Polarization ────────────────────── */

typedef struct {
    double polarization_order;    /* OP: 0 = consensus, 1 = complete polarization */
    double control_param;         /* p: echo chamber strength */
    double critical_value;        /* pc: estimated critical echo chamber */
    bool   is_polarized;
    int*   opinion_counts;
    int    n_bins;
} OEPOpinionPolarization;

OEPOpinionPolarization* oep_app_opinion_polarization_analyze(const double* opinions,
                                                               int n_people) {
    /* Analyze opinion polarization as an order parameter phenomenon.
     *
     * The distribution of opinions can be unimodal (consensus) or
     * bimodal (polarized). The bifurcation is driven by the
     * "echo chamber" strength — the tendency of people to interact
     * only with like-minded others.
     *
     * Model:
     *   dξ/dt = α(p)ξ - βξ³
     * where ξ is the polarization order parameter (distance from median),
     * α(p) ∝ p - p_c (p = echo chamber strength).
     *
     * References:
     *   Dandekar et al. (2013) — Proc. Natl. Acad. Sci. 110:5791
     *   Baumann et al. (2020) — Nature Human Behaviour 4:1091 */
    if (!opinions || n_people < 10) return NULL;

    OEPOpinionPolarization* pol = calloc(1, sizeof(OEPOpinionPolarization));
    if (!pol) return NULL;

    /* Compute polarization as normalized variance */
    double mu = oep_ts_mean(opinions, n_people);
    double var = oep_ts_variance(opinions, n_people);

    /* Bimodality coefficient:
     *   BC = (skewness² + 1) / (kurtosis + 3·(n-1)²/((n-2)(n-3)))
     * BC > 0.555 suggests bimodality. */
    double skew = 0.0, kurt = 0.0;
    if (var > OEP_EPSILON) {
        double std = sqrt(var);
        for (int i = 0; i < n_people; i++) {
            double z = (opinions[i] - mu) / std;
            double z2 = z * z;
            skew += z * z2;
            kurt += z2 * z2;
        }
        skew /= (double)n_people;
        kurt /= (double)n_people;
    }
    double n_d = (double)n_people;
    double bc = (skew * skew + 1.0)
                / (kurt + 3.0 * (n_d - 1.0) * (n_d - 1.0)
                   / ((n_d - 2.0) * (n_d - 3.0)));

    pol->polarization_order = bc;
    pol->is_polarized = (bc > 0.555);

    /* Estimate echo chamber strength from polarization */
    pol->control_param = bc;
    pol->critical_value = 0.555;

    return pol;
}

void oep_app_opinion_free(OEPOpinionPolarization* pol) {
    if (!pol) return;
    free(pol->opinion_counts);
    free(pol);
}

/* ── Application 4: Ecosystem Resilience ────────────────────── */

typedef struct {
    double resilience_order;   /* OP: distance to tipping point */
    double recovery_rate;      /* λ: linear recovery rate (λ < 0 → unstable) */
    double autocorr_lag1;      /* AR(1) coefficient (1 → critical) */
    double variance;           /* Increasing variance → flickering */
    bool   is_resilient;
    bool   is_tipping;
} OEPEcosystemResilience;

OEPEcosystemResilience* oep_app_ecosystem_resilience_analyze(const double* biomass,
                                                               int n_obs) {
    /* Analyze ecosystem resilience using early warning signals (EWS)
     * derived from order parameter theory.
     *
     * As an ecosystem approaches a tipping point (e.g., desertification,
     * fishery collapse, lake eutrophication), the dominant eigenvalue
     * of the linearized dynamics approaches zero from below:
     *
     *   λ → 0⁻ (critical slowing down)
     *
     * EWS indicators:
     *   1. Increasing AR(1) coefficient → approaching 1
     *   2. Increasing variance (flickering)
     *   3. Increasing autocorrelation
     *
     * References:
     *   Scheffer et al. (2009) — Nature 461:53
     *   Dakos et al. (2008) — Proc. Natl. Acad. Sci. 105:14308 */
    if (!biomass || n_obs < 30) return NULL;

    OEPEcosystemResilience* eco = calloc(1, sizeof(OEPEcosystemResilience));
    if (!eco) return NULL;

    /* Detrend */
    double* detrended = malloc((size_t)n_obs * sizeof(double));
    if (!detrended) { free(eco); return NULL; }
    memcpy(detrended, biomass, (size_t)n_obs * sizeof(double));
    oep_ts_detrend_linear(detrended, n_obs);

    /* AR(1) coefficient (lag-1 autocorrelation) */
    eco->autocorr_lag1 = oep_ts_autocorrelation(detrended, n_obs, 1);

    /* Variance */
    eco->variance = oep_ts_variance(detrended, n_obs);

    /* Recovery rate from AR(1):
     *   λ = log(AR1) / Δt  (discrete-time approximation)
     * λ < 0: stable (returns to equilibrium)
     * λ → 0: critical (approaching tipping point) */
    if (eco->autocorr_lag1 > 0.0 && eco->autocorr_lag1 < 1.0)
        eco->recovery_rate = log(eco->autocorr_lag1);
    else if (eco->autocorr_lag1 <= 0.0)
        eco->recovery_rate = -1e10;  /* Very stable */
    else
        eco->recovery_rate = 0.0;  /* Critical */

    /* Resilience order parameter:
     *   R = -λ (positive = resilient, zero = critical) */
    eco->resilience_order = -eco->recovery_rate;
    eco->is_resilient = (eco->recovery_rate < -0.05);
    eco->is_tipping = (eco->autocorr_lag1 > 0.8) && (eco->variance > 2.0 * 1.0);

    free(detrended);
    return eco;
}

void oep_app_ecosystem_free(OEPEcosystemResilience* eco) { free(eco); }

/* ── Application 5: Market Regime Detection ─────────────────── */

typedef struct {
    double regime_order;       /* OP: market regime (trending = high, mean-reverting = low) */
    double hurst_exponent;     /* H: H > 0.5 = trending, H < 0.5 = mean-reverting */
    double volatility;         /* σ: volatility */
    int    regime;             /* 0 = calm, 1 = volatile, 2 = crash */
    bool   is_critical;        /* Near a regime shift */
} OEPMarketRegime;

OEPMarketRegime* oep_app_market_regime_analyze(const double* returns, int n_obs) {
    /* Analyze financial market regime as an order parameter phenomenon.
     *
     * Market regimes can be characterized by:
     *   1. Hurst exponent H (persistence / mean reversion)
     *   2. Volatility clustering (GARCH-like dynamics)
     *   3. Correlation structure breakdown near crashes
     *
     * The order parameter is a composite indicator of "market coherence":
     * high values → trending (herding behavior synchronizes traders)
     * low values → mean-reverting (diverse opinions, efficient market)
     *
     * References:
     *   Sornette (2003) — Why Stock Markets Crash
     *   Preis et al. (2011) — Phys. Rev. E 84:056106 */
    if (!returns || n_obs < 100) return NULL;

    OEPMarketRegime* mr = calloc(1, sizeof(OEPMarketRegime));
    if (!mr) return NULL;

    /* Volatility */
    mr->volatility = sqrt(oep_ts_variance(returns, n_obs));

    /* Hurst exponent via R/S analysis (simplified) */
    int max_lag = n_obs / 4;
    double* rs_values = malloc((size_t)(max_lag - 1) * sizeof(double));
    double* lags = malloc((size_t)(max_lag - 1) * sizeof(double));
    if (rs_values && lags) {
        for (int lag = 10; lag < max_lag; lag++) {
            int n_segments = n_obs / lag;
            double rs_avg = 0.0;
            for (int seg = 0; seg < n_segments; seg++) {
                int start = seg * lag;
                double mean_seg = 0.0;
                for (int i = 0; i < lag; i++) mean_seg += returns[start + i];
                mean_seg /= (double)lag;

                double* cum_dev = malloc((size_t)lag * sizeof(double));
                if (!cum_dev) continue;
                double cum = 0.0;
                double R = -1e100, r_min = 1e100;
                for (int i = 0; i < lag; i++) {
                    cum += returns[start + i] - mean_seg;
                    cum_dev[i] = cum;
                    if (cum > R) R = cum;
                    if (cum < r_min) r_min = cum;
                }

                double S = 0.0;
                for (int i = 0; i < lag; i++)
                    S += (returns[start + i] - mean_seg) * (returns[start + i] - mean_seg);
                S = sqrt(S / (double)lag);

                if (S > OEP_EPSILON) rs_avg += (R - r_min) / S;
                free(cum_dev);
            }
            rs_avg /= (double)n_segments;
            rs_values[lag - 10] = log(rs_avg);
            lags[lag - 10] = log((double)lag);
        }

        /* Linear regression: log(R/S) = H * log(lag) + const */
        double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
        int n_pts = max_lag - 10;
        for (int i = 0; i < n_pts; i++) {
            sx += lags[i]; sy += rs_values[i];
            sxx += lags[i] * lags[i];
            sxy += lags[i] * rs_values[i];
        }
        double denom = (double)n_pts * sxx - sx * sx;
        if (fabs(denom) > OEP_EPSILON) {
            mr->hurst_exponent = ((double)n_pts * sxy - sx * sy) / denom;
        } else {
            mr->hurst_exponent = 0.5;
        }
    }
    free(rs_values); free(lags);

    /* Regime classification */
    if (mr->hurst_exponent > 0.6) {
        mr->regime = 1;  /* Trending (herding/persistent) */
    } else if (mr->hurst_exponent < 0.4) {
        mr->regime = 2;  /* Mean-reverting */
    } else {
        mr->regime = 0;  /* Efficient (random walk) */
    }

    mr->regime_order = mr->hurst_exponent;
    mr->is_critical = (mr->hurst_exponent > 0.7);  /* Approaching criticality */

    return mr;
}

void oep_app_market_free(OEPMarketRegime* mr) { free(mr); }

/* ── Application 6: Traffic Flow Order Parameter ────────────── */

typedef struct {
    double density;            /* ρ: vehicle density */
    double flow;               /* J = ρ·v: traffic flow */
    double order_param;        /* Free flow = low, congested = high */
    double critical_density;   /* ρ_c: critical density for jamming transition */
    bool   is_congested;
    bool   is_jammed;
} OEPTrafficFlow;

OEPTrafficFlow* oep_app_traffic_flow_analyze(const double* speeds,
                                              const double* densities,
                                              int n_obs) {
    /* Analyze traffic flow as a nonequilibrium phase transition.
     *
     * The traffic fundamental diagram J(ρ) exhibits:
     *   ρ < ρ_c: free flow (vehicles move freely, low correlation)
     *   ρ > ρ_c: congested flow (synchronized, stop-and-go waves)
     *   ρ ≫ ρ_c: jammed state (flow = 0 at jam density ρ_max)
     *
     * This is analogous to a second-order phase transition:
     *   Order parameter = velocity deviation from free-flow speed
     *   Control parameter = density
     *   ρ_c = critical density
     *
     * References:
     *   Kerner (2004) — The Physics of Traffic
     *   Nagel & Schreckenberg (1992) — J. Phys. I France 2:2221 */
    if (!speeds || !densities || n_obs < 10) return NULL;

    OEPTrafficFlow* tf = calloc(1, sizeof(OEPTrafficFlow));
    if (!tf) return NULL;

    /* Fit fundamental diagram: J = ρ·v */
    /* Find free-flow speed v_free from low-density data */
    double v_free = 0.0;
    int n_free = 0;
    for (int i = 0; i < n_obs; i++) {
        if (densities[i] < 10.0) {  /* Low density */
            v_free += speeds[i];
            n_free++;
        }
    }
    if (n_free > 0) v_free /= (double)n_free;
    if (v_free < OEP_EPSILON) v_free = 100.0;  /* Default (km/h) */

    /* Find critical density where speed starts dropping */
    double max_deriv = 0.0;
    double rho_c_est = 20.0;
    for (int i = 1; i < n_obs; i++) {
        double drho = densities[i] - densities[i-1];
        if (fabs(drho) < OEP_EPSILON) continue;
        double dv = fabs(speeds[i] - speeds[i-1]) / fabs(drho);
        if (dv > max_deriv) {
            max_deriv = dv;
            rho_c_est = densities[i];
        }
    }
    tf->critical_density = rho_c_est;

    /* Current state */
    tf->density = densities[n_obs - 1];
    tf->flow = tf->density * speeds[n_obs - 1];
    tf->order_param = 1.0 - speeds[n_obs - 1] / v_free;  /* 0 = free, 1 = jammed */
    if (tf->order_param < 0.0) tf->order_param = 0.0;
    if (tf->order_param > 1.0) tf->order_param = 1.0;

    tf->is_congested = (tf->density > tf->critical_density);
    tf->is_jammed = (speeds[n_obs - 1] < 5.0);  /* Almost stopped */

    return tf;
}

void oep_app_traffic_free(OEPTrafficFlow* tf) { free(tf); }

/* ── Application 7: Pattern Formation in BZ Reaction ────────── */

typedef struct {
    OEPSpatialField* concentration;  /* Chemical concentration field */
    double  order_param;             /* Amplitude of pattern (0 = uniform, >0 = patterned) */
    double  control_param;           /* Distance from Hopf bifurcation */
    bool    is_patterned;
    double  dominant_wavelength;     /* λ: characteristic wavelength of Turing pattern */
    int     nx, ny;
} OEPBZPattern;

OEPBZPattern* oep_app_bz_pattern_analyze(const double* conc_data,
                                          int nx, int ny) {
    /* Analyze Belousov-Zhabotinsky reaction pattern formation.
     *
     * The BZ reaction exhibits Turing patterns when the system is
     * driven far from equilibrium. The order parameter is the
     * amplitude of the concentration modulation.
     *
     * References:
     *   Turing (1952) — Phil. Trans. R. Soc. B 237:37
     *   Zaikin & Zhabotinsky (1970) — Nature 225:535 */
    if (!conc_data || nx < 10 || ny < 10) return NULL;

    OEPBZPattern* bz = calloc(1, sizeof(OEPBZPattern));
    if (!bz) return NULL;
    bz->nx = nx; bz->ny = ny;
    bz->concentration = oep_field_create(nx, ny, 1.0, 1.0);
    if (!bz->concentration) { free(bz); return NULL; }
    memcpy(bz->concentration->data, conc_data, (size_t)(nx * ny) * sizeof(double));

    /* Order parameter: spatial variance (Turing amplitude) */
    double mean = oep_field_spatial_mean(bz->concentration);
    double var = oep_field_spatial_variance(bz->concentration);
    bz->order_param = var;

    /* Determine if patterned: check structure factor for peaks */
    double* S_k = malloc((size_t)(nx / 2) * sizeof(double));
    if (S_k) {
        oep_structure_factor(bz->concentration->data, nx, ny, nx / 2, S_k);

        /* Find dominant wavelength from peak of S(k) */
        double max_S = 0.0;
        int max_k = 0;
        for (int k = 0; k < nx / 2; k++) {
            if (S_k[k] > max_S) { max_S = S_k[k]; max_k = k; }
        }
        if (max_k > 0) {
            bz->dominant_wavelength = (double)nx / (double)max_k;
            bz->is_patterned = (max_S > 2.0 * S_k[0]);
        } else {
            bz->dominant_wavelength = 0.0;
            bz->is_patterned = false;
        }
        free(S_k);
    }

    /* Control parameter: distance from uniform state threshold */
    bz->control_param = var;
    return bz;
}

void oep_app_bz_free(OEPBZPattern* bz) {
    if (!bz) return;
    oep_field_free(bz->concentration);
    free(bz);
}
