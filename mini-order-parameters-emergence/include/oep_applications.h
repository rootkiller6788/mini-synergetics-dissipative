#ifndef OEP_APPLICATIONS_H
#define OEP_APPLICATIONS_H
#include "oep_core.h"
#include "oep_order_parameter.h"
#include "oep_landau_ginzburg.h"

/* ── EN SO Analysis ────────────────────────────────────── */
typedef struct {
    double* sst_anomalies;
    int     n_obs;
    double  order_param;
    double  threshold;
    double  transition_param;
    double  critical_value;
    bool    is_el_nino;
    bool    is_la_nina;
} OEPENSOAnalysis;

OEPENSOAnalysis* oep_app_enso_analyze(const double* sst_data, int n_obs);
void             oep_app_enso_free(OEPENSOAnalysis* enso);

/* ── Brain Coherence ───────────────────────────────────── */
typedef struct {
    double* eeg_channels;
    int     n_channels;
    int     n_samples;
    double  kuramoto_order;
    double  phi_integrated;
    double  complexity;
    bool    is_conscious;
} OEPBrainCoherence;

OEPBrainCoherence* oep_app_brain_coherence_analyze(const double* eeg_data,
                                                    int n_channels, int n_samples);
void               oep_app_brain_free(OEPBrainCoherence* bc);

/* ── Opinion Polarization ──────────────────────────────── */
typedef struct {
    double polarization_order;
    double control_param;
    double critical_value;
    bool   is_polarized;
    int*   opinion_counts;
    int    n_bins;
} OEPOpinionPolarization;

OEPOpinionPolarization* oep_app_opinion_polarization_analyze(const double* opinions,
                                                               int n_people);
void                    oep_app_opinion_free(OEPOpinionPolarization* pol);

/* ── Ecosystem Resilience ──────────────────────────────── */
typedef struct {
    double resilience_order;
    double recovery_rate;
    double autocorr_lag1;
    double variance;
    bool   is_resilient;
    bool   is_tipping;
} OEPEcosystemResilience;

OEPEcosystemResilience* oep_app_ecosystem_resilience_analyze(const double* biomass,
                                                               int n_obs);
void                    oep_app_ecosystem_free(OEPEcosystemResilience* eco);

/* ── Market Regime ─────────────────────────────────────── */
typedef struct {
    double regime_order;
    double hurst_exponent;
    double volatility;
    int    regime;
    bool   is_critical;
} OEPMarketRegime;

OEPMarketRegime* oep_app_market_regime_analyze(const double* returns, int n_obs);
void             oep_app_market_free(OEPMarketRegime* mr);

/* ── Traffic Flow ──────────────────────────────────────── */
typedef struct {
    double density;
    double flow;
    double order_param;
    double critical_density;
    bool   is_congested;
    bool   is_jammed;
} OEPTrafficFlow;

OEPTrafficFlow* oep_app_traffic_flow_analyze(const double* speeds,
                                              const double* densities,
                                              int n_obs);
void            oep_app_traffic_free(OEPTrafficFlow* tf);

/* ── BZ Pattern Formation ──────────────────────────────── */
typedef struct {
    OEPSpatialField* concentration;
    double  order_param;
    double  control_param;
    bool    is_patterned;
    double  dominant_wavelength;
    int     nx, ny;
} OEPBZPattern;

OEPBZPattern* oep_app_bz_pattern_analyze(const double* conc_data,
                                          int nx, int ny);
void          oep_app_bz_free(OEPBZPattern* bz);

#endif
