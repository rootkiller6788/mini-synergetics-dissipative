#include "oep_emergence_metrics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Emergence Quantification Metrics Implementation
 *
 * Quantifying emergence — the appearance of macroscopic patterns
 * irreducible to microscopic rules — is a central challenge in
 * complex systems theory.
 *
 * This module implements several complementary metrics:
 *   1. Information-theoretic: excess entropy, mutual information
 *   2. Integrated information Φ (IIT framework)
 *   3. Causal emergence (Erik Hoel's framework)
 *   4. Downward causation measures
 *
 * References:
 *   Crutchfield & Young (1989) — Computational mechanics
 *   Tononi (2004) — Integrated Information Theory
 *   Hoel, Albantakis, Tononi (2013) — Causal emergence
 *   Prokopenko, Boschetti, Ryan (2009) — Information-theoretic emergence
 * ============================================================== */

/* ── Information-Theoretic Basics ──────────────────────────── */

double oep_em_entropy(const double* prob, int n_states) {
    /* Shannon entropy H = -Σ p_i log(p_i)
     * Returns 0 for degenerate distributions.
     * Uses natural log (nats). Multiply by log2(e) to convert to bits. */
    if (!prob || n_states <= 0) return 0.0;
    double H = 0.0;
    for (int i = 0; i < n_states; i++) {
        if (prob[i] > OEP_EPSILON)
            H -= prob[i] * log(prob[i]);
    }
    return H;
}

double oep_em_joint_entropy(const OEPMatrix* joint_prob,
                             int n_states_x, int n_states_y) {
    /* Joint entropy H(X, Y) = -Σ Σ p(x,y) log p(x,y)
     * joint_prob must be n_states_x × n_states_y. */
    if (!joint_prob || n_states_x <= 0 || n_states_y <= 0) return 0.0;
    double H = 0.0;
    for (int i = 0; i < n_states_x; i++) {
        for (int j = 0; j < n_states_y; j++) {
            double p = oep_mat_get(joint_prob, i, j);
            if (p > OEP_EPSILON) H -= p * log(p);
        }
    }
    return H;
}

double oep_em_mutual_information(const OEPMatrix* joint_prob,
                                  int n_states_x, int n_states_y) {
    /* Mutual information I(X; Y) = H(X) + H(Y) - H(X, Y)
     * = Σ Σ p(x,y) log[p(x,y) / (p(x)·p(y))] */
    if (!joint_prob || n_states_x <= 0 || n_states_y <= 0) return 0.0;

    /* Compute marginals */
    double* px = calloc((size_t)n_states_x, sizeof(double));
    double* py = calloc((size_t)n_states_y, sizeof(double));
    if (!px || !py) { free(px); free(py); return 0.0; }

    for (int i = 0; i < n_states_x; i++) {
        for (int j = 0; j < n_states_y; j++) {
            double p = oep_mat_get(joint_prob, i, j);
            px[i] += p;
            py[j] += p;
        }
    }

    double I = 0.0;
    for (int i = 0; i < n_states_x; i++) {
        for (int j = 0; j < n_states_y; j++) {
            double pxy = oep_mat_get(joint_prob, i, j);
            if (pxy > OEP_EPSILON && px[i] > OEP_EPSILON && py[j] > OEP_EPSILON)
                I += pxy * log(pxy / (px[i] * py[j]));
        }
    }

    free(px);
    free(py);
    return I;
}

double oep_em_conditional_entropy(const OEPMatrix* joint_prob,
                                   int n_states_x, int n_states_y) {
    /* H(X|Y) = H(X, Y) - H(Y) */
    double H_xy = oep_em_joint_entropy(joint_prob, n_states_x, n_states_y);

    /* Compute marginal p(y) */
    double* py = calloc((size_t)n_states_y, sizeof(double));
    if (!py) return H_xy;
    for (int i = 0; i < n_states_x; i++)
        for (int j = 0; j < n_states_y; j++)
            py[j] += oep_mat_get(joint_prob, i, j);

    double H_y = oep_em_entropy(py, n_states_y);
    free(py);
    return H_xy - H_y;
}

double oep_em_transfer_entropy(const double* x, const double* y,
                                int n, int delay, int n_bins) {
    /* Transfer entropy TE_{Y→X} = I(X_{t+delay}; Y_t | X_t)
     *
     * Quantifies the directed information flow from Y to X,
     * controlling for the auto-predictability of X.
     *
     * Implementation: discretize into bins, compute conditional
     * mutual information from the joint histogram.
     *
     * Reference: Schreiber (2000) Phys. Rev. Lett. 85:461 */
    if (!x || !y || n < delay + 2 || n_bins < 2) return 0.0;

    /* Find ranges */
    double xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0];
    for (int i = 1; i < n; i++) {
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i];
        if (y[i] > ymax) ymax = y[i];
    }
    double xrng = xmax - xmin, yrng = ymax - ymin;
    if (xrng < OEP_EPSILON) xrng = 1.0;
    if (yrng < OEP_EPSILON) yrng = 1.0;

    /* Build 3D histogram: p(x_future, x_past, y_past) */
    int h_size = n_bins * n_bins * n_bins;
    double* hist = calloc((size_t)h_size, sizeof(double));
    if (!hist) return 0.0;

    int N = n - delay;
    for (int t = 0; t < N; t++) {
        int bx_past = (int)((x[t] - xmin) / xrng * (double)(n_bins - 1));
        int by_past = (int)((y[t] - ymin) / yrng * (double)(n_bins - 1));
        int bx_fut  = (int)((x[t + delay] - xmin) / xrng * (double)(n_bins - 1));
        if (bx_past < 0) bx_past = 0;
        if (bx_past >= n_bins) bx_past = n_bins - 1;
        if (by_past < 0) by_past = 0;
        if (by_past >= n_bins) by_past = n_bins - 1;
        if (bx_fut < 0) bx_fut = 0;
        if (bx_fut >= n_bins) bx_fut = n_bins - 1;
        hist[(bx_fut * n_bins + bx_past) * n_bins + by_past] += 1.0;
    }
    for (int i = 0; i < h_size; i++) hist[i] /= (double)N;

    /* H(X_future | X_past) */
    double* p_fut_past = calloc((size_t)(n_bins * n_bins), sizeof(double));
    double* p_past = calloc((size_t)n_bins, sizeof(double));
    if (!p_fut_past || !p_past) { free(hist); free(p_fut_past); free(p_past); return 0.0; }

    for (int k = 0; k < n_bins; k++) {
        for (int i = 0; i < n_bins; i++) {
            for (int j = 0; j < n_bins; j++) {
                double p = hist[(k * n_bins + i) * n_bins + j];
                p_fut_past[k * n_bins + i] += p;
                p_past[i] += p;
            }
        }
    }

    double H_fut_given_past = 0.0;
    for (int i = 0; i < n_bins; i++) {
        for (int k = 0; k < n_bins; k++) {
            double p_joint = p_fut_past[k * n_bins + i];
            if (p_joint > OEP_EPSILON && p_past[i] > OEP_EPSILON)
                H_fut_given_past -= p_joint * log(p_joint / p_past[i]);
        }
    }

    /* H(X_future | X_past, Y_past) */
    double* p_fut_past_ypast = calloc((size_t)h_size, sizeof(double));
    double* p_xy = calloc((size_t)(n_bins * n_bins), sizeof(double));
    if (!p_fut_past_ypast || !p_xy) {
        free(hist); free(p_fut_past); free(p_past);
        free(p_fut_past_ypast); free(p_xy);
        return 0.0;
    }
    memcpy(p_fut_past_ypast, hist, (size_t)h_size * sizeof(double));
    for (int i = 0; i < n_bins; i++) {
        for (int j = 0; j < n_bins; j++) {
            for (int k = 0; k < n_bins; k++) {
                p_xy[i * n_bins + j] += hist[(k * n_bins + i) * n_bins + j];
            }
        }
    }

    double H_fut_given_past_y = 0.0;
    for (int i = 0; i < n_bins; i++) {
        for (int j = 0; j < n_bins; j++) {
            for (int k = 0; k < n_bins; k++) {
                double p_joint3 = hist[(k * n_bins + i) * n_bins + j];
                double p_cond = p_xy[i * n_bins + j];
                if (p_joint3 > OEP_EPSILON && p_cond > OEP_EPSILON)
                    H_fut_given_past_y -= p_joint3 * log(p_joint3 / p_cond);
            }
        }
    }

    free(hist); free(p_fut_past); free(p_past);
    free(p_fut_past_ypast); free(p_xy);

    return H_fut_given_past - H_fut_given_past_y;
}

/* ── Predictive Information ────────────────────────────────── */

double oep_em_predictive_information(const double* data, int n,
                                      int past_len, int future_len) {
    /* Predictive information:
     *   I_pred = I(X_past; X_future)
     *          = H(X_future) - H(X_future | X_past)
     *
     * For a Markov process of order k, I_pred plateaus after past_len > k.
     * High I_pred indicates structured, non-random dynamics.
     *
     * This is estimated using block entropy from the time series. */
    if (!data || n < past_len + future_len + 2) return 0.0;

    /* Simple approach: use Gaussian approximation for real-valued data.
     * For Gaussian processes:
     *   I_pred = (1/2) log(det(Σ_future) / det(Σ_future|past))
     *
     * For computational efficiency, we use the correlation-based estimate:
     *   I_pred ≈ -(1/2) log(1 - ρ²)
     * where ρ is the maximum correlation between past and future windows. */
    int N = n - past_len - future_len;
    if (N < 1) return 0.0;

    double max_corr = 0.0;
    for (int start = 0; start < N; start++) {
        double* past = (double*)data + start;
        double* future = (double*)data + start + past_len;
        /* Compute correlation between past and future means */
        double mu_p = 0.0, mu_f = 0.0;
        for (int i = 0; i < past_len; i++) mu_p += past[i];
        mu_p /= (double)past_len;
        for (int i = 0; i < future_len; i++) mu_f += future[i];
        mu_f /= (double)future_len;

        double cov = 0.0, var_p = 0.0, var_f = 0.0;
        int min_len = (past_len < future_len) ? past_len : future_len;
        for (int i = 0; i < min_len; i++) {
            double dp = past[i] - mu_p;
            double df = future[i] - mu_f;
            cov += dp * df;
            var_p += dp * dp;
            var_f += df * df;
        }
        if (var_p > OEP_EPSILON && var_f > OEP_EPSILON) {
            double corr = fabs(cov) / sqrt(var_p * var_f);
            if (corr > max_corr) max_corr = corr;
        }
    }

    if (max_corr >= 1.0 - OEP_EPSILON) return 1e10;  /* Perfect prediction */
    if (max_corr <= OEP_EPSILON) return 0.0;
    return -0.5 * log(1.0 - max_corr * max_corr);
}

double oep_em_excess_entropy(const double* data, int n,
                              int max_block_len, double entropy_rate) {
    /* Excess entropy E = lim_{L→∞} [H(L) - L·h_μ]
     *
     * Where H(L) is the block entropy of length L and h_μ is the
     * entropy rate (per symbol) of the process.
     *
     * E measures the "stored" information — how much the past informs
     * the future beyond the entropy rate.
     *
     * For a periodic process: E = log(period) (large)
     * For a random process: E = 0
     * For a chaotic process: E is finite and related to the complexity. */
    if (!data || n < max_block_len + 2) return 0.0;

    double max_H_L = 0.0;
    int max_L = 0;
    for (int L = 1; L <= max_block_len && L < n; L++) {
        /* Estimate block entropy H(L) from data using Gaussian approximation */
        int N_blocks = n - L + 1;
        double* block_means = malloc((size_t)N_blocks * sizeof(double));
        if (!block_means) break;
        for (int i = 0; i < N_blocks; i++) {
            double s = 0.0;
            for (int j = 0; j < L; j++) s += data[i + j];
            block_means[i] = s / (double)L;
        }
        double mu = oep_ts_mean(block_means, N_blocks);
        double var = oep_ts_variance(block_means, N_blocks);

        /* Gaussian entropy: H(L) ≈ 1/2 log(2πe σ²) */
        double H_L;
        if (var > OEP_EPSILON)
            H_L = 0.5 * log(2.0 * OEP_PI * OEP_E * var);
        else
            H_L = 0.0;

        double excess = H_L - (double)L * entropy_rate;
        if (excess > max_H_L) { max_H_L = excess; max_L = L; }
        free(block_means);
    }
    (void)max_L;
    return max_H_L;
}

void oep_em_complexity_entropy_pair(const double* data, int n,
                                     double* entropy_rate,
                                     double* excess_entropy) {
    /* Compute the (entropy_rate, excess_entropy) pair for the
     * Crutchfield-Young complexity-entropy diagram.
     *
     * Entropy rate (h_μ): measures randomness per time step.
     * Excess entropy (E): measures structural complexity.
     *
     * On the C-E plane:
     *   - Periodic processes: (h_μ = 0, E finite)
     *   - Random processes: (h_μ = H_max, E = 0)
     *   - Complex processes: (h_μ moderate, E finite)
     *
     * "Complexity" in this framework = E (excess entropy). */
    if (!data || !entropy_rate || !excess_entropy || n < 10) {
        if (entropy_rate) *entropy_rate = 0.0;
        if (excess_entropy) *excess_entropy = 0.0;
        return;
    }

    /* Estimate entropy rate from the variance of differences (Δx) */
    double* diff = malloc((size_t)(n - 1) * sizeof(double));
    if (!diff) { *entropy_rate = 0.0; *excess_entropy = 0.0; return; }
    for (int i = 0; i < n - 1; i++) diff[i] = data[i + 1] - data[i];
    double var_diff = oep_ts_variance(diff, n - 1);
    free(diff);

    /* Gaussian entropy rate: h_μ ≈ 1/2 log(2πe σ²_{Δx}) */
    if (var_diff > OEP_EPSILON)
        *entropy_rate = 0.5 * log(2.0 * OEP_PI * OEP_E * var_diff);
    else
        *entropy_rate = 0.0;

    *excess_entropy = oep_em_excess_entropy(data, n,
        (n > 20) ? 20 : n / 2, *entropy_rate);
}

/* ── Effective Information ─────────────────────────────────── */

double oep_em_effective_information(const OEPMatrix* transition_matrix,
                                     int n_states) {
    /* Effective Information EI = I(X_t; X_{t+1}) with uniform input.
     *
     * EI = Σ_x p_u(x) D_KL[ p(x_{t+1}|x_t=x) || p_u(x_{t+1}) ]
     *
     * where p_u is uniform distribution.
     *
     * EI measures the causal power: how much does knowing the current
     * state constrain the future state, under uniform intervention?
     *
     * High EI → deterministic dynamics
     * Low EI → noise-dominated dynamics */
    if (!transition_matrix || n_states <= 1) return 0.0;

    double u = 1.0 / (double)n_states;
    double EI = 0.0;

    /* Compute the "effect distribution" (average of rows) */
    double* avg_effect = calloc((size_t)n_states, sizeof(double));
    if (!avg_effect) return 0.0;
    for (int j = 0; j < n_states; j++) {
        for (int i = 0; i < n_states; i++)
            avg_effect[j] += oep_mat_get(transition_matrix, i, j);
        avg_effect[j] *= u;
    }

    /* EI = Σ_i u · Σ_j T_{ij} log(T_{ij} / avg_effect[j]) */
    for (int i = 0; i < n_states; i++) {
        for (int j = 0; j < n_states; j++) {
            double Tij = oep_mat_get(transition_matrix, i, j);
            if (Tij > OEP_EPSILON && avg_effect[j] > OEP_EPSILON)
                EI += u * Tij * log(Tij / avg_effect[j]);
        }
    }

    free(avg_effect);
    return EI;
}

double oep_em_integrated_information_bipartition(
    const OEPMatrix* full_transition,
    int n1_states, int n2_states) {
    /* Integrated information Φ for bipartition into parts of sizes n1, n2.
     *
     * Φ = EI_whole - EI_part1 - EI_part2
     *
     * where EI_whole is the effective information of the full system,
     * and EI_parti is the EI of each part when disconnected.
     *
     * This implements the IIT 2.0 bipartition measure.
     * Reference: Balduzzi & Tononi (2008) PLoS Comput. Biol. */
    if (!full_transition || n1_states <= 0 || n2_states <= 0) return 0.0;
    int n_full = n1_states * n2_states;
    if (full_transition->rows != n_full || full_transition->cols != n_full) return 0.0;

    /* EI of full system */
    double EI_whole = oep_em_effective_information(full_transition, n_full);

    /* EI of part 1 (marginalize over part 2 states with uniform input to part 1) */
    OEPMatrix* T1 = oep_mat_create(n1_states, n1_states);
    if (!T1) return 0.0;
    for (int i1 = 0; i1 < n1_states; i1++) {
        for (int j1 = 0; j1 < n1_states; j1++) {
            double sum = 0.0;
            for (int i2 = 0; i2 < n2_states; i2++) {
                for (int j2 = 0; j2 < n2_states; j2++) {
                    int from_idx = i1 * n2_states + i2;
                    int to_idx = j1 * n2_states + j2;
                    sum += oep_mat_get(full_transition, from_idx, to_idx);
                }
            }
            /* Marginalize: average over part 2 states */
            oep_mat_set(T1, i1, j1, sum / (double)n2_states);
        }
    }
    double EI_1 = oep_em_effective_information(T1, n1_states);
    oep_mat_free(T1);

    /* EI of part 2 */
    OEPMatrix* T2 = oep_mat_create(n2_states, n2_states);
    if (!T2) { return EI_whole - EI_1; }
    for (int i2 = 0; i2 < n2_states; i2++) {
        for (int j2 = 0; j2 < n2_states; j2++) {
            double sum = 0.0;
            for (int i1 = 0; i1 < n1_states; i1++) {
                for (int j1 = 0; j1 < n1_states; j1++) {
                    int from_idx = i1 * n2_states + i2;
                    int to_idx = j1 * n2_states + j2;
                    sum += oep_mat_get(full_transition, from_idx, to_idx);
                }
            }
            oep_mat_set(T2, i2, j2, sum / (double)n1_states);
        }
    }
    double EI_2 = oep_em_effective_information(T2, n2_states);
    oep_mat_free(T2);

    return EI_whole - EI_1 - EI_2;
}

double oep_em_min_information_bipartition(
    const OEPMatrix* full_transition,
    int n_vars, int states_per_var) {
    /* Minimum Information Bipartition (MIB): find the cut that
     * minimizes the integrated information Φ.
     *
     * For systems with n_vars binary variables (states_per_var = 2):
     * There are 2^{n_vars-1} possible bipartitions.
     *
     * We use a greedy heuristic: start with all variables in part 1,
     * iteratively move variables to part 2 to minimize Φ.
     *
     * Returns the minimum Φ value. */
    if (!full_transition || n_vars <= 1 || states_per_var <= 1) return 0.0;

    /* For simplicity, approximate: try all single-variable-part-2 cuts */
    double min_phi = 1e100;
    int n_full = 1;
    for (int i = 0; i < n_vars; i++) n_full *= states_per_var;

    for (int var = 0; var < n_vars; var++) {
        int n1 = 1;
        for (int i = 0; i < var; i++) n1 *= states_per_var;
        int n2 = n_full / n1;

        double phi = oep_em_integrated_information_bipartition(
            full_transition, n1, n2);
        if (phi < min_phi) min_phi = phi;
    }

    return (min_phi < 1e99) ? min_phi : 0.0;
}

/* ── Causal Emergence ──────────────────────────────────────── */

OEPMatrix* oep_em_coarse_grain_transition(const OEPMatrix* micro_trans,
                                           int n_micro_states,
                                           const int* coarse_mapping,
                                           int n_macro_states) {
    /* Coarse-grain a micro-state transition probability matrix.
     *
     * The mapping coarse_mapping[i] = j means micro state i maps to macro state j.
     *
     * The coarse-grained transition probability is computed using the
     * time-averaged intervention distribution:
     *
     *   P(M_{t+1}=y | M_t=x) ∝ Σ_{i∈x} Σ_{j∈y} P(micro_{t+1}=j | micro_t=i)
     *
     * with the micro state uniformly weighted within each macro state. */
    if (!micro_trans || !coarse_mapping || n_micro_states <= 0 || n_macro_states <= 0)
        return NULL;

    OEPMatrix* macro_trans = oep_mat_create(n_macro_states, n_macro_states);
    if (!macro_trans) return NULL;

    /* Count micro states per macro state */
    int* micro_count = calloc((size_t)n_macro_states, sizeof(int));
    if (!micro_count) { oep_mat_free(macro_trans); return NULL; }
    for (int i = 0; i < n_micro_states; i++) {
        int m = coarse_mapping[i];
        if (m >= 0 && m < n_macro_states) micro_count[m]++;
    }

    /* Compute coarse-grained transition probabilities */
    for (int x = 0; x < n_macro_states; x++) {
        if (micro_count[x] == 0) continue;
        for (int y = 0; y < n_macro_states; y++) {
            double sum = 0.0;
            for (int i = 0; i < n_micro_states; i++) {
                if (coarse_mapping[i] != x) continue;
                for (int j = 0; j < n_micro_states; j++) {
                    if (coarse_mapping[j] != y) continue;
                    sum += oep_mat_get(micro_trans, i, j);
                }
            }
            /* Average over micro states in source macro state */
            oep_mat_set(macro_trans, x, y, sum / (double)micro_count[x]);
        }
    }

    free(micro_count);
    return macro_trans;
}

double oep_em_causal_emergence(const OEPMatrix* micro_trans,
                                int n_micro_states,
                                const int* coarse_mapping,
                                int n_macro_states) {
    /* Compute causal emergence:
     *   CE = EI_macro - EI_micro
     *
     * If CE > 0: the macro description has MORE causal power.
     * If CE < 0: the micro description is more causally informative.
     *
     * CE > 0 is the signature of genuine emergence (Hoel et al., 2013). */
    if (!micro_trans || !coarse_mapping || n_micro_states <= 0 || n_macro_states <= 0)
        return 0.0;

    double EI_micro = oep_em_effective_information(micro_trans, n_micro_states);

    OEPMatrix* macro_trans = oep_em_coarse_grain_transition(
        micro_trans, n_micro_states, coarse_mapping, n_macro_states);
    if (!macro_trans) return 0.0;

    double EI_macro = oep_em_effective_information(macro_trans, n_macro_states);
    oep_mat_free(macro_trans);

    return EI_macro - EI_micro;
}

double oep_em_optimal_causal_emergence(const OEPMatrix* micro_trans,
                                        int n_micro_states,
                                        int max_macro_states,
                                        int* best_mapping) {
    /* Find the coarse-graining that maximizes CE.
     *
     * We search over partitions of micro states into ≤ max_macro_states
     * groups. This is NP-hard in general; we use a greedy merging
     * heuristic: start with identity mapping, iteratively merge the
     * two macro states that produce the largest CE increase.
     *
     * Reference: Hoel et al. (2013) — Causal emergence */
    if (!micro_trans || n_micro_states <= 1 || max_macro_states < 1) {
        if (best_mapping) for (int i = 0; i < n_micro_states; i++) best_mapping[i] = i;
        return 0.0;
    }

    /* Initialize with identity mapping */
    int n_macro = n_micro_states;
    int* mapping = malloc((size_t)n_micro_states * sizeof(int));
    if (!mapping) return 0.0;
    for (int i = 0; i < n_micro_states; i++) mapping[i] = i;

    double best_CE = oep_em_effective_information(micro_trans, n_micro_states);

    /* Greedy minimization of EI (we invert to maximize CE) */
    /* For each reduction step */
    while (n_macro > max_macro_states && n_macro > 1) {
        /* Find best pair to merge */
        int best_a = -1, best_b = -1;
        double best_ei = -1e100;

        for (int a = 0; a < n_micro_states; a++) {
            for (int b = a + 1; b < n_micro_states; b++) {
                if (mapping[a] == mapping[b]) continue;

                /* Try merging groups mapping[a] and mapping[b] */
                int* trial_mapping = malloc((size_t)n_micro_states * sizeof(int));
                if (!trial_mapping) continue;
                memcpy(trial_mapping, mapping, (size_t)n_micro_states * sizeof(int));

                int old_b = mapping[b];
                for (int i = 0; i < n_micro_states; i++)
                    if (trial_mapping[i] == old_b)
                        trial_mapping[i] = mapping[a];

                /* Compute macro EI */
                OEPMatrix* macro = oep_em_coarse_grain_transition(
                    micro_trans, n_micro_states, trial_mapping, n_macro - 1);
                if (macro) {
                    double ei = oep_em_effective_information(macro, n_macro - 1);
                    if (ei > best_ei) { best_ei = ei; best_a = a; best_b = b; }
                    oep_mat_free(macro);
                }
                free(trial_mapping);
            }
        }

        if (best_a < 0 || best_b < 0) break;

        /* Merge groups */
        int old_b = mapping[best_b];
        for (int i = 0; i < n_micro_states; i++)
            if (mapping[i] == old_b) mapping[i] = mapping[best_a];
        n_macro--;

        if (best_ei > best_CE) best_CE = best_ei;
    }

    /* Normalize mapping to 0..n_macro-1 */
    int next_id = 0;
    int* remap = calloc((size_t)n_micro_states, sizeof(int));
    for (int i = 0; i < n_micro_states; i++) remap[i] = -1;
    for (int i = 0; i < n_micro_states; i++) {
        if (remap[mapping[i]] < 0) remap[mapping[i]] = next_id++;
        mapping[i] = remap[mapping[i]];
    }
    free(remap);

    if (best_mapping) memcpy(best_mapping, mapping, (size_t)n_micro_states * sizeof(int));
    free(mapping);

    /* CE = EI(best_macro) - EI(micro) */
    return best_CE - oep_em_effective_information(micro_trans, n_micro_states);
}

/* ── Downward Causation ────────────────────────────────────── */

double oep_em_downward_causation(const OEPMatrix* micro_data,
                                  const OEPVector* macro_state,
                                  int n_obs, int n_micro_vars,
                                  int n_macro_states, int n_bins) {
    /* Downward causation: DC = H(micro_{t+1} | micro_t) - H(micro_{t+1} | micro_t, macro_t)
     *
     * This measures how much the macro (order parameter) state reduces
     * uncertainty about future micro states beyond what the micro
     * history alone tells us.
     *
     * DC > 0: macro state carries causal information not in micro states
     * (genuine downward causation). */
    if (!micro_data || !macro_state || n_obs < 3 || n_micro_vars < 1) return 0.0;

    /* Discretize the first micro variable and macro state */
    double* micro_var1 = malloc((size_t)n_obs * sizeof(double));
    if (!micro_var1) return 0.0;
    for (int t = 0; t < n_obs; t++)
        micro_var1[t] = oep_mat_get(micro_data, t, 0);

    /* Transfer entropy from macro to micro as a proxy for DC */
    double TE = oep_em_transfer_entropy(micro_var1, macro_state->data,
                                         n_obs, 1,
                                         (n_bins > 0) ? n_bins : 10);
    free(micro_var1);
    return TE;
    (void)n_macro_states;  /* Reserved for discretized version */
}

/* ── Multiscale Complexity ─────────────────────────────────── */

double oep_em_multiscale_complexity(const double* data, int n,
                                     int n_scales) {
    /* Multiscale complexity: sum of excess entropies across scales.
     *
     * For each scale s, we coarse-grain the data by taking means over
     * windows of size 2^s, then compute the excess entropy.
     *
     * Total complexity = Σ_s E_s
     *
     * This captures how much structure exists at each scale. */
    if (!data || n < 2 || n_scales <= 0) return 0.0;

    double total_E = 0.0;
    double* working = malloc((size_t)n * sizeof(double));
    if (!working) return 0.0;

    int current_n = n;
    memcpy(working, data, (size_t)n * sizeof(double));

    for (int s = 0; s < n_scales; s++) {
        if (current_n < 4) break;

        /* Compute excess entropy at this scale */
        double rate, E;
        oep_em_complexity_entropy_pair(working, current_n, &rate, &E);
        total_E += E;

        /* Coarse-grain: take means of adjacent pairs */
        int next_n = current_n / 2;
        for (int i = 0; i < next_n; i++)
            working[i] = (working[2 * i] + working[2 * i + 1]) * 0.5;
        current_n = next_n;
    }

    free(working);
    return total_E;
}

/* ── Holistic Emergence Analysis ───────────────────────────── */

OEPEmergenceResult* oep_em_analyze_emergence(const OEPMatrix* micro_data,
                                              const OEPVector* macro_data,
                                              int n_obs, int n_micro_vars,
                                              int n_bins_micro, int n_bins_macro) {
    if (!micro_data || !macro_data || n_obs < 10 || n_micro_vars < 1) return NULL;

    OEPEmergenceResult* r = calloc(1, sizeof(OEPEmergenceResult));
    if (!r) return NULL;

    /* Compute mutual information between macro and micro */
    double* micro_var1 = malloc((size_t)n_obs * sizeof(double));
    if (!micro_var1) { free(r); return NULL; }
    for (int t = 0; t < n_obs; t++)
        micro_var1[t] = oep_mat_get(micro_data, t, 0);

    /* Mutual information: I(macro; micro)
     * Approximate from correlation ρ: I ≈ -1/2 log(1-ρ²) (Gaussian) */
    double corr = 0.0;
    double m_micro = oep_ts_mean(micro_var1, n_obs);
    double m_macro = oep_ts_mean(macro_data->data, n_obs);
    double cov = 0.0, var_mic = 0.0, var_mac = 0.0;
    for (int t = 0; t < n_obs; t++) {
        double dm = micro_var1[t] - m_micro;
        double dM = macro_data->data[t] - m_macro;
        cov += dm * dM;
        var_mic += dm * dm;
        var_mac += dM * dM;
    }
    if (var_mic > OEP_EPSILON && var_mac > OEP_EPSILON)
        corr = cov / sqrt(var_mic * var_mac);
    if (fabs(corr) < 1.0 - OEP_EPSILON)
        r->mutual_info_macro_micro = -0.5 * log(1.0 - corr * corr);

    /* Compression ratio */
    r->compression_ratio = (double)n_micro_vars;  /* micro dim / macro dim = n_vars / 1 */

    /* Coherence: autocorrelation of the macro at lag 1 */
    r->coherence = oep_ts_autocorrelation(macro_data->data, n_obs, 1);

    /* Downward causation */
    r->downward_causation = oep_em_downward_causation(
        micro_data, macro_data, n_obs, n_micro_vars, n_bins_macro, n_bins_micro);

    /* Default values for metrics requiring full transition matrix */
    r->excess_entropy = 0.0;
    r->integrated_information = 0.0;
    r->causal_emergence = 0.0;
    r->novelty = 0.0;
    r->is_genuinely_emergent = false;

    free(micro_var1);
    return r;
}

void oep_em_result_free(OEPEmergenceResult* result) { free(result); }

bool oep_em_is_genuine_emergence(const OEPEmergenceResult* r,
                                   double ce_threshold,
                                   double phi_threshold) {
    if (!r) return false;
    return (r->causal_emergence > ce_threshold) ||
           (r->integrated_information > phi_threshold) ||
           (r->downward_causation > 0.1);
}

void oep_em_result_print(const OEPEmergenceResult* r) {
    if (!r) return;
    printf("=== Emergence Analysis ===\n");
    printf("Excess Entropy:           %12.6f\n", r->excess_entropy);
    printf("Integrated Information Φ: %12.6f\n", r->integrated_information);
    printf("Causal Emergence CE:      %12.6f\n", r->causal_emergence);
    printf("Downward Causation:       %12.6f\n", r->downward_causation);
    printf("Mutual Info (Macro/Micro):%12.6f\n", r->mutual_info_macro_micro);
    printf("Compression Ratio:        %12.6f\n", r->compression_ratio);
    printf("Coherence:                %12.6f\n", r->coherence);
    printf("Novelty:                  %12.6f\n", r->novelty);
    printf("Genuinely Emergent:       %s\n",
           r->is_genuinely_emergent ? "YES" : "NO");
}
