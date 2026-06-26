/**
 * @file fluctuation_theorems.h
 * @brief Fluctuation theorems in non-equilibrium statistical mechanics.
 *
 * Fluctuation theorems quantify the probability of observing
 * violations of the Second Law of Thermodynamics on small scales
 * and over short time intervals. They provide a deep connection
 * between equilibrium fluctuations and non-equilibrium response.
 *
 * Key theorems:
 *   - Evans-Searles Fluctuation Theorem (1994)
 *   - Crooks Fluctuation Theorem (1999)
 *   - Jarzynski Equality (1997)
 *   - Gallavotti-Cohen Fluctuation Theorem (1995)
 *   - Fluctuation-Dissipation Theorem (Callen-Welton, 1951)
 *
 * For a system driven out of equilibrium over time tau, the
 * probability of observing entropy production Sigma vs -Sigma:
 *
 *   P(Sigma) / P(-Sigma) = exp(Sigma / k_B)
 *
 * This means the Second Law (Sigma > 0) is overwhelmingly probable
 * for macroscopic systems, but violations are possible for small
 * systems observed over short times.
 *
 * References:
 *   - Evans, D.J. & Searles, D.J. (1994) Phys. Rev. E 50, 1645
 *   - Jarzynski, C. (1997) Phys. Rev. Lett. 78, 2690
 *   - Crooks, G.E. (1999) Phys. Rev. E 60, 2721
 *   - Gallavotti, G. & Cohen, E.G.D. (1995) Phys. Rev. Lett. 74, 2694
 *   - Seifert, U. (2012) Rep. Prog. Phys. 75, 126001
 */

#ifndef FLUCTUATION_THEOREMS_H
#define FLUCTUATION_THEOREMS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FT_MAX_SAMPLES 1000000
#define FT_MAX_BINS 1024

/* ---------------------------------------------------------------------------
 * L1: Definitions -- Fluctuation theorem structures
 * ------------------------------------------------------------------------- */

/**
 * @brief Stochastic trajectory of a small thermodynamic system.
 *
 * Records the time evolution of a system subject to thermal
 * fluctuations. For each time step, we track the state and
 * the incremental entropy production.
 */
typedef struct {
    double *time;              /**< Time points [s] (length n_steps) */
    double *position;          /**< Position x(t) [m] (for a Brownian particle) */
    double *velocity;          /**< Velocity v(t) [m/s] */
    double *entropy_production; /**< Cumulative entropy production Sigma(t) [J/K] */
    double *work_done;         /**< Cumulative work done on system W(t) [J] */
    double *heat_dissipated;   /**< Cumulative heat Q(t) [J] */
    int     n_steps;           /**< Number of time steps */
    double  dt;                /**< Time step [s] */
    double  temperature;       /**< Bath temperature [K] */
    double  k_B;               /**< Boltzmann constant [J/K] */
} StochasticTrajectory;

/**
 * @brief Histogram of entropy production for fluctuation theorem verification.
 *
 * Collects statistics of entropy production from many trajectories
 * to verify P(Sigma)/P(-Sigma) = exp(Sigma/k_B).
 */
typedef struct {
    double *bins;              /**< Bin centers for Sigma values */
    double *pdf;               /**< Probability density for each bin */
    double *pdf_reverse;       /**< PDF for time-reversed trajectories */
    int     num_bins;          /**< Number of histogram bins */
    double  bin_width;         /**< Width of each bin [J/K] */
    double  sigma_min;         /**< Minimum Sigma in histogram */
    double  sigma_max;         /**< Maximum Sigma in histogram */
    int     total_samples;     /**< Total number of trajectories sampled */
    double  k_B;               /**< Boltzmann constant [J/K] */
    double  asymmetry_ratio;   /**< Measured log[P(Sigma)/P(-Sigma)] */
    double  ft_slope;          /**< Slope of asymmetry vs Sigma (should = 1/k_B) */
    double  ft_r2;             /**< R^2 of the fluctuation theorem fit */
} EntropyProductionHistogram;

/**
 * @brief Work distribution for Jarzynski and Crooks relations.
 *
 * The Jarzynski equality relates the free energy difference Delta F
 * to the exponential average of non-equilibrium work:
 *
 *   <exp(-W/k_B T)> = exp(-Delta F / k_B T)
 *
 * where the average is over many realizations of a process that
 * drives the system between two equilibrium states.
 */
typedef struct {
    double *work_values;       /**< Work measured in each realization [J] */
    double *forward_work;      /**< Work in forward direction [J] */
    double *reverse_work;      /**< Work in reverse direction [J] */
    int     num_realizations;  /**< Number of experimental realizations */
    double  temperature;       /**< Bath temperature [K] */
    double  k_B;               /**< Boltzmann constant [J/K] */
    double  delta_F_true;      /**< True free energy difference [J] (if known) */
    double  delta_F_jarzynski; /**< Delta F estimated from Jarzynski equality [J] */
    double  delta_F_crooks;    /**< Delta F estimated from Crooks relation [J] */
    int     crooks_verified;   /**< 1 if Crooks fluctuation theorem holds */
} WorkDistribution;

/* ---------------------------------------------------------------------------
 * L2: Core Concepts -- Fluctuation theorem verification
 * ------------------------------------------------------------------------- */

/**
 * @brief Simulate a Brownian particle in a harmonic trap.
 *
 * Langevin equation:
 *   m dv/dt = -gamma v - k x + sqrt(2 gamma k_B T) xi(t)
 *
 * where xi(t) is Gaussian white noise with <xi(t)> = 0 and
 * <xi(t) xi(t')> = delta(t - t').
 *
 * @param traj        [out] Stochastic trajectory
 * @param mass        Particle mass m [kg]
 * @param friction    Friction coefficient gamma [kg/s]
 * @param spring_k    Trap stiffness k [N/m]
 * @param temperature Bath temperature [K]
 * @param k_B         Boltzmann constant [J/K]
 * @param x0          Initial position [m]
 * @param v0          Initial velocity [m/s]
 * @param dt          Time step [s]
 * @param n_steps     Number of steps
 * @return            0 on success
 *
 * Reference: Sekimoto (2010) "Stochastic Energetics"
 */
int ft_simulate_brownian_harmonic(StochasticTrajectory *traj,
                                   double mass, double friction,
                                   double spring_k, double temperature,
                                   double k_B, double x0, double v0,
                                   double dt, int n_steps);

/**
 * @brief Simulate a driven colloidal particle in a time-dependent trap.
 *
 * For Jarzynski/Crooks verification: the trap center is moved from
 * position 0 to position L over time tau:
 *
 *   U(x, t) = (1/2) k (x - lambda(t))^2
 *
 * where lambda(t) = (t/tau) L is the protocol.
 *
 * The work done on the particle is:
 *   W = integral_0^tau (partial U / partial lambda) dlambda
 *
 * @param traj        [out] Trajectory
 * @param mass        Mass [kg]
 * @param friction    Friction [kg/s]
 * @param spring_k    Trap stiffness [N/m]
 * @param L           Total displacement of trap center [m]
 * @param tau         Protocol duration [s]
 * @param temperature Bath temperature [K]
 * @param k_B         Boltzmann constant [J/K]
 * @param dt          Time step [s]
 * @return            0 on success
 */
int ft_simulate_jarzynski_protocol(StochasticTrajectory *traj,
                                    double mass, double friction,
                                    double spring_k, double L, double tau,
                                    double temperature, double k_B, double dt);

/**
 * @brief Compute the entropy production along a stochastic trajectory.
 *
 * For a Brownian particle, the stochastic entropy production is:
 *   Sigma(t) = Q(t) / T + Delta S_sys(t)
 *
 * where Q(t) is the heat dissipated and Delta S_sys is the change
 * in system entropy. The entropy production satisfies the
 * integral fluctuation theorem:
 *   <exp(-Sigma/k_B)> = 1
 *
 * @param traj   Trajectory [in/out: entropy_production[] filled]
 * @return       0 on success
 */
int ft_compute_trajectory_entropy(StochasticTrajectory *traj);

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures -- Fluctuation relations
 * ------------------------------------------------------------------------- */

/**
 * @brief Build an entropy production histogram from trajectory data.
 *
 * Collects entropy production values Sigma from many trajectories
 * into a histogram for verification of the fluctuation theorem.
 *
 * @param hist          [out] Histogram
 * @param sigma_values  Entropy production values from trajectories
 * @param n_samples     Number of samples
 * @param num_bins      Number of histogram bins
 * @param k_B           Boltzmann constant [J/K]
 * @return              0 on success
 */
int ft_build_sigma_histogram(EntropyProductionHistogram *hist,
                              const double sigma_values[],
                              int n_samples, int num_bins, double k_B);

/**
 * @brief Verify the Evans-Searles fluctuation theorem.
 *
 * For each bin i, compute:
 *   A(Sigma_i) = ln[P(Sigma_i) / P(-Sigma_i)]
 *
 * The theorem predicts A(Sigma) = Sigma / k_B.
 * Performs linear regression A ~ Sigma to find the slope
 * (which should equal 1/k_B) and R^2 goodness of fit.
 *
 * @param hist   Entropy production histogram [in/out: ft_slope, ft_r2 set]
 * @return       0 on success
 *
 * Reference: Evans & Searles (1994)
 */
int ft_verify_evans_searles(EntropyProductionHistogram *hist);

/**
 * @brief Verify the Crooks fluctuation theorem.
 *
 * For a process with time-reversed counterpart:
 *   P_F(W) / P_R(-W) = exp((W - Delta F) / k_B T)
 *
 * where P_F is the forward work distribution and P_R is the reverse.
 * Computes Delta F from the crossing point where P_F(W) = P_R(-W).
 *
 * @param wd   Work distribution [in/out: delta_F_crooks set]
 * @return     0 on success
 *
 * Reference: Crooks (1999)
 */
int ft_verify_crooks(WorkDistribution *wd);

/**
 * @brief Compute the Jarzynski free energy estimate.
 *
 * Delta F = -k_B T ln(<exp(-W/k_B T)>)
 *
 * The exponential average is over N forward realizations.
 * This function also computes the bias and variance of the
 * estimator using bootstrap resampling.
 *
 * @param wd   Work distribution [in/out: delta_F_jarzynski set]
 * @return     0 on success
 *
 * Reference: Jarzynski (1997)
 */
int ft_compute_jarzynski_deltaF(WorkDistribution *wd);

/* ---------------------------------------------------------------------------
 * L4: Fundamental Laws -- Fluctuation-Dissipation theorem
 * ------------------------------------------------------------------------- */

/**
 * @brief Verify the classical Fluctuation-Dissipation Theorem (FDT).
 *
 * The FDT relates the linear response function chi(t) to the
 * equilibrium correlation function C(t):
 *
 *   chi(t) = -(1/k_B T) dC(t)/dt  for t > 0
 *
 * where C(t) = <x(0) x(t)>_eq and chi(t) defines the response
 * <x(t)> = integral chi(t-s) f(s) ds to a weak perturbation f.
 *
 * @param traj   Equilibrium trajectory (no perturbation)
 * @param chi_t  [out] Response function values
 * @param C_t    [out] Correlation function values
 * @param tau    Maximum lag time [s]
 * @param n_lags Number of lag times
 * @return       Maximum relative deviation |chi(t) + dC/dt/kT| / |chi(t)|
 *
 * Reference: Callen & Welton (1951), Kubo (1966)
 */
double ft_verify_fdt(const StochasticTrajectory *traj,
                     double chi_t[], double C_t[],
                     double tau, int n_lags);

/* ---------------------------------------------------------------------------
 * L5: Algorithms -- Stochastic integration
 * ------------------------------------------------------------------------- */

/**
 * @brief Generate a Gaussian random number (Box-Muller method).
 *
 * @param mean   Mean of the distribution
 * @param stddev Standard deviation
 * @return       Gaussian random number
 */
double ft_rand_gaussian(double mean, double stddev);

/**
 * @brief Integrate a Langevin equation using the Euler-Maruyama scheme.
 *
 * dx = a(x, t) dt + b(x, t) dW
 *
 * where dW ~ sqrt(dt) * N(0, 1) is the Wiener increment.
 *
 * @param x       [in/out] State variable
 * @param drift   Drift function a(x, t)
 * @param noise   Noise amplitude b(x, t)
 * @param t       Current time [s]
 * @param dt      Time step [s]
 * @param params  Parameters for drift and noise functions
 */
void ft_euler_maruyama_step(double *x, double (*drift)(double, double, const void*),
                            double (*noise)(double, double, const void*),
                            double t, double dt, const void *params);

/**
 * @brief Bootstrap resample a dataset with replacement.
 *
 * Used for estimating confidence intervals on Jarzynski/Crooks
 * free energy estimates.
 *
 * @param data          Original data (length n)
 * @param n             Number of data points
 * @param resampled     [out] Resampled data (length n)
 * @param seed          Random seed
 */
void ft_bootstrap_resample(const double data[], int n,
                           double resampled[], unsigned int seed);

/**
 * @brief Run multiple independent trajectories to build work distribution.
 *
 * For Jarzynski/Crooks verification, runs N realizations of the
 * same protocol with different noise seeds.
 *
 * @param wd               [out] Work distribution
 * @param mass             Particle mass [kg]
 * @param friction         Friction coefficient [kg/s]
 * @param spring_k         Trap stiffness [N/m]
 * @param L                Protocol displacement [m]
 * @param tau              Protocol duration [s]
 * @param temperature      Bath temperature [K]
 * @param k_B              Boltzmann constant [J/K]
 * @param dt               Time step [s]
 * @param num_realizations Number of trajectories
 * @return                 0 on success
 */
int ft_build_work_distribution(WorkDistribution *wd,
                                double mass, double friction,
                                double spring_k, double L, double tau,
                                double temperature, double k_B,
                                double dt, int num_realizations);

/* Utility */

void ft_init_trajectory(StochasticTrajectory *traj, int n_steps,
                        double dt, double temperature, double k_B);
void ft_free_trajectory(StochasticTrajectory *traj);
void ft_init_histogram(EntropyProductionHistogram *hist, int num_bins,
                       double sigma_min, double sigma_max, double k_B);
void ft_free_histogram(EntropyProductionHistogram *hist);
void ft_init_work_distribution(WorkDistribution *wd, int num_realizations,
                               double temperature, double k_B);
void ft_free_work_distribution(WorkDistribution *wd);

#ifdef __cplusplus
}
#endif

#endif /* FLUCTUATION_THEOREMS_H */