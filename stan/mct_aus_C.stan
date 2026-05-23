// Australian Multivariate Core Trend (MCT) inflation — Variant C
// (Combined common-trend RW + common-cycle white noise, tighter SV priors).
//
// Conceptual analogue of NY Fed MCT's tau_c + eps_c structure — sectors load
// onto BOTH a slow-moving common trend AND a transient common cycle, rather
// than the model being forced to choose between A (cycle only) and B (trend
// only). The data partition the common component into trend vs cycle.
//
// Variant C nests A and B as restrictions:
//   - sigma_hec → 0 / lambda_eps → 0  ⇒  Variant B's common-trend-only structure
//   - sigma_hc  → 0 / lambda_tau → 0  ⇒  Variant A's common-cycle-only structure
//
// Differences from a faithful NY Fed port: MA(3) sector errors and the
// Student-t / outlier-scale-mixture likelihood are NOT included. Pilot tests
// showed HMC mixing breaks down with all four SV processes (common-trend,
// common-cycle, sector-trend, sector-noise) PLUS an MA filter PLUS a robust
// likelihood simultaneously active on simulated data — too richly
// parameterised for NUTS even with tight priors. NY Fed sidesteps this with
// a custom Gibbs sampler (slice/Kim-Shephard SV mixture updates per block).
// Re-introducing items 2 (MA) and 3 (robust likelihood) would require
// either centred SV parameterisations and very aggressive bounds, or moving
// off Stan entirely. Documented as v2 work.
//
// Loadings: BOTH alpha_tau (on tau_c) and alpha_eps (on c) are
// time-invariant. The NY Fed uses TVP loadings with very tight priors,
// allowing loadings to drift only ~0.08 over the whole sample. In our
// HMC setup, even modestly-prior'd TVPs on top of two common factors with
// SV creates multi-modal warmup behaviour — 50% of chains get stuck in
// alternate factor-allocation modes on simulated data. Fixing the
// loadings in time removes the structural ambiguity entirely while
// preserving cross-sector differences in factor sensitivity.
//
// Mixed-frequency observation handling matches A and B: quarterly
// observations placed at the last month of the quarter equal the simple
// average of three latent monthly inflations.
//
// Identification: alpha_tau[t, ref] = 1 for all t; alpha_eps[ref] = 1.
// All random walks parameterised non-centrally.

data {
  int<lower=1> T;
  int<lower=2> N;
  int<lower=1, upper=N> ref;
  vector<lower=0>[N] w;

  int<lower=0> n_obs_m;
  array[n_obs_m] int<lower=1, upper=T> t_m;
  array[n_obs_m] int<lower=1, upper=N> i_m;
  vector[n_obs_m] y_m;

  int<lower=0> n_obs_q;
  array[n_obs_q] int<lower=3, upper=T> t_q;
  array[n_obs_q] int<lower=1, upper=N> i_q;
  vector[n_obs_q] y_q;
}

parameters {
  // Log-variance RW hyperparams
  real mu_hc;                              // common-trend SV level
  real mu_hec;                             // common-cycle SV level
  vector[N] mu_hs;
  vector[N] mu_he;

  real<lower=0> sigma_hc;
  real<lower=0> sigma_hec;
  vector<lower=0>[N] sigma_hs;
  vector<lower=0>[N] sigma_he;

  // Initial conditions and time-invariant loadings on both common factors
  real tau_c_init;
  vector[N] s_init;
  vector[N - 1] lambda_tau_fixed;          // alpha_tau[i] for non-ref sectors
  vector[N - 1] lambda_eps_fixed;          // alpha_eps[i] for non-ref sectors

  // Non-centred RW innovations
  vector[T] z_hc;
  vector[T] z_hec;
  matrix[T, N] z_hs;
  matrix[T, N] z_he;
  matrix[T, N] z_s;
  vector[T] z_tau_c;
  vector[T] z_c;                           // common-cycle white-noise innovations
}

transformed parameters {
  vector[T] h_c;
  vector[T] h_ec;
  matrix[T, N] h_s;
  matrix[T, N] h_e;
  vector[T] tau_c;                         // common stochastic trend (RW)
  vector[T] c;                             // common cycle (i.i.d. with SV)
  matrix[T, N] s_trend;
  vector[N] lambda_tau;                    // time-invariant loadings on tau_c
  vector[N] lambda_eps;                    // time-invariant loadings on c

  // Common-trend log-vol RW
  h_c[1] = mu_hc + sigma_hc * z_hc[1];
  for (t in 2:T) h_c[t] = h_c[t - 1] + sigma_hc * z_hc[t];

  // Common-cycle log-vol RW
  h_ec[1] = mu_hec + sigma_hec * z_hec[1];
  for (t in 2:T) h_ec[t] = h_ec[t - 1] + sigma_hec * z_hec[t];

  // Sector trend / noise log-vol RWs
  for (i in 1:N) {
    h_s[1, i] = mu_hs[i] + sigma_hs[i] * z_hs[1, i];
    h_e[1, i] = mu_he[i] + sigma_he[i] * z_he[1, i];
    for (t in 2:T) {
      h_s[t, i] = h_s[t - 1, i] + sigma_hs[i] * z_hs[t, i];
      h_e[t, i] = h_e[t - 1, i] + sigma_he[i] * z_he[t, i];
    }
  }

  // Sector trends (RW with SV)
  for (i in 1:N) {
    s_trend[1, i] = s_init[i];
    for (t in 2:T)
      s_trend[t, i] = s_trend[t - 1, i] + exp(h_s[t, i] / 2) * z_s[t, i];
  }

  // Loadings: both time-invariant, ref pinned at 1
  lambda_tau[ref] = 1.0;
  lambda_eps[ref] = 1.0;
  {
    int k = 1;
    for (i in 1:N) {
      if (i != ref) {
        lambda_tau[i] = lambda_tau_fixed[k];
        lambda_eps[i] = lambda_eps_fixed[k];
        k += 1;
      }
    }
  }

  // Common factors
  //   tau_c: random walk anchored at tau_c_init, with stochastic volatility
  //   c:     i.i.d. white noise, scaled by its time-varying volatility
  tau_c[1] = tau_c_init;
  for (t in 2:T) tau_c[t] = tau_c[t - 1] + exp(h_c[t] / 2) * z_tau_c[t];
  for (t in 1:T) c[t] = exp(h_ec[t] / 2) * z_c[t];
}

model {
  // Hyperpriors — tighter than A/B (item 5). NY Fed's s2_gam = 1/720 implies
  // SD ~ 0.04 on the log-vol RW step; A/B use 0.1 which is ~3x looser.
  mu_hc  ~ normal(0, 1);
  mu_hec ~ normal(0, 1);
  mu_hs  ~ normal(0, 1);
  mu_he  ~ normal(0, 1);

  sigma_hc          ~ normal(0, 0.04);
  sigma_hec         ~ normal(0, 0.04);
  sigma_hs          ~ normal(0, 0.04);
  sigma_he          ~ normal(0, 0.04);

  tau_c_init       ~ normal(0, 2);
  s_init           ~ normal(0, 0.5);
  lambda_tau_fixed ~ normal(1, 0.3);
  lambda_eps_fixed ~ normal(1, 0.3);

  z_hc ~ std_normal();
  z_hec ~ std_normal();
  to_vector(z_hs) ~ std_normal();
  to_vector(z_he) ~ std_normal();
  to_vector(z_s) ~ std_normal();
  z_tau_c ~ std_normal();
  z_c ~ std_normal();

  // Monthly likelihood — Gaussian noise (no MA, no robust likelihood)
  if (n_obs_m > 0) {
    vector[n_obs_m] mu_m;
    vector[n_obs_m] sigma_m;
    for (k in 1:n_obs_m) {
      int t = t_m[k];
      int i = i_m[k];
      mu_m[k] = lambda_tau[i] * tau_c[t]
              + lambda_eps[i]    * c[t]
              + s_trend[t, i];
      sigma_m[k] = exp(h_e[t, i] / 2);
    }
    y_m ~ normal(mu_m, sigma_m);
  }

  // Quarterly likelihood — 3-month average of latent monthly inflations
  if (n_obs_q > 0) {
    vector[n_obs_q] mu_q;
    vector[n_obs_q] sigma_q;
    for (k in 1:n_obs_q) {
      int t = t_q[k];
      int i = i_q[k];
      mu_q[k] = (
          (lambda_tau[i] * tau_c[t - 2] + lambda_eps[i] * c[t - 2] + s_trend[t - 2, i])
        + (lambda_tau[i] * tau_c[t - 1] + lambda_eps[i] * c[t - 1] + s_trend[t - 1, i])
        + (lambda_tau[i] * tau_c[t]     + lambda_eps[i] * c[t]     + s_trend[t,     i])
      ) / 3.0;
      sigma_q[k] = sqrt(
        (exp(h_e[t - 2, i]) + exp(h_e[t - 1, i]) + exp(h_e[t, i])) / 9.0
      );
    }
    y_q ~ normal(mu_q, sigma_q);
  }
}

generated quantities {
  vector[T] trend;
  vector[T] common_share;
  vector[T] common_trend_level;            // level: sum_i w_i alpha_tau[t,i] tau_c[t]
  matrix[T, N] sector_trend;               // per-sector total trend = alpha_tau*tau_c + s_trend
  vector[T] common_transitory = c;
  vector[n_obs_m + n_obs_q] log_lik;

  for (t in 1:T)
    for (i in 1:N)
      sector_trend[t, i] = lambda_tau[i] * tau_c[t] + s_trend[t, i];

  // Pre-compute weighted-sum loading on tau_c (time-invariant)
  real wlam_tau_const = 0;
  for (i in 1:N) wlam_tau_const += w[i] * lambda_tau[i];

  // Weighted-aggregate trend and its common-vs-sector level decomposition
  for (t in 1:T) {
    real wlam_tau = 0;
    real ws = 0;
    for (i in 1:N) {
      wlam_tau += w[i] * lambda_tau[i];
      ws       += w[i] * s_trend[t, i];
    }
    common_trend_level[t] = wlam_tau * tau_c[t];
    trend[t] = common_trend_level[t] + ws;
  }

  // Common share as a variance ratio (matches A/B). Common-cycle variance
  // is NOT included — the "common share" question is about trend
  // co-movement, and tau_c is the common trend in C.
  for (t in 1:T) {
    real wlam_tau = 0;
    real sec_var = 0;
    for (i in 1:N) {
      wlam_tau += w[i] * lambda_tau[i];
      sec_var  += w[i] * w[i] * exp(h_s[t, i]);
    }
    real com_var = wlam_tau * wlam_tau * exp(h_c[t]);
    common_share[t] = com_var / (com_var + sec_var);
  }

  // Pointwise log-likelihood for LOO/WAIC
  for (k in 1:n_obs_m) {
    int t = t_m[k];
    int i = i_m[k];
    real mu = lambda_tau[i] * tau_c[t] + lambda_eps[i] * c[t] + s_trend[t, i];
    log_lik[k] = normal_lpdf(y_m[k] | mu, exp(h_e[t, i] / 2));
  }
  for (k in 1:n_obs_q) {
    int t = t_q[k];
    int i = i_q[k];
    real mu = (
        (lambda_tau[i] * tau_c[t - 2] + lambda_eps[i] * c[t - 2] + s_trend[t - 2, i])
      + (lambda_tau[i] * tau_c[t - 1] + lambda_eps[i] * c[t - 1] + s_trend[t - 1, i])
      + (lambda_tau[i] * tau_c[t]     + lambda_eps[i] * c[t]     + s_trend[t,     i])
    ) / 3.0;
    real var_obs = (
        exp(h_e[t - 2, i]) + exp(h_e[t - 1, i]) + exp(h_e[t, i])
    ) / 9.0;
    log_lik[n_obs_m + k] = normal_lpdf(y_q[k] | mu, sqrt(var_obs));
  }
}
