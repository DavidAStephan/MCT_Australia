// Australian Multivariate Core Trend (MCT) inflation — Variant Ac
// (Variant A with time-INVARIANT loadings; constant lambda).
//
// Identical to Variant A except lambda[i] is a single scalar per non-
// reference sector — no time variation. This removes ~T*(N-1) = 4,350
// sampled parameters at T=435 N=11 (vs Variant A's z_lambda matrix), and
// removes the corresponding RW reconstruction work in transformed
// parameters.
//
// Justification: the existing model has sigma_lambda ~ N(0, 0.01) which
// already implies negligible drift over T=435 months. Effectively the
// loadings barely move under the prior; making them constant costs little
// modelling flexibility but a lot of sampling time.
//
// All other structure (mixed-frequency obs, AR(1) common with SV, RW
// sector trends with SV, time-varying log-vols) is unchanged from
// stan/mct_aus_A.stan.
//
// Note: profile() blocks were initially added for #4 timing diagnosis but
// removed after a 5x per-iter slowdown was observed in the first benchmark.
// Stan's profile() runs once per leapfrog step (not per iteration); with
// max_treedepth=12 the overhead compounds severely. Re-add only for
// targeted diagnostic runs, never for production fits.

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
  real mu_hc;
  vector[N] mu_hs;
  vector[N] mu_he;

  real<lower=0> sigma_hc;
  vector<lower=0>[N] sigma_hs;
  vector<lower=0>[N] sigma_he;

  real<lower=-1, upper=1> rho;

  vector[N] s_init;
  vector[N - 1] lambda_const;   // time-invariant non-reference loadings

  vector[T] z_hc;
  matrix[T, N] z_hs;
  matrix[T, N] z_he;
  matrix[T, N] z_s;
  vector[T] z_c;
}

transformed parameters {
  vector[T] h_c;
  matrix[T, N] h_s;
  matrix[T, N] h_e;
  matrix[T, N] s_trend;
  vector[N] lambda;             // constant across t — N-vector
  vector[T] c;

  h_c[1] = mu_hc + sigma_hc * z_hc[1];
  for (t in 2:T) h_c[t] = h_c[t - 1] + sigma_hc * z_hc[t];

  for (i in 1:N) {
    h_s[1, i] = mu_hs[i] + sigma_hs[i] * z_hs[1, i];
    h_e[1, i] = mu_he[i] + sigma_he[i] * z_he[1, i];
    for (t in 2:T) {
      h_s[t, i] = h_s[t - 1, i] + sigma_hs[i] * z_hs[t, i];
      h_e[t, i] = h_e[t - 1, i] + sigma_he[i] * z_he[t, i];
    }
  }

  for (i in 1:N) {
    s_trend[1, i] = s_init[i];
    for (t in 2:T)
      s_trend[t, i] = s_trend[t - 1, i] + exp(h_s[t, i] / 2) * z_s[t, i];
  }

  // Time-invariant loadings: scalar per non-ref sector, ref pinned at 1.
  lambda[ref] = 1.0;
  {
    int k = 1;
    for (i in 1:N) {
      if (i != ref) {
        lambda[i] = lambda_const[k];
        k += 1;
      }
    }
  }

  c[1] = exp(h_c[1] / 2) / sqrt(1 - rho * rho) * z_c[1];
  for (t in 2:T) c[t] = rho * c[t - 1] + exp(h_c[t] / 2) * z_c[t];
}

model {
  // Hyperpriors (same as Variant A)
  mu_hc ~ normal(0, 1);
  mu_hs ~ normal(0, 1);
  mu_he ~ normal(0, 1);
  sigma_hc ~ normal(0, 0.1);
  sigma_hs ~ normal(0, 0.1);
  sigma_he ~ normal(0, 0.1);
  rho ~ normal(0.5, 0.3);

  s_init ~ normal(0, 0.5);
  lambda_const ~ normal(1, 0.3);

  z_hc ~ std_normal();
  to_vector(z_hs) ~ std_normal();
  to_vector(z_he) ~ std_normal();
  to_vector(z_s) ~ std_normal();
  z_c ~ std_normal();

  if (n_obs_m > 0) {
    vector[n_obs_m] mu_m;
    vector[n_obs_m] sigma_m;
    for (k in 1:n_obs_m) {
      mu_m[k] = lambda[i_m[k]] * c[t_m[k]] + s_trend[t_m[k], i_m[k]];
      sigma_m[k] = exp(h_e[t_m[k], i_m[k]] / 2);
    }
    y_m ~ normal(mu_m, sigma_m);
  }

  if (n_obs_q > 0) {
    vector[n_obs_q] mu_q;
    vector[n_obs_q] sigma_q;
    for (k in 1:n_obs_q) {
      int t = t_q[k];
      int i = i_q[k];
      mu_q[k] = (
          lambda[i] * c[t - 2] + s_trend[t - 2, i]
        + lambda[i] * c[t - 1] + s_trend[t - 1, i]
        + lambda[i] * c[t]     + s_trend[t,     i]
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
  matrix[T, N] sector_trend = s_trend;
  vector[T] common_transitory = c;
  vector[n_obs_m + n_obs_q] log_lik;

  // Trend = weighted average of sector-specific trends (unchanged from A).
  for (t in 1:T) {
    real tt = 0;
    for (i in 1:N) tt += w[i] * s_trend[t, i];
    trend[t] = tt;
  }

  // Common share: lambda is constant across t now, so wlam is too — pull
  // out of the loop. SV variance evolution still gives a time-varying ratio.
  {
    real wlam = 0;
    for (i in 1:N) wlam += w[i] * lambda[i];
    real wlam2 = wlam * wlam;
    for (t in 1:T) {
      real sec_var = 0;
      for (i in 1:N) sec_var += w[i] * w[i] * exp(h_s[t, i]);
      real com_var = wlam2 * exp(h_c[t]);
      common_share[t] = com_var / (com_var + sec_var);
    }
  }

  // Pointwise log-likelihood for LOO. lambda[i] is scalar (no t index).
  for (k in 1:n_obs_m) {
    real mu = lambda[i_m[k]] * c[t_m[k]] + s_trend[t_m[k], i_m[k]];
    log_lik[k] = normal_lpdf(y_m[k] | mu, exp(h_e[t_m[k], i_m[k]] / 2));
  }
  for (k in 1:n_obs_q) {
    int t = t_q[k];
    int i = i_q[k];
    real mu = (
        lambda[i] * c[t - 2] + s_trend[t - 2, i]
      + lambda[i] * c[t - 1] + s_trend[t - 1, i]
      + lambda[i] * c[t]     + s_trend[t,     i]
    ) / 3.0;
    real var_obs = (
        exp(h_e[t - 2, i]) + exp(h_e[t - 1, i]) + exp(h_e[t, i])
    ) / 9.0;
    log_lik[n_obs_m + k] = normal_lpdf(y_q[k] | mu, sqrt(var_obs));
  }
}
