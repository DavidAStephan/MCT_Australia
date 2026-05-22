// Australian Multivariate Core Trend (MCT) inflation — Variant B
// (Random-walk common trend). Stock-Watson 2016 / NY Fed MCT analogue —
// this is the closest direct analogue: c is the common stochastic TREND,
// not a transitory factor.
//
// Mixed-frequency observation: a quarterly observation placed at the last
// month of the quarter equals the simple average of the three latent monthly
// inflations in months t-2, t-1, t.
//
// Identification: lambda[t, ref] = 1 for all t (ref = Housing).
// Common factor: random walk; diffuse initial distribution c[1] ~ N(0, 2).
//
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
  real mu_hc;
  vector[N] mu_hs;
  vector[N] mu_he;

  real<lower=0> sigma_hc;
  vector<lower=0>[N] sigma_hs;
  vector<lower=0>[N] sigma_he;
  vector<lower=0>[N - 1] sigma_lambda;

  // No rho in Variant B (common factor is RW, not AR(1)).

  vector[N] s_init;
  vector[N - 1] lambda_init;

  vector[T] z_hc;
  matrix[T, N] z_hs;
  matrix[T, N] z_he;
  matrix[T, N] z_s;
  matrix[T, N - 1] z_lambda;
  vector[T] z_c;
}

transformed parameters {
  vector[T] h_c;
  matrix[T, N] h_s;
  matrix[T, N] h_e;
  matrix[T, N] s_trend;
  matrix[T, N] lambda;
  vector[T] c;

  // Log-volatilities
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

  // Sector trends
  for (i in 1:N) {
    s_trend[1, i] = s_init[i];
    for (t in 2:T)
      s_trend[t, i] = s_trend[t - 1, i] + exp(h_s[t, i] / 2) * z_s[t, i];
  }

  // Loadings
  for (t in 1:T) lambda[t, ref] = 1.0;
  {
    int k = 1;
    for (i in 1:N) {
      if (i != ref) {
        lambda[1, i] = lambda_init[k];
        for (t in 2:T)
          lambda[t, i] = lambda[t - 1, i] + sigma_lambda[k] * z_lambda[t, k];
        k += 1;
      }
    }
  }

  // Common factor: random walk with diffuse initial c[1] ~ N(0, 2^2).
  c[1] = 2.0 * z_c[1];
  for (t in 2:T) c[t] = c[t - 1] + exp(h_c[t] / 2) * z_c[t];
}

model {
  // Hyperpriors. Tightened from the brief's defaults — see mct_aus_A.stan
  // for diagnosis & rationale.
  mu_hc ~ normal(0, 1);
  mu_hs ~ normal(0, 1);
  mu_he ~ normal(0, 1);
  sigma_hc ~ normal(0, 0.1);
  sigma_hs ~ normal(0, 0.1);
  sigma_he ~ normal(0, 0.1);
  sigma_lambda ~ normal(0, 0.01);

  s_init ~ normal(0, 0.5);
  lambda_init ~ normal(1, 0.3);

  z_hc ~ std_normal();
  to_vector(z_hs) ~ std_normal();
  to_vector(z_he) ~ std_normal();
  to_vector(z_s) ~ std_normal();
  to_vector(z_lambda) ~ std_normal();
  z_c ~ std_normal();

  if (n_obs_m > 0) {
    vector[n_obs_m] mu_m;
    vector[n_obs_m] sigma_m;
    for (k in 1:n_obs_m) {
      mu_m[k] = lambda[t_m[k], i_m[k]] * c[t_m[k]] + s_trend[t_m[k], i_m[k]];
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
          lambda[t - 2, i] * c[t - 2] + s_trend[t - 2, i]
        + lambda[t - 1, i] * c[t - 1] + s_trend[t - 1, i]
        + lambda[t,     i] * c[t]     + s_trend[t,     i]
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
  vector[T] common_transitory;       // = 0 in Variant B (no transitory common)
  vector[n_obs_m + n_obs_q] log_lik;

  // Variant B trend = (w'lambda_t) * c_t + w' s_trend_t.
  // The common factor enters trend directly (it IS the common trend).
  for (t in 1:T) {
    real wlam = 0;
    real ws = 0;
    for (i in 1:N) {
      wlam += w[i] * lambda[t, i];
      ws += w[i] * s_trend[t, i];
    }
    trend[t] = wlam * c[t] + ws;
    common_transitory[t] = 0;
  }

  // Common share of trend INNOVATION variance — canonical Stock-Watson.
  for (t in 1:T) {
    real wlam = 0;
    real sec_var = 0;
    for (i in 1:N) {
      wlam += w[i] * lambda[t, i];
      sec_var += w[i] * w[i] * exp(h_s[t, i]);
    }
    real com_var = wlam * wlam * exp(h_c[t]);
    common_share[t] = com_var / (com_var + sec_var);
  }

  for (k in 1:n_obs_m) {
    real mu = lambda[t_m[k], i_m[k]] * c[t_m[k]] + s_trend[t_m[k], i_m[k]];
    log_lik[k] = normal_lpdf(y_m[k] | mu, exp(h_e[t_m[k], i_m[k]] / 2));
  }
  for (k in 1:n_obs_q) {
    int t = t_q[k];
    int i = i_q[k];
    real mu = (
        lambda[t - 2, i] * c[t - 2] + s_trend[t - 2, i]
      + lambda[t - 1, i] * c[t - 1] + s_trend[t - 1, i]
      + lambda[t,     i] * c[t]     + s_trend[t,     i]
    ) / 3.0;
    real var_obs = (
        exp(h_e[t - 2, i]) + exp(h_e[t - 1, i]) + exp(h_e[t, i])
    ) / 9.0;
    log_lik[n_obs_m + k] = normal_lpdf(y_q[k] | mu, sqrt(var_obs));
  }
}
