// Australian Multivariate Core Trend (MCT) inflation — Variant Act
// (Variant Ac + within-chain parallelism via reduce_sum on the
// likelihood). Phase 2 #3 in plan.md.
//
// Identical to mct_aus_Ac.stan except the monthly and quarterly
// likelihood loops are wrapped in `reduce_sum_lpdf` partial sums.
// Stan distributes the partial sums across `threads_per_chain` worker
// threads at runtime.
//
// Compile with cpp_options = list(stan_threads = TRUE).
// Sample with threads_per_chain = K (typically 2-4 on M1).
//
// The monthly and quarterly obs are conditionally independent given the
// latent paths (c, s_trend, h_e, lambda), so two independent reduce_sum
// calls is correct. Grainsize defaults to 1 — Stan auto-tunes the
// partition size.

functions {
  // Sum monthly observation log-likelihoods for a slice of observation
  // indices in slice_k (each value is an index into y_m / t_m / i_m).
  // Function name has no _lpdf suffix because the first arg is an int
  // array (the slice), not the variate. Stan reserves _lpdf for functions
  // whose first arg is the real-valued variate.
  real partial_sum_monthly(
      array[] int slice_k, int start, int end,
      vector y_m, array[] int t_m, array[] int i_m,
      vector lambda, vector c, matrix s_trend, matrix h_e) {
    real lp = 0;
    for (n in 1:size(slice_k)) {
      int k = slice_k[n];
      int t = t_m[k];
      int i = i_m[k];
      real mu = lambda[i] * c[t] + s_trend[t, i];
      real sigma = exp(h_e[t, i] / 2);
      lp += normal_lpdf(y_m[k] | mu, sigma);
    }
    return lp;
  }

  // Same shape for quarterly — three-month accumulator inside the loop.
  real partial_sum_quarterly(
      array[] int slice_k, int start, int end,
      vector y_q, array[] int t_q, array[] int i_q,
      vector lambda, vector c, matrix s_trend, matrix h_e) {
    real lp = 0;
    for (n in 1:size(slice_k)) {
      int k = slice_k[n];
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
      lp += normal_lpdf(y_q[k] | mu, sqrt(var_obs));
    }
    return lp;
  }
}

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

transformed data {
  // Pre-built index arrays for the reduce_sum calls. Stan slices these
  // and passes each chunk to the partial_sum_* functions.
  array[n_obs_m] int idx_m;
  array[n_obs_q] int idx_q;
  for (k in 1:n_obs_m) idx_m[k] = k;
  for (k in 1:n_obs_q) idx_q[k] = k;
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
  vector[N - 1] lambda_const;

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
  vector[N] lambda;
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
  // Hyperpriors (identical to Ac)
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

  // Parallel-summed likelihood. grainsize=1 lets Stan auto-tune the
  // chunk size at compile/runtime.
  if (n_obs_m > 0) {
    target += reduce_sum(partial_sum_monthly, idx_m, 1,
                         y_m, t_m, i_m, lambda, c, s_trend, h_e);
  }
  if (n_obs_q > 0) {
    target += reduce_sum(partial_sum_quarterly, idx_q, 1,
                         y_q, t_q, i_q, lambda, c, s_trend, h_e);
  }
}

generated quantities {
  vector[T] trend;
  vector[T] common_share;
  matrix[T, N] sector_trend = s_trend;
  vector[T] common_transitory = c;
  vector[n_obs_m + n_obs_q] log_lik;

  for (t in 1:T) {
    real tt = 0;
    for (i in 1:N) tt += w[i] * s_trend[t, i];
    trend[t] = tt;
  }

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
