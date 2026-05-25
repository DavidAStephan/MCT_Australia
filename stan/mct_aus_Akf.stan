// Australian Multivariate Core Trend (MCT) inflation — Variant Akf
// (Phase 3: Kalman-filter marginalization of the linear latent states).
//
// **DRAFT — compile-clean target, math not yet validated.** See plan.md
// "Phase 3" for the literature anchor (Sørensen 2026, NUTS-Kalman pattern).
//
// Sampled parameters: hyperparams + log-vol RW innovations
// (z_hc, z_hs, z_he) + lambda_const. Total ~1,400 sampled at T=435 N=11
// — down from Ac's ~15,300. (z_c and z_s are MARGINALIZED, not sampled.)
//
// Marginalized via Kalman filter: the joint linear state
//   x_t = (c_t, s_{1,t}, ..., s_{N,t})
// is integrated out analytically, conditional on the SV paths
// (h_c, h_s, h_e) and the time-invariant loadings lambda.
//
// State augmentation: to handle the quarterly-average observation
// y^Q_{i,t} = (y_{i,t-2} + y_{i,t-1} + y_{i,t})/3 in a single-step
// observation equation, the state is augmented with 2 lags of every
// component. Augmented state dim:
//   D = 3 * (N + 1) = 36 at N = 11
//
// Indexing convention (1-based):
//   x[1] = c_t,    x[2] = c_{t-1},   x[3] = c_{t-2}
//   x[3i + 1] = s_{i,t},   x[3i + 2] = s_{i,t-1},   x[3i + 3] = s_{i,t-2}
//   for i = 1..N. So positions 4..36 cover the sector trend block.
//
// Transition (predicting x_{t+1} from x_t):
//   c_{t+1}      = rho * c_t + N(0, exp(h_c[t+1]))      [row 1]
//   c_t          carried (row 2 = old row 1)
//   c_{t-1}      carried (row 3 = old row 2)
//   s_{i,t+1}    = s_{i,t} + N(0, exp(h_s[i,t+1]))
//   s_{i,t}      carried
//   s_{i,t-1}    carried
//
// Monthly obs at (t, i):   y = lambda[i] * c_t + s_{i,t} + N(0, exp(h_e[i,t]))
//   Z row: lambda[i] at col 1, 1 at col 3i+1, else 0
//   H scalar: exp(h_e[i, t])
//
// Quarterly obs at (t, i): y^Q = lambda[i]*(c_{t-2}+c_{t-1}+c_t)/3
//                                + (s_{i,t-2}+s_{i,t-1}+s_{i,t})/3 + noise
//   Z^Q row: lambda[i]/3 at cols 1,2,3; 1/3 at cols 3i+1, 3i+2, 3i+3
//   H^Q scalar: (exp(h_e[i,t-2]) + exp(h_e[i,t-1]) + exp(h_e[i,t])) / 9
//
// **PRE-CONDITION**: observation arrays must be sorted by t ascending.
// Enforce in build_stan_data() before calling fit_mct(..., variant="Akf").

functions {
  // Forward Kalman filter returning the marginal log-likelihood.
  // Walks t = 1..T_, processing all observations at each t.
  real kalman_log_lik(
      int T_, int N, int ref,
      real rho,
      vector h_c, matrix h_s, matrix h_e,
      vector lambda, vector s_init,
      // Sorted-by-t observation arrays.
      int n_obs_m, array[] int t_m, array[] int i_m, vector y_m,
      int n_obs_q, array[] int t_q, array[] int i_q, vector y_q
  ) {
    int D = 3 * (N + 1);
    vector[D] x = rep_vector(0, D);     // state mean
    matrix[D, D] P = rep_matrix(0, D, D); // state covariance
    real loglik = 0;
    real log_2pi = 1.8378770664093453;

    // Initial state mean and covariance at t = 1.
    // c_1 ~ N(0, exp(h_c[1]) / (1 - rho^2))  (stationary AR(1) start)
    // c_0, c_{-1}: diffuse, var = 100 (acts as "no info" prior)
    // s_{i,1} = s_init[i], var = 0.25 (matching prior s_init ~ N(0, 0.5))
    // s_{i,0}, s_{i,-1}: diffuse, var = 100
    P[1, 1] = exp(h_c[1]) / (1.0 - rho * rho);
    P[2, 2] = 100.0;
    P[3, 3] = 100.0;
    for (i in 1:N) {
      x[3 * i + 1] = s_init[i];
      P[3 * i + 1, 3 * i + 1] = 0.25;
      P[3 * i + 2, 3 * i + 2] = 100.0;
      P[3 * i + 3, 3 * i + 3] = 100.0;
    }

    // Track pointers through the (sorted) observation arrays.
    int km = 1; // next monthly obs to process
    int kq = 1; // next quarterly obs to process

    for (t in 1:T_) {
      // --- Predict step (t > 1) ------------------------------------
      if (t > 1) {
        // Build transition. Sparse; only ~2(N+1)+1 non-zeros. Apply by
        // direct assignment rather than full matrix multiply to save ops.
        vector[D] x_pred = rep_vector(0, D);
        x_pred[1] = rho * x[1];
        x_pred[2] = x[1];
        x_pred[3] = x[2];
        for (i in 1:N) {
          x_pred[3 * i + 1] = x[3 * i + 1];           // s_{i,t} -> s_{i,t+1} (RW)
          x_pred[3 * i + 2] = x[3 * i + 1];
          x_pred[3 * i + 3] = x[3 * i + 2];
        }
        // T * P * T' — also sparse; build P_pred via cell-by-cell formula
        // to avoid the full O(D^3) matrix multiply. For now (sketch) use
        // dense ops; optimise later.
        matrix[D, D] Tm = rep_matrix(0, D, D);
        Tm[1, 1] = rho;
        Tm[2, 1] = 1;
        Tm[3, 2] = 1;
        for (i in 1:N) {
          Tm[3 * i + 1, 3 * i + 1] = 1;
          Tm[3 * i + 2, 3 * i + 1] = 1;
          Tm[3 * i + 3, 3 * i + 2] = 1;
        }
        matrix[D, D] P_pred = quad_form_sym(P, Tm');

        // Add process noise (diagonal, N+1 non-zero entries).
        P_pred[1, 1] += exp(h_c[t]);
        for (i in 1:N) {
          int p = 3 * i + 1;
          P_pred[p, p] += exp(h_s[t, i]);
        }
        x = x_pred;
        P = P_pred;
      }

      // --- Update with monthly obs at this t -----------------------
      while (km <= n_obs_m && t_m[km] == t) {
        int i = i_m[km];
        row_vector[D] Z = rep_row_vector(0, D);
        Z[1] = lambda[i];
        Z[3 * i + 1] = 1;
        real H = exp(h_e[t, i]);

        real v = y_m[km] - Z * x;                  // innovation
        real F = Z * P * Z' + H;                  // innovation variance (scalar)
        loglik += -0.5 * (log_2pi + log(F) + v * v / F);

        // Standard (non-Joseph) update; cheap, may need Joseph for stability.
        vector[D] K = P * Z' / F;
        x = x + K * v;
        P = P - K * (Z * P);

        km += 1;
      }

      // --- Update with quarterly obs at this t (needs t >= 3) ------
      while (kq <= n_obs_q && t_q[kq] == t) {
        int i = i_q[kq];
        row_vector[D] Z = rep_row_vector(0, D);
        Z[1] = lambda[i] / 3.0;
        Z[2] = lambda[i] / 3.0;
        Z[3] = lambda[i] / 3.0;
        Z[3 * i + 1] = 1.0 / 3.0;
        Z[3 * i + 2] = 1.0 / 3.0;
        Z[3 * i + 3] = 1.0 / 3.0;
        real H = (exp(h_e[t, i]) + exp(h_e[t - 1, i]) + exp(h_e[t - 2, i])) / 9.0;

        real v = y_q[kq] - Z * x;
        real F = Z * P * Z' + H;
        loglik += -0.5 * (log_2pi + log(F) + v * v / F);

        vector[D] K = P * Z' / F;
        x = x + K * v;
        P = P - K * (Z * P);

        kq += 1;
      }
    }

    return loglik;
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

  // SV innovations — same as Ac. z_c and z_s are GONE (marginalized).
  vector[T] z_hc;
  matrix[T, N] z_hs;
  matrix[T, N] z_he;
}

transformed parameters {
  vector[T] h_c;
  matrix[T, N] h_s;
  matrix[T, N] h_e;
  vector[N] lambda;

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
}

model {
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

  target += kalman_log_lik(
    T, N, ref, rho,
    h_c, h_s, h_e, lambda, s_init,
    n_obs_m, t_m, i_m, y_m,
    n_obs_q, t_q, i_q, y_q
  );
}

generated quantities {
  // TODO Phase 3.5: Kalman smoother to recover posterior of (c, s_trend)
  // for trend, common_share, sector_contributions, log_lik. Until then,
  // downstream targets (trend_Akf, common_share_Akf, ...) cannot be wired.
  //
  // Implementation sketch (Rauch-Tung-Striebel smoother):
  //   1. Re-run the forward filter, storing predictive and filtered
  //      means/covariances at every t.
  //   2. Backward recursion: for t = T-1, T-2, ..., 1:
  //        smoothed_t = filtered_t + J_t * (smoothed_{t+1} - predicted_{t+1})
  //        where J_t = filtered_P_t * T' * inv(predicted_P_{t+1})
  //   3. Extract c_t from smoothed[1], s_{i,t} from smoothed[3i+1].
  //   4. Then compute trend, common_share, log_lik as in mct_aus_Ac.stan
  //      generated quantities.

  // Placeholder so cmdstanr doesn't complain about missing outputs.
  vector[T] trend = rep_vector(0, T);
  vector[T] common_share = rep_vector(0, T);
}
