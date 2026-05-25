// Australian Multivariate Core Trend (MCT) inflation — Variant Akfs
// (Akf with SPARSE Kalman filter implementation).
//
// Same model as mct_aus_Akf.stan; same marginalization; same sampled
// parameters. The difference is that the per-step KF operations exploit
// the sparsity of the transition matrix T and observation matrix Z:
//
//   Akf  (dense):  quad_form_sym(P, T') = O(D^3) = ~47K ops at D=36
//   Akfs (sparse): cell-by-cell predict  = O(D^2) = ~1.3K ops at D=36
//                  → ~36x cheaper per filter step (autodiff tape too)
//
// Why this works: every row of T has at most 1 non-zero entry. So
// (T * P * T')[a, b] = T[a, u_a] * P[u_a, u_b] * T[b, u_b] where
// (u_a, t_a = T[a, u_a]) is precomputed in transformed data.
//
// Observation matrix Z for monthly obs at (t, i) has 2 non-zeros
// (lambda[i] at col 1, 1 at col 3i+1); for quarterly, 6 non-zeros.
// Z * P * Z', P * Z', and (I - K*Z) * P are all computed directly
// without ever materialising Z as a dense row vector.

functions {
  // Sparse forward Kalman filter. Returns marginal log-likelihood.
  // Pre-condition: obs arrays t_m / t_q must be sorted ascending.
  real kalman_log_lik(
      int T_, int N, int ref,
      real rho,
      vector h_c, matrix h_s, matrix h_e,
      vector lambda, vector s_init,
      array[] int u_a, vector t_a,   // sparse-T encoding
      int n_obs_m, array[] int t_m, array[] int i_m, vector y_m,
      int n_obs_q, array[] int t_q, array[] int i_q, vector y_q
  ) {
    int D = 3 * (N + 1);
    vector[D] x = rep_vector(0, D);
    matrix[D, D] P = rep_matrix(0, D, D);
    real loglik = 0;
    real log_2pi = 1.8378770664093453;

    // Initial state — same as Akf for direct comparability.
    P[1, 1] = exp(h_c[1]) / (1.0 - rho * rho);
    P[2, 2] = 100.0;
    P[3, 3] = 100.0;
    for (i in 1:N) {
      x[3 * i + 1] = s_init[i];
      P[3 * i + 1, 3 * i + 1] = 0.25;
      P[3 * i + 2, 3 * i + 2] = 100.0;
      P[3 * i + 3, 3 * i + 3] = 100.0;
    }

    int km = 1;
    int kq = 1;

    for (t in 1:T_) {
      // --- Predict step (t > 1) using sparse T -----------------------
      if (t > 1) {
        // x_pred[a] = T[a, u_a] * x[u_a] = t_a[a] * x[u_a[a]]
        vector[D] x_pred;
        for (a in 1:D) x_pred[a] = t_a[a] * x[u_a[a]];

        // P_pred[a, b] = t_a[a] * t_a[b] * P[u_a[a], u_a[b]]
        // This is the central O(D^2) win vs quad_form_sym's O(D^3).
        matrix[D, D] P_pred;
        for (a in 1:D)
          for (b in 1:D)
            P_pred[a, b] = t_a[a] * t_a[b] * P[u_a[a], u_a[b]];

        // Add process noise to diagonal (only N+1 non-zero entries).
        P_pred[1, 1] += exp(h_c[t]);
        for (i in 1:N) {
          int p = 3 * i + 1;
          P_pred[p, p] += exp(h_s[t, i]);
        }
        x = x_pred;
        P = P_pred;
      }

      // --- Update with monthly obs at this t -------------------------
      // Z has 2 non-zeros for monthly at (t, i): Z[1]=lambda[i], Z[3i+1]=1.
      // Compute v, F, K, and P-update directly without materialising Z.
      while (km <= n_obs_m && t_m[km] == t) {
        int i = i_m[km];
        int p = 3 * i + 1;     // position of s_{i,t} in the state
        real li = lambda[i];

        // Innovation: v = y - lambda[i] * x[1] - x[p]
        real v = y_m[km] - li * x[1] - x[p];

        // F = Z * P * Z' + H. With Z = [li 0 0 ... 0 1 0 0 ...]:
        //   F = li^2 * P[1,1] + 2 * li * P[1, p] + P[p, p] + H
        real F = li * li * P[1, 1] + 2.0 * li * P[1, p] + P[p, p]
                 + exp(h_e[t, i]);

        // K[a] = (P * Z')[a] / F = (li * P[a, 1] + P[a, p]) / F
        vector[D] K;
        for (a in 1:D) K[a] = (li * P[a, 1] + P[a, p]) / F;

        loglik += -0.5 * (log_2pi + log(F) + v * v / F);

        // x_new = x + K * v
        x = x + K * v;

        // P_new = P - K * (Z * P)
        // (Z * P)[b] = li * P[1, b] + P[p, b]
        // P_new[a, b] = P[a, b] - K[a] * (li * P[1, b] + P[p, b])
        matrix[D, D] P_new;
        for (a in 1:D)
          for (b in 1:D)
            P_new[a, b] = P[a, b] - K[a] * (li * P[1, b] + P[p, b]);
        P = P_new;

        km += 1;
      }

      // --- Update with quarterly obs at this t -----------------------
      // Z has 6 non-zeros for quarterly at (t, i):
      //   Z[1]=Z[2]=Z[3]=li/3; Z[p]=Z[p+1]=Z[p+2]=1/3 (p = 3i+1)
      while (kq <= n_obs_q && t_q[kq] == t) {
        int i = i_q[kq];
        int p = 3 * i + 1;
        real li3 = lambda[i] / 3.0;
        real one3 = 1.0 / 3.0;

        // ZX = sum of Z[a] * x[a] for the 6 non-zero positions
        real ZX = li3 * (x[1] + x[2] + x[3])
                  + one3 * (x[p] + x[p + 1] + x[p + 2]);
        real v = y_q[kq] - ZX;

        // (P * Z')[a] = sum of Z[u] * P[a, u]
        //            = li3 * (P[a,1] + P[a,2] + P[a,3])
        //              + one3 * (P[a,p] + P[a,p+1] + P[a,p+2])
        vector[D] PZ;
        for (a in 1:D)
          PZ[a] = li3 * (P[a, 1] + P[a, 2] + P[a, 3])
                  + one3 * (P[a, p] + P[a, p + 1] + P[a, p + 2]);

        // F = Z * (P * Z') + H
        real ZPZ = li3 * (PZ[1] + PZ[2] + PZ[3])
                   + one3 * (PZ[p] + PZ[p + 1] + PZ[p + 2]);
        real H = (exp(h_e[t, i]) + exp(h_e[t - 1, i]) + exp(h_e[t - 2, i]))
                 / 9.0;
        real F = ZPZ + H;

        vector[D] K = PZ / F;

        loglik += -0.5 * (log_2pi + log(F) + v * v / F);

        x = x + K * v;

        // P_new = P - K * (Z * P)
        // (Z * P)[b] = li3 * (P[1,b] + P[2,b] + P[3,b])
        //              + one3 * (P[p,b] + P[p+1,b] + P[p+2,b])
        // Since P is symmetric, (Z * P) = (P * Z')' = PZ', so use PZ.
        matrix[D, D] P_new;
        for (a in 1:D)
          for (b in 1:D)
            P_new[a, b] = P[a, b] - K[a] * PZ[b];
        P = P_new;

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

transformed data {
  // Pre-encode the sparse transition matrix T. Each row has at most one
  // non-zero entry. u_a[a] = column of the non-zero in row a; t_a[a] = the
  // value (rho for row 1, 1 elsewhere). With rho moved here, we'd lose
  // its derivative — instead, t_a[1] is the rho-row marker. Read further:
  //   For row 1 (c_t = rho * c_{t-1}), we set u_a[1] = 1 and put `rho`
  //   into t_a[1] inside the filter (which depends on the parameter).
  // To keep transformed data parameter-free, encode just u_a here and
  // build t_a inside the filter using rho.
  int D = 3 * (N + 1);
  array[D] int u_a;

  u_a[1] = 1; // c_t row: depends on c_{t-1}, coeff rho (set in filter)
  u_a[2] = 1; // c_{t-1} carry
  u_a[3] = 2; // c_{t-2} carry
  for (i in 1:N) {
    u_a[3 * i + 1] = 3 * i + 1; // s_{i,t} carry (RW)
    u_a[3 * i + 2] = 3 * i + 1; // s_{i,t-1} = old s_{i,t}
    u_a[3 * i + 3] = 3 * i + 2; // s_{i,t-2} = old s_{i,t-1}
  }
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

  // Build t_a now that we have rho. Cheap; D=3(N+1) scalar assignments.
  vector[3 * (N + 1)] t_a;
  t_a[1] = rho;        // c_t = rho * c_{t-1}
  t_a[2] = 1;
  t_a[3] = 1;
  for (i in 1:N) {
    t_a[3 * i + 1] = 1; // s_{i,t} = s_{i,t-1} (RW)
    t_a[3 * i + 2] = 1;
    t_a[3 * i + 3] = 1;
  }

  target += kalman_log_lik(
    T, N, ref, rho,
    h_c, h_s, h_e, lambda, s_init,
    u_a, t_a,
    n_obs_m, t_m, i_m, y_m,
    n_obs_q, t_q, i_q, y_q
  );
}

generated quantities {
  // TODO Phase 3c: RTS smoother. Placeholder for now so cmdstanr accepts
  // the model and downstream targets at least don't crash.
  vector[T] trend = rep_vector(0, T);
  vector[T] common_share = rep_vector(0, T);
}
