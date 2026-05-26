# Mixed-frequency state-space-model constructor for the AU MCT Variant A.
#
# Step 9a of gibbs_port_plan.md. Adapts build_mct_ssm() to accept a
# per-(t, i) obs_type indicator (0 = missing, 1 = monthly, 2 = quarterly
# average) and emits a state-space model with the joint state augmented
# by 2 lags so the quarterly-average observation y^Q_{i,t} =
# (y_{i,t-2} + y_{i,t-1} + y_{i,t})/3 fits in a single-step linear obs.
#
# State (dim D = 3(N+1) = 36 at N=11):
#   x[1..3]         = (c_t, c_{t-1}, c_{t-2})
#   x[3i+1..3i+3]   = (s_{i,t}, s_{i,t-1}, s_{i,t-2})   for i = 1..N
#
# Transition (sparse, ~3(N+1)+1 non-zeros):
#   c_{t+1}     = rho * c_t + noise (sigma_c[t+1])
#   c_t carry; c_{t-1} carry
#   s_{i,t+1}   = s_{i,t} + noise (sigma_s[t+1, i])
#   s_{i,t}, s_{i,t-1} carry
#
# Observation (time-varying H):
#   monthly  at (t,i): y = lambda[i] * c_t + s_{i,t} + noise
#                      noise SD = sigma_eps[t, i]
#   quarterly at (t,i): y = lambda[i] * (c_t + c_{t-1} + c_{t-2})/3 +
#                           (s_{i,t} + s_{i,t-1} + s_{i,t-2})/3 + noise
#                      noise variance = (sigma_eps[t,i]^2 +
#                                         sigma_eps[t-1,i]^2 +
#                                         sigma_eps[t-2,i]^2) / 9
#   missing at (t,i): y_{i,t} = NA; H[i, :, t] is don't-care.

#' Build a mixed-frequency SSM for the AU MCT Variant A model.
#'
#' @param rho       Scalar AR(1) coef on c.
#' @param lambda    Length-N vector of time-invariant loadings.
#' @param sigma_c   Length-T vector of common-factor innovation SDs.
#' @param sigma_s   T x N matrix of sector-trend innovation SDs.
#' @param sigma_eps T x N matrix of measurement-noise SDs.
#' @param obs_type  T x N integer matrix; entries are 0 (missing),
#'                  1 (monthly), or 2 (quarterly average; requires t >= 3).
#' @param T_        Number of time periods.
#' @param N         Number of sectors.
#' @param var_s_init Prior variance on s_{i,1}. Default 0.25.
#' @return List of SSM matrices with augmented state, ready to pass
#'   to kalman_filter / fast_smoother / simulation_smoother.
build_mct_ssm_mixed <- function(rho, lambda, sigma_c, sigma_s, sigma_eps,
                                obs_type, T_, N, var_s_init = 0.25) {
  stopifnot(length(lambda) == N,
            length(sigma_c) == T_,
            all(dim(sigma_s) == c(T_, N)),
            all(dim(sigma_eps) == c(T_, N)),
            all(dim(obs_type) == c(T_, N)))

  D <- 3L * (N + 1L)
  K <- N + 1L                # number of independent innovations

  # --- Transition F (sparse) -----------------------------------------
  F_ <- matrix(0, D, D)
  F_[1L, 1L] <- rho
  F_[2L, 1L] <- 1
  F_[3L, 2L] <- 1
  for (i in seq_len(N)) {
    p <- 3L * i + 1L
    F_[p,     p] <- 1
    F_[p + 1L, p] <- 1
    F_[p + 2L, p + 1L] <- 1
  }

  # --- Disturbance loadings G (D x K) --------------------------------
  # Only the "current-t" positions receive fresh noise.
  G <- matrix(0, D, K)
  G[1L, 1L] <- 1
  for (i in seq_len(N)) G[3L * i + 1L, 1L + i] <- 1

  # --- Process noise covariance (K x K x T-1) ------------------------
  Sigma_eta <- array(0, c(K, K, T_ - 1L))
  for (t in seq_len(T_ - 1L)) {
    Sigma_eta[, , t] <- diag(c(sigma_c[t + 1L]^2, sigma_s[t + 1L, ]^2))
  }

  # --- Time-varying observation matrix H (N x D x T) -----------------
  # And per-time obs-noise covariance (N x N x T)
  H <- array(0, c(N, D, T_))
  Sigma_eps_3d <- array(0, c(N, N, T_))
  for (t in seq_len(T_)) {
    for (i in seq_len(N)) {
      ot <- obs_type[t, i]
      if (ot == 1L) {                            # monthly
        H[i, 1L, t]              <- lambda[i]
        H[i, 3L * i + 1L, t]     <- 1
        Sigma_eps_3d[i, i, t]    <- sigma_eps[t, i]^2
      } else if (ot == 2L && t >= 3L) {          # quarterly avg
        H[i, 1L:3L, t]           <- lambda[i] / 3
        H[i, 3L * i + 1L:3L, t]  <- 1 / 3
        Sigma_eps_3d[i, i, t]    <- (sigma_eps[t,     i]^2 +
                                     sigma_eps[t - 1L, i]^2 +
                                     sigma_eps[t - 2L, i]^2) / 9
      } else {
        # missing — H[i, :, t] stays zero; Sigma_eps_3d[i, i, t] = small
        # nonzero to keep S = HPH' + Sigma_eps invertible when the obs
        # at this (t, i) is NA. The kalman_filter ignores Sigma_eps
        # entries on missing rows, so the value here doesn't matter
        # numerically — but we set a tiny floor for safety.
        Sigma_eps_3d[i, i, t] <- 1e-12
      }
    }
  }

  # --- Initial state mean + covariance -------------------------------
  # Diffuse priors for c lags and s lags (no t < 1 info); tight prior
  # on s_{i,1} (matches Stan variant Ac s_init ~ N(0, 0.5^2)).
  mu_1 <- rep(0, D)
  Sigma_1 <- diag(D)
  Sigma_1[1L, 1L] <- 100
  Sigma_1[2L, 2L] <- 100
  Sigma_1[3L, 3L] <- 100
  for (i in seq_len(N)) {
    Sigma_1[3L * i + 1L, 3L * i + 1L] <- var_s_init
    Sigma_1[3L * i + 2L, 3L * i + 2L] <- 100
    Sigma_1[3L * i + 3L, 3L * i + 3L] <- 100
  }

  list(
    H         = H,
    F         = F_,
    G         = G,
    Sigma_eps = Sigma_eps_3d,
    Sigma_eta = Sigma_eta,
    mu_1      = mu_1,
    Sigma_1   = Sigma_1
  )
}
