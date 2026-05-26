# State-augmented SSM for MCT Variant A with MA(q) measurement errors.
#
# Monthly-only mode for v1 (Step 12). Mixed-frequency + MA combined
# is a future extension — the state dim would be 3(N+1) + (q+1)*N
# which gets heavy.
#
# Model:
#   y_{i,t} = lambda_i * c_t + s_{i,t} +
#             sigma_eps_{i,t} * (eps_{i,t} + theta_1 eps_{i,t-1}
#                                + ... + theta_q eps_{i,t-q})
#   c_t        = rho * c_{t-1} + sigma_c[t] * nu_c
#   s_{i,t}    = s_{i,t-1} + sigma_s[t,i] * nu_s_i
#   eps_{i,t}  iid N(0, 1)
#
# State (dim D = 1 + N + (q+1)*N at monthly):
#   x[1]                                              = c_t
#   x[1 + i]                                          = s_{i,t}    for i=1..N
#   x[1 + N + (q+1)*(i-1) + 1..1 + N + (q+1)*i]       = eps_{i,t}, eps_{i,t-1}, ..., eps_{i,t-q}
#
# By folding the MA structure into the state, the obs noise becomes
# (numerically near-)zero — the entire eps-MA expansion is part of
# the H matrix.

#' Build the MA(q)-augmented SSM (monthly-obs only).
#'
#' @param rho       Scalar AR(1) coef on c.
#' @param lambda    Length-N vector of time-invariant loadings.
#' @param sigma_c   Length-T vector of common-factor innovation SDs.
#' @param sigma_s   T x N matrix of sector-trend innovation SDs.
#' @param sigma_eps T x N matrix of measurement-noise SDs.
#' @param theta     N x q matrix of MA coefficients (per sector).
#' @param T_        Number of time periods.
#' @param N         Number of sectors.
#' @param q_MA      MA order (number of MA lags).
#' @param var_s_init Prior variance on s_{i,1}. Default 0.25.
#' @param obs_noise_floor Tiny diagonal added to Sigma_eps for numerical
#'   stability (defaults 1e-8); the true obs noise is folded into the state.
#' @return List of SSM matrices ready for kalman_filter / simulation_smoother.
build_mct_ssm_ma <- function(rho, lambda, sigma_c, sigma_s, sigma_eps,
                             theta, T_, N, q_MA = 3L,
                             var_s_init = 0.25,
                             obs_noise_floor = 1e-8) {
  stopifnot(length(lambda) == N,
            length(sigma_c) == T_,
            all(dim(sigma_s) == c(T_, N)),
            all(dim(sigma_eps) == c(T_, N)),
            all(dim(theta) == c(N, q_MA)))

  eps_block <- q_MA + 1L                  # current + q lags per sector
  D <- 1L + N + eps_block * N             # total state dim
  K <- 1L + N + N                         # K disturbances: c, s_1..N, eps_1..N

  # Helper: starting position of sector i's eps block in the state vector
  eps_pos_start <- function(i) 1L + N + eps_block * (i - 1L) + 1L

  # --- Transition F (D x D, sparse) ---------------------------------
  F_ <- matrix(0, D, D)
  F_[1L, 1L] <- rho                       # c_t AR(1)
  for (i in seq_len(N)) F_[1L + i, 1L + i] <- 1   # s_i RW
  for (i in seq_len(N)) {
    p <- eps_pos_start(i)
    # eps_{i,t+1} (NEW position p) gets fresh noise only — F row stays 0.
    # eps_{i,t+1-l} = OLD eps_{i,t+1-l+1} = OLD x[p + l - 1] for l=1..q
    for (l in seq_len(q_MA)) F_[p + l, p + l - 1L] <- 1
  }

  # --- Disturbance loading G (D x K) --------------------------------
  G <- matrix(0, D, K)
  G[1L, 1L] <- 1                                  # c gets disturbance 1
  for (i in seq_len(N)) G[1L + i, 1L + i] <- 1    # s_i gets disturbance 1+i
  for (i in seq_len(N)) {
    G[eps_pos_start(i), 1L + N + i] <- 1          # eps_{i,t} fresh noise
  }

  # --- Process noise covariance (K x K x T-1) ----------------------
  # NY Fed convention: eps state innovations have variance sigma_eps^2.
  # That makes the eps state's value DIRECTLY interpretable as the
  # "scaled noise" so the obs equation H entry is just theta_l (no
  # extra sigma factor), and the SV update on sigma_eps takes the eps
  # state directly as x in x_t = sigma_t * unit_eps_t.
  Sigma_eta <- array(0, c(K, K, T_ - 1L))
  for (t in seq_len(T_ - 1L)) {
    Sigma_eta[, , t] <- diag(c(
      sigma_c[t + 1L]^2,                          # c noise
      sigma_s[t + 1L, ]^2,                        # s_i noise
      sigma_eps[t + 1L, ]^2                       # eps noise (sigma^2 var)
    ))
  }

  # --- Observation H (N x D x T) -----------------------------------
  # Monthly obs at (t, i):
  #   y_{i,t} = lambda_i * c_t + s_{i,t} + eps_{i,t} + theta_1*eps_{i,t-1}
  #             + ... + theta_q*eps_{i,t-q}
  # where eps_{i,*} already has variance sigma_eps^2 (folded into state).
  H <- array(0, c(N, D, T_))
  Sigma_eps_3d <- array(0, c(N, N, T_))
  for (t in seq_len(T_)) {
    for (i in seq_len(N)) {
      H[i, 1L, t]     <- lambda[i]
      H[i, 1L + i, t] <- 1
      p <- eps_pos_start(i)
      H[i, p, t]      <- 1                         # eps_{i,t}
      for (l in seq_len(q_MA)) {
        H[i, p + l, t] <- theta[i, l]              # theta_l * eps_{i,t-l}
      }
      Sigma_eps_3d[i, i, t] <- obs_noise_floor
    }
  }

  # --- Initial state mean + covariance -----------------------------
  mu_1 <- rep(0, D)
  Sigma_1 <- diag(D)
  Sigma_1[1L, 1L] <- sigma_c[1L]^2 / max(1 - rho^2, 1e-6)
  for (i in seq_len(N)) Sigma_1[1L + i, 1L + i] <- var_s_init
  # eps states start with variance sigma_eps[1, i]^2 (matches the new
  # convention where eps state's marginal variance = sigma_eps^2).
  # Lagged eps at t=1 are unobserved priors at the same scale.
  for (i in seq_len(N)) {
    p <- eps_pos_start(i)
    for (l in 0:q_MA) Sigma_1[p + l, p + l] <- sigma_eps[1L, i]^2
  }

  list(
    H         = H,
    F         = F_,
    G         = G,
    Sigma_eps = Sigma_eps_3d,
    Sigma_eta = Sigma_eta,
    mu_1      = mu_1,
    Sigma_1   = Sigma_1,
    # Bookkeeping: positions of c, s, eps for downstream extraction
    .meta     = list(
      c_pos    = 1L,
      s_pos    = 1L + seq_len(N),
      eps_start = as.integer(vapply(seq_len(N), eps_pos_start, numeric(1L))),
      eps_block = eps_block,
      q_MA     = q_MA,
      D        = D
    )
  )
}
