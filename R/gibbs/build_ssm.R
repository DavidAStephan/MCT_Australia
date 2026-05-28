# Build the state-space-model (SSM) struct for the Australian MCT
# Variant A spec, given current parameter draws. Adapted from the
# NY Fed estimate_MCT.m SSM construction block to our simpler model.
#
# Our SSM (Variant A, monthly-only obs):
#
#   State x_t = (c_t, s_{1,t}, ..., s_{N,t})         dim (N+1)
#   y_{i,t}   = lambda[i] * c_t + s_{i,t} + eps_{i,t}, eps ~ N(0, sigma_eps[i,t]^2)
#   c_{t+1}   = rho * c_t + nu_c_{t+1},               nu_c  ~ N(0, sigma_c[t+1]^2)
#   s_{i,t+1} = s_{i,t} + nu_s_{i,t+1},               nu_si ~ N(0, sigma_s[i,t+1]^2)
#   c_1       ~ N(0, sigma_c[1]^2 / (1 - rho^2))   (stationary AR(1))
#   s_{i,1}   ~ N(0, var_s_init)                    (diffuse-ish prior)
#
# Mixed-frequency observations not yet supported in this version; will
# be added in Step 9. For v1 the model only accepts monthly obs.

#' Build the SSM list for the MCT model (non-mixed-frequency path).
#'
#' Supports Variant A (AR(1) common factor, default) and Variant B
#' (RW common factor — the canonical Stock-Watson common trend).
#'
#' @param rho       Scalar AR(1) coefficient on c. Ignored when
#'                  `variant = "B"` (transition coefficient on c is 1).
#' @param lambda    Length-N vector of time-invariant loadings.
#' @param sigma_c   Length-T vector of common-factor innovation SDs.
#' @param sigma_s   T x N matrix of sector-trend innovation SDs.
#' @param sigma_eps T x N matrix of measurement-noise SDs.
#' @param T_        Number of time periods.
#' @param N         Number of sectors.
#' @param var_s_init Prior variance on s_{i,1}. Default 100 (diffuse).
#' @param variant   "A" (AR(1) common; default) or "B" (RW common trend).
#'                  Affects F[1,1] and the c_1 prior variance: A uses
#'                  diffuse Sigma_1[1,1] = 100; B uses the Stan-Variant-B
#'                  prior c[1] ~ N(0, 2^2).
#' @return List of SSM matrices ready to pass to simulate_ssm,
#'   kalman_filter, fast_smoother, simulation_smoother.
build_mct_ssm <- function(rho, lambda, sigma_c, sigma_s, sigma_eps,
                          T_, N, var_s_init = 100,
                          variant = c("A", "B")) {
  variant <- match.arg(variant)
  rho_eff <- if (variant == "B") 1 else rho

  stopifnot(length(lambda) == N,
            length(sigma_c) == T_,
            all(dim(sigma_s) == c(T_, N)),
            all(dim(sigma_eps) == c(T_, N)))

  M <- N + 1L          # state dim
  K <- N + 1L          # disturbance dim

  # Transition F: block-diagonal-ish (rho for c, 1 for each s_i).
  # Variant A: F[1,1] = rho; Variant B: F[1,1] = 1 (RW).
  F_ <- diag(c(rho_eff, rep(1, N)))

  # Disturbance loadings G: identity (each disturbance hits its own state).
  G <- diag(M)

  # Observation H: y_{i,t} = lambda[i] c + s_i + noise. Time-invariant.
  H <- cbind(lambda, diag(N))   # N x (N+1)

  # Time-varying process-noise covariances: 3D array (M x M x T_)
  # MATLAB-style: column-major time axis on Sigma_eta is for t=1..T-1
  # (the disturbance between t and t+1). But our simulate_ssm uses
  # Sigma_eta at index t-1 for the t -> t+1 transition. We'll pass a
  # (M x M x T_-1) array; init step uses Sigma_1 separately.
  Sigma_eta <- array(0, c(M, M, T_ - 1L))
  for (t in seq_len(T_ - 1L)) {
    Sigma_eta[, , t] <- diag(c(sigma_c[t + 1L]^2,
                               sigma_s[t + 1L, ]^2))
  }

  # Time-varying observation noise: N x N x T_
  Sigma_eps <- array(0, c(N, N, T_))
  for (t in seq_len(T_)) {
    Sigma_eps[, , t] <- diag(sigma_eps[t, ]^2, nrow = N)
  }

  # Initial state.
  #   Variant A: diffuse Sigma_1[1,1] = 100 on c_1 (the stationary
  #     AR(1) variance would create a self-fulfilling bias when rho is
  #     small early in the chain; let FFBS learn the right scale).
  #   Variant B: Sigma_1[1,1] = 4 matches the Stan Variant B prior
  #     c[1] ~ N(0, 2^2) — the canonical Stock-Watson diffuse start.
  mu_1 <- rep(0, M)
  Sigma_1 <- diag(M)
  Sigma_1[1, 1] <- if (variant == "B") 4 else 100
  for (i in 1:N) Sigma_1[1 + i, 1 + i] <- var_s_init

  list(
    H         = H,
    F         = F_,
    G         = G,
    Sigma_eps = Sigma_eps,
    Sigma_eta = Sigma_eta,
    mu_1      = mu_1,
    Sigma_1   = Sigma_1
  )
}
