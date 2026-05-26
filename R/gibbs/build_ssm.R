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

#' Build the SSM list for our Variant A MCT model.
#'
#' @param rho       Scalar AR(1) coefficient on c.
#' @param lambda    Length-N vector of time-invariant loadings.
#' @param sigma_c   Length-T vector of common-factor innovation SDs.
#' @param sigma_s   T x N matrix of sector-trend innovation SDs.
#' @param sigma_eps T x N matrix of measurement-noise SDs.
#' @param T_        Number of time periods.
#' @param N         Number of sectors.
#' @param var_s_init Prior variance on s_{i,1}. Default 100 (diffuse).
#' @return List of SSM matrices ready to pass to simulate_ssm,
#'   kalman_filter, fast_smoother, simulation_smoother.
build_mct_ssm <- function(rho, lambda, sigma_c, sigma_s, sigma_eps,
                          T_, N, var_s_init = 100) {
  stopifnot(length(lambda) == N,
            length(sigma_c) == T_,
            all(dim(sigma_s) == c(T_, N)),
            all(dim(sigma_eps) == c(T_, N)))

  M <- N + 1L          # state dim
  K <- N + 1L          # disturbance dim

  # Transition F: block-diagonal-ish (rho for c, 1 for each s_i).
  F_ <- diag(c(rho, rep(1, N)))

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

  # Initial state. Originally we used the stationary AR(1) variance
  # for c_1: Sigma_1[1,1] = sigma_c^2 / (1 - rho^2). That created a
  # self-fulfilling bias: when rho is small early in the chain, c_1
  # is forced near 0 (small stationary variance), which biases c-draws
  # small, which biases rho small. Use a diffuse prior on c_1 instead;
  # the FFBS will learn the right scale from the data.
  mu_1 <- rep(0, M)
  Sigma_1 <- diag(M)
  Sigma_1[1, 1] <- 100      # diffuse on c_1
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
