# R port of MCT-NYFED/functions/simulate_SSM.m
#
# Generic state-space-model simulator following Durbin & Koopman (2012):
#
#   y_t     = D_t + H_t * x_t + eps_t      (measurement)
#   x_{t+1} = C_t + F_t * x_t + G_t * eta_t (transition)
#   x_1     ~ N(mu_1, Sigma_1)
#   eps_t   ~ N(0, Sigma_eps_t)
#   eta_t   ~ N(0, Sigma_eta_t)
#
# State-space matrices passed in `ssm` may be constant or time-varying.
# Time-varying matrices use the 3rd array dim for time. Default values
# for omitted matrices mirror the MATLAB version (zeros / identity).
#
# Cross-checks against the MATLAB original (functions/simulate_SSM.m):
#   - mvnrnd in MATLAB ↔ MASS::mvrnorm in R; same parameterisation
#     (mean vector + full covariance), tiny numerical drift expected
#     in chain of rng calls but distribution is identical.
#   - Indexing is 1-based in both languages.

#' Symmetrise a matrix by averaging with its transpose. Used before
#' multivariate-normal draws to defend against numerical asymmetry in
#' covariance matrices passed through long arithmetic chains.
.symmetrize <- function(A) (A + t(A)) / 2

#' Draw a single sample from MVN(mu, Sigma). Returns a length-d numeric
#' vector. Defends against tiny non-PD asymmetry by symmetrising first.
.mvn_draw <- function(mu, Sigma) {
  d <- length(mu)
  if (d == 0L) return(numeric(0))
  s <- .symmetrize(Sigma)
  drop(MASS::mvrnorm(n = 1, mu = mu, Sigma = s))
}

#' Get the time-t slice of a state-space matrix that may be constant
#' or time-varying. For 2D matrices that vary by column (D_t, C_t) the
#' time index is the column. For 3D arrays (H_t, F_t, etc.) the time
#' index is the third array dim.
#'
#' @param X Matrix or 3D array.
#' @param t 1-based time index.
#' @param time_dim Which dimension of X is the time axis. 2 for D, C
#'   (column-varying); 3 for H, F, G, Sigma_*.
.slice_t <- function(X, t, time_dim) {
  if (time_dim == 2L) {
    # Column-major time-varying (e.g. intercepts)
    if (ncol(X) == 1L) return(X[, 1L, drop = FALSE])
    return(X[, t, drop = FALSE])
  }
  # 3D array, time_dim == 3
  if (length(dim(X)) < 3L || dim(X)[3L] == 1L) {
    # Constant — drop time dim if 3D, return as-is if 2D
    if (length(dim(X)) == 3L) return(X[, , 1L, drop = TRUE])
    return(X)
  }
  X[, , t, drop = TRUE]
}

#' Simulate from a (possibly time-varying) linear-Gaussian state-space
#' model.
#'
#' @param ssm List with fields: H (NxM or NxMxT), F (MxM or MxMx(T-1)),
#'   G (MxK or MxKx(T-1)), mu_1 (M), Sigma_1 (MxM). Optional: D (Nx1 or
#'   NxT, default zeros), C (Mx1 or Mx(T-1), default zeros), Sigma_eps
#'   (NxN or NxNxT, default zeros), Sigma_eta (KxK or KxKx(T-1), default
#'   identity).
#' @param T_ Number of periods to simulate.
#' @param need_states If TRUE, also return the latent state path.
#' @param need_disturbances If TRUE, also return the drawn disturbances.
#' @return List with fields:
#'   - Y: N x T matrix of observations (columns are y_t)
#'   - states (optional): M x T matrix of latent states (columns are x_t)
#'   - disturbances (optional): list with m_errors (N x T) and
#'     shocks (K x (T-1))
simulate_ssm <- function(ssm, T_, need_states = FALSE,
                         need_disturbances = FALSE) {
  H <- ssm$H
  F_ <- ssm$F
  G <- ssm$G

  N <- dim(H)[1L]
  M <- dim(G)[1L]
  K <- dim(G)[2L]

  # Defaults for omitted state-space matrices
  D <- if (!is.null(ssm$D)) ssm$D else matrix(0, N, 1)
  Sigma_eps <- if (!is.null(ssm$Sigma_eps)) ssm$Sigma_eps else
    array(0, c(N, N, 1))
  C <- if (!is.null(ssm$C)) ssm$C else matrix(0, M, 1)
  Sigma_eta <- if (!is.null(ssm$Sigma_eta)) ssm$Sigma_eta else
    array(diag(K), c(K, K, 1))

  # Ensure 3D arrays for H/F/G/Sigma_* (even if constant)
  ensure_3d <- function(X, expected_dims) {
    if (length(dim(X)) == 2L) {
      return(array(X, c(expected_dims, 1L)))
    }
    X
  }
  H         <- ensure_3d(H,         c(N, M))
  F_        <- ensure_3d(F_,        c(M, M))
  G         <- ensure_3d(G,         c(M, K))
  Sigma_eps <- ensure_3d(Sigma_eps, c(N, N))
  Sigma_eta <- ensure_3d(Sigma_eta, c(K, K))

  # Storage
  Y <- matrix(0, N, T_)
  states   <- if (need_states)       matrix(0, M, T_) else NULL
  m_errors <- if (need_disturbances) matrix(0, N, T_) else NULL
  shocks   <- if (need_disturbances) matrix(0, K, T_ - 1L) else NULL

  # t = 1: initial condition
  x_t <- .mvn_draw(ssm$mu_1, ssm$Sigma_1)
  H_t <- .slice_t(H, 1L, 3L)
  D_t <- .slice_t(D, 1L, 2L)[, 1L]
  Seps_t <- .slice_t(Sigma_eps, 1L, 3L)
  eps_t <- .mvn_draw(rep(0, N), Seps_t)
  Y[, 1L] <- D_t + H_t %*% x_t + eps_t

  if (need_states) {
    states[, 1L] <- x_t
    if (need_disturbances) m_errors[, 1L] <- eps_t
  }

  # Subsequent steps
  for (t in 2:T_) {
    C_t <- .slice_t(C, t - 1L, 2L)[, 1L]
    F_t <- .slice_t(F_, t - 1L, 3L)
    G_t <- .slice_t(G, t - 1L, 3L)
    Seta_t <- .slice_t(Sigma_eta, t - 1L, 3L)

    eta_t <- .mvn_draw(rep(0, K), Seta_t)
    x_t   <- C_t + F_t %*% x_t + G_t %*% eta_t
    x_t   <- as.numeric(x_t)

    D_t <- .slice_t(D, t, 2L)[, 1L]
    H_t <- .slice_t(H, t, 3L)
    Seps_t <- .slice_t(Sigma_eps, t, 3L)

    eps_t <- .mvn_draw(rep(0, N), Seps_t)
    Y[, t] <- D_t + H_t %*% x_t + eps_t

    if (need_states) {
      states[, t] <- x_t
      if (need_disturbances) {
        shocks[, t - 1L]   <- eta_t
        m_errors[, t]      <- eps_t
      }
    }
  }

  out <- list(Y = Y)
  if (need_states) out$states <- states
  if (need_disturbances) out$disturbances <- list(m_errors = m_errors,
                                                  shocks = shocks)
  out
}
