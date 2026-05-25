# R port of MCT-NYFED/functions/Kalman_filter.m
#
# Standard Kalman filter for the state-space model:
#
#   y_t     = D_t + H_t * x_t + eps_t,   eps_t ~ N(0, Sigma_eps_t)
#   x_{t+1} = C_t + F_t * x_t + G_t * eta_t,  eta_t ~ N(0, Sigma_eta_t)
#   x_1     ~ N(mu_1, Sigma_1)
#
# Handles missing data via NA-coded entries in Y (per-column missing
# pattern). All state-space matrices may be constant or time-varying.
#
# Cross-checks (Phase 3 from gibbs_port_plan.md):
#   - log-likelihood should agree with the Akf Stan filter at
#     `stan/mct_aus_Akf.stan` on the same SSM and data — within
#     numerical tolerance.
#   - Sequence of filtered means should match what's computed inside
#     the simulation_smoother when seeded with the same inputs.

#' Kalman filter recursion. Returns log-likelihood and optionally the
#' one-step prediction / filtered state objects needed by the smoother.
#'
#' @param Y  N x T matrix; NA entries are treated as missing.
#' @param ssm List of state-space matrices (see simulate_ssm() for shape).
#' @param need_prediction If TRUE, also return per-t prediction error,
#'   inverse-MSE, and Kalman gain matrices. Required for the smoother.
#' @param need_filter If TRUE (requires need_prediction), also return
#'   the filtered state means + covariances at each t.
#' @return List with `log_likelihood` (scalar) and optionally
#'   `prediction` and `filter` sublists matching the MATLAB output.
kalman_filter <- function(Y, ssm,
                          need_prediction = FALSE,
                          need_filter     = FALSE) {
  stopifnot(!need_filter || need_prediction)

  H  <- ssm$H
  F_ <- ssm$F
  G  <- ssm$G
  N  <- nrow(Y)
  T_ <- ncol(Y)
  M  <- dim(G)[1L]
  K  <- dim(G)[2L]

  # Subtract D from Y (per MATLAB; D defaults to 0)
  if (!is.null(ssm$D)) {
    Y <- Y - if (ncol(ssm$D) == 1L) {
      matrix(ssm$D, N, T_)
    } else {
      ssm$D
    }
  }

  Sigma_eps <- if (!is.null(ssm$Sigma_eps)) ssm$Sigma_eps else
    array(0, c(N, N, 1))
  C  <- if (!is.null(ssm$C)) ssm$C else matrix(0, M, 1)
  Sigma_eta <- if (!is.null(ssm$Sigma_eta)) ssm$Sigma_eta else
    array(diag(K), c(K, K, 1))

  # Ensure 3D shape for matrices that may be supplied as 2D
  ensure_3d <- function(X, dims) {
    if (length(dim(X)) == 2L) array(X, c(dims, 1L)) else X
  }
  H         <- ensure_3d(H,         c(N, M))
  F_        <- ensure_3d(F_,        c(M, M))
  G         <- ensure_3d(G,         c(M, K))
  Sigma_eps <- ensure_3d(Sigma_eps, c(N, N))
  Sigma_eta <- ensure_3d(Sigma_eta, c(K, K))

  isTV_H         <- dim(H)[3L]         > 1L
  isTV_Sigma_eps <- dim(Sigma_eps)[3L] > 1L
  isTV_C         <- ncol(C)            > 1L
  isTV_F         <- dim(F_)[3L]        > 1L
  isTV_G         <- dim(G)[3L]         > 1L
  isTV_Sigma_eta <- dim(Sigma_eta)[3L] > 1L

  # Initialise accumulator. Note: the MATLAB original uses a fixed
  # `-(T*N/2)*log(2*pi)` term that ignores missing data. We instead
  # accumulate `-0.5 * n_obs_t * log(2*pi)` per step so the returned
  # log-likelihood is the joint log-density of OBSERVED data only
  # (the standard convention; matches dmvnorm of the observed entries).
  log_lik <- 0

  prediction_error  <- if (need_prediction)
    matrix(NA_real_, N, T_) else NULL
  prediction_invMSE <- if (need_prediction)
    array(NA_real_, c(N, N, T_)) else NULL
  prediction_gain   <- if (need_prediction)
    array(NA_real_, c(M, N, T_)) else NULL
  filter_mu         <- if (need_filter) matrix(0, M, T_) else NULL
  filter_Sigma      <- if (need_filter) array(0, c(M, M, T_)) else NULL

  # Helper: slice a 3D array at time t, preserving the 2D matrix
  # shape (R drops 1x1 to scalar otherwise).
  slice2 <- function(X, t) {
    matrix(X[, , t], nrow = dim(X)[1L], ncol = dim(X)[2L])
  }

  # Auxiliary slices for t = 1
  Ct          <- C[, 1L]
  Ft          <- slice2(F_, 1L)
  Gt          <- slice2(G,  1L)
  Sigma_eta_t <- slice2(Sigma_eta, 1L)
  mu          <- ssm$mu_1
  Sigma       <- .symmetrize(ssm$Sigma_1)
  Ht          <- slice2(H, 1L)
  Sigma_eps_t <- slice2(Sigma_eps, 1L)

  # Helper: symmetric PD-matrix inverse via Cholesky (matches MATLAB
  # linsolve with SYM/POSDEF options for numerical stability).
  sym_inv <- function(A) chol2inv(chol(A))

  # First iteration (t = 1)
  nonmiss <- !is.na(Y[, 1L])
  if (any(nonmiss)) {
    Y_aux   <- Y[nonmiss, 1L]
    Ht_aux  <- Ht[nonmiss, , drop = FALSE]
    Sep_aux <- Sigma_eps_t[nonmiss, nonmiss, drop = FALSE]
    e       <- as.numeric(Y_aux - Ht_aux %*% mu)
    S       <- .symmetrize(Ht_aux %*% Sigma %*% t(Ht_aux) + Sep_aux)
    S_inv   <- sym_inv(S)
    Kt      <- Sigma %*% t(Ht_aux) %*% S_inv
    log_lik <- log_lik -
      0.5 * (length(e) * log(2 * pi) + log(det(S)) +
             as.numeric(t(e) %*% S_inv %*% e))
  } else {
    # All missing at t=1: no likelihood contribution, Kt is N(M, 0)-shaped 0
    Ht_aux <- Ht[integer(0), , drop = FALSE]
    e      <- numeric(0)
    Kt     <- matrix(0, M, 0)
  }

  if (need_prediction) {
    if (any(nonmiss)) {
      prediction_error[nonmiss, 1L]            <- e
      prediction_invMSE[nonmiss, nonmiss, 1L]  <- S_inv
      prediction_gain[, nonmiss, 1L]           <- Kt
    }
    if (need_filter) {
      filter_mu[, 1L]      <- mu
      filter_Sigma[, , 1L] <- Sigma
    }
  }

  # Subsequent iterations
  for (t in 2:T_) {
    if (isTV_C)         Ct          <- C[, t - 1L]
    if (isTV_F)         Ft          <- slice2(F_, t - 1L)
    if (isTV_G)         Gt          <- slice2(G,  t - 1L)
    if (isTV_Sigma_eta) Sigma_eta_t <- slice2(Sigma_eta, t - 1L)

    # Filtered → predicted for the next step
    # mu_pred    = Ct + Ft * (mu_filt + Kt*e)
    # Sigma_pred = Ft * (Sigma_filt - Kt*Ht_aux*Sigma_filt') * Ft' + Gt*Sigma_eta*Gt'
    if (length(e) > 0L) {
      mu <- as.numeric(Ct + Ft %*% (mu + Kt %*% e))
      Sigma_filtered_correction <- Kt %*% Ht_aux %*% t(Sigma)
    } else {
      mu <- as.numeric(Ct + Ft %*% mu)
      Sigma_filtered_correction <- matrix(0, M, M)
    }
    Sigma <- .symmetrize(
      Ft %*% (Sigma - Sigma_filtered_correction) %*% t(Ft) +
        Gt %*% Sigma_eta_t %*% t(Gt)
    )

    if (isTV_H)         Ht          <- slice2(H, t)
    if (isTV_Sigma_eps) Sigma_eps_t <- slice2(Sigma_eps, t)

    nonmiss <- !is.na(Y[, t])
    if (any(nonmiss)) {
      Y_aux   <- Y[nonmiss, t]
      Ht_aux  <- Ht[nonmiss, , drop = FALSE]
      Sep_aux <- Sigma_eps_t[nonmiss, nonmiss, drop = FALSE]
      e       <- as.numeric(Y_aux - Ht_aux %*% mu)
      S       <- .symmetrize(Ht_aux %*% Sigma %*% t(Ht_aux) + Sep_aux)
      S_inv   <- sym_inv(S)
      Kt      <- Sigma %*% t(Ht_aux) %*% S_inv
      log_lik <- log_lik -
      0.5 * (length(e) * log(2 * pi) + log(det(S)) +
             as.numeric(t(e) %*% S_inv %*% e))
    } else {
      Ht_aux <- Ht[integer(0), , drop = FALSE]
      e      <- numeric(0)
      Kt     <- matrix(0, M, 0)
    }

    if (need_prediction) {
      if (any(nonmiss)) {
        prediction_error[nonmiss, t]           <- e
        prediction_invMSE[nonmiss, nonmiss, t] <- S_inv
        prediction_gain[, nonmiss, t]          <- Kt
      }
      if (need_filter) {
        filter_mu[, t]      <- mu
        filter_Sigma[, , t] <- Sigma
      }
    }
  }

  out <- list(log_likelihood = log_lik)
  if (need_prediction) {
    out$prediction <- list(error  = prediction_error,
                           invMSE = prediction_invMSE,
                           gain   = prediction_gain)
    if (need_filter) {
      out$filter <- list(mu = filter_mu, Sigma = filter_Sigma)
    }
  }
  out
}
