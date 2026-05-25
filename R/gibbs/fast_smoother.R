# R port of MCT-NYFED/functions/fast_smoother.m
#
# Disturbance + state smoother following Durbin & Koopman (2012). Given
# the filtered prediction outputs from kalman_filter(), runs a backward
# recursion to compute smoothed disturbances E[eta_t|y_{1:T}] and
# E[eps_t|y_{1:T}], then a forward recursion to reconstruct smoothed
# states E[x_t|y_{1:T}].
#
# Notation (matches the MATLAB original):
#   e_t      = prediction error at t (from KF)
#   S_t^{-1} = inverse innovation variance at t (from KF, sym + PD)
#   K_t      = Kalman gain at t (from KF)
#   u_t = S_t^{-1} e_t - (F_t K_t)' r_t   ("scaled smoothing error")
#   r_{t-1} = H_t' u_t + F_t' r_t          (backward recursion)
#
# r_T initialised to zero. After the backward pass:
#   smoothed shock at t-1:  eta_hat_{t-1} = Sigma_eta_{t-1} G_{t-1}' r_t
#   smoothed obs noise at t: eps_hat_t    = Sigma_eps_t S_t^{-1} u_t
#   smoothed state at 1:    x_hat_1       = mu_1 + Sigma_1 r_1
#   forward propagation:    x_hat_t       = C_{t-1} + F_{t-1} x_hat_{t-1}
#                                            + G_{t-1} eta_hat_{t-1}
#
# **NOT** ported in this pass: MSE recursions (smoothed variances). Add
# if/when we need uncertainty bands on smoothed states (Gibbs sampler
# doesn't need them — it just uses point smooths).

#' Fast disturbance + state smoother.
#'
#' @param Y N x T matrix of observations (NA for missing).
#' @param ssm State-space matrices (see simulate_ssm() / kalman_filter()).
#' @param need_states If TRUE, also run the forward state recursion.
#' @return List with:
#'   - disturbances: list of m_errors (N x T) and shocks (K x (T-1))
#'   - states (if need_states): M x T matrix of E[x_t | y_{1:T}]
fast_smoother <- function(Y, ssm, need_states = FALSE) {
  H  <- ssm$H
  F_ <- ssm$F
  G  <- ssm$G
  N  <- nrow(Y)
  T_ <- ncol(Y)
  M  <- dim(G)[1L]
  K  <- dim(G)[2L]

  Sigma_eps <- if (!is.null(ssm$Sigma_eps)) ssm$Sigma_eps else
    array(0, c(N, N, 1))
  Sigma_eta <- if (!is.null(ssm$Sigma_eta)) ssm$Sigma_eta else
    array(diag(K), c(K, K, 1))

  ensure_3d <- function(X, dims) {
    if (length(dim(X)) == 2L) array(X, c(dims, 1L)) else X
  }
  H         <- ensure_3d(H,         c(N, M))
  F_        <- ensure_3d(F_,        c(M, M))
  G         <- ensure_3d(G,         c(M, K))
  Sigma_eps <- ensure_3d(Sigma_eps, c(N, N))
  Sigma_eta <- ensure_3d(Sigma_eta, c(K, K))

  slice2 <- function(X, t) {
    matrix(X[, , t], nrow = dim(X)[1L], ncol = dim(X)[2L])
  }

  isTV_H         <- dim(H)[3L]         > 1L
  isTV_Sigma_eps <- dim(Sigma_eps)[3L] > 1L
  isTV_F         <- dim(F_)[3L]        > 1L
  isTV_G         <- dim(G)[3L]         > 1L
  isTV_Sigma_eta <- dim(Sigma_eta)[3L] > 1L

  # Run KF with need_prediction = TRUE to get e, S^-1, K at every t.
  kf <- kalman_filter(Y, ssm, need_prediction = TRUE)
  e     <- kf$prediction$error      # N x T
  S_inv <- kf$prediction$invMSE     # N x N x T
  Kt    <- kf$prediction$gain       # M x N x T

  # Storage
  m_errors <- matrix(0, N, T_)
  shocks   <- matrix(0, K, T_ - 1L)

  # --- Backward recursion (initialise at t = T) ----------------------
  Ht          <- if (isTV_H) slice2(H, T_) else slice2(H, 1L)
  Sigma_eps_t <- if (isTV_Sigma_eps) slice2(Sigma_eps, T_) else
    slice2(Sigma_eps, 1L)

  nonmiss <- !is.na(e[, T_])
  if (any(nonmiss)) {
    e_aux           <- e[nonmiss, T_]
    S_inv_aux       <- S_inv[nonmiss, nonmiss, T_, drop = FALSE]
    S_inv_aux       <- matrix(S_inv_aux,
                              nrow = sum(nonmiss),
                              ncol = sum(nonmiss))
    Ht_aux          <- Ht[nonmiss, , drop = FALSE]
    Sigma_eps_t_aux <- Sigma_eps_t[, nonmiss, drop = FALSE]
    u_aux           <- as.numeric(S_inv_aux %*% e_aux)
    r_aux           <- as.numeric(t(Ht_aux) %*% u_aux)
    m_errors[, T_]  <- as.numeric(Sigma_eps_t_aux %*% u_aux)
  } else {
    r_aux           <- rep(0, M)
    # m_errors[, T_] stays at 0
  }

  for (t in (T_ - 1L):1L) {
    if (isTV_F)         Ft          <- slice2(F_, t)         else if (t == T_ - 1L) Ft <- slice2(F_, 1L)
    if (isTV_G)         Gt          <- slice2(G,  t)         else if (t == T_ - 1L) Gt <- slice2(G, 1L)
    if (isTV_Sigma_eta) Sigma_eta_t <- slice2(Sigma_eta, t)  else if (t == T_ - 1L) Sigma_eta_t <- slice2(Sigma_eta, 1L)

    # Smoothed shock for the transition INTO time t+1 (i.e. eta_t in
    # MATLAB indexing where shocks(:, t) is the shock between t and t+1)
    shocks[, t] <- as.numeric(Sigma_eta_t %*% t(Gt) %*% r_aux)

    if (isTV_H)         Ht          <- slice2(H, t)          else if (t == T_ - 1L) Ht <- slice2(H, 1L)
    if (isTV_Sigma_eps) Sigma_eps_t <- slice2(Sigma_eps, t)  else if (t == T_ - 1L) Sigma_eps_t <- slice2(Sigma_eps, 1L)

    nonmiss <- !is.na(e[, t])
    if (any(nonmiss)) {
      e_aux           <- e[nonmiss, t]
      S_inv_aux       <- S_inv[nonmiss, nonmiss, t, drop = FALSE]
      S_inv_aux       <- matrix(S_inv_aux,
                                nrow = sum(nonmiss),
                                ncol = sum(nonmiss))
      Kt_aux          <- Kt[, nonmiss, t, drop = FALSE]
      Kt_aux          <- matrix(Kt_aux, nrow = M, ncol = sum(nonmiss))
      Ht_aux          <- Ht[nonmiss, , drop = FALSE]
      Sigma_eps_t_aux <- Sigma_eps_t[, nonmiss, drop = FALSE]

      u_aux <- as.numeric(S_inv_aux %*% e_aux -
                          t(Ft %*% Kt_aux) %*% r_aux)
      r_aux <- as.numeric(t(Ht_aux) %*% u_aux + t(Ft) %*% r_aux)
      m_errors[, t] <- as.numeric(Sigma_eps_t_aux %*% u_aux)
    } else {
      # No obs at t: u_t = -(Ft Kt)' r_t but Kt has 0 cols, so u_t = 0.
      # Then r_{t-1} = Ft' r_t (state-space backward recursion only).
      r_aux <- as.numeric(t(Ft) %*% r_aux)
      # m_errors[, t] stays at 0.
    }
  }

  out <- list(disturbances = list(m_errors = m_errors, shocks = shocks))

  # --- Optional forward state recursion ------------------------------
  if (need_states) {
    C  <- if (!is.null(ssm$C)) ssm$C else matrix(0, M, 1)
    isTV_C <- ncol(C) > 1L
    states <- matrix(0, M, T_)

    Ct <- C[, 1L]
    x  <- as.numeric(ssm$mu_1 + ssm$Sigma_1 %*% r_aux)
    states[, 1L] <- x

    for (t in 2:T_) {
      if (isTV_C) Ct <- C[, t - 1L]
      if (isTV_F) Ft <- slice2(F_, t - 1L)
      if (isTV_G) Gt <- slice2(G,  t - 1L)

      x <- as.numeric(Ct + Ft %*% x + Gt %*% shocks[, t - 1L])
      states[, t] <- x
    }
    out$states <- states
  }

  out
}
