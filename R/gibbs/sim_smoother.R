# R port of MCT-NYFED/functions/simulation_smoother.m
#
# Durbin & Koopman (2002) simulation smoother. Draws ONE realization of
# the latent state path conditional on the data: x ~ p(x | y, ssm).
#
# Algorithm:
#   1. Build aux SSM with intercepts and initial mean zeroed out.
#   2. Simulate (y_sim, x_sim) ~ p(y, x | aux_ssm).
#   3. Run fast_smoother on (y - y_sim) under the ORIGINAL ssm,
#      giving smoothed states x_smooth.
#   4. Return x = x_sim + x_smooth. This is a valid draw from p(x | y).
#
# This is the workhorse for the Gibbs sampler — the FFBS step that
# updates the latent factor and sector trends each iteration.

#' Draw one realization of the latent state path from p(x | y).
#'
#' @param Y N x T matrix of observations (NA for missing).
#' @param ssm State-space matrices (see simulate_ssm() / kalman_filter()).
#' @param need_disturbances If TRUE, also draw the disturbances.
#' @return List with:
#'   - states: M x T matrix, a draw from p(x | y)
#'   - disturbances (if requested): list of m_errors + shocks draws
simulation_smoother <- function(Y, ssm, need_disturbances = FALSE) {
  T_ <- ncol(Y)
  M  <- dim(ssm$G)[1L]

  # Build aux SSM: zero out intercepts and initial mean.
  ssm_aux <- ssm
  ssm_aux$mu_1 <- rep(0, M)
  if (!is.null(ssm$C)) ssm_aux$C <- matrix(0, nrow(ssm$C), ncol(ssm$C))
  if (!is.null(ssm$D)) ssm_aux$D <- matrix(0, nrow(ssm$D), ncol(ssm$D))

  # 1. Simulate from aux SSM
  sim <- simulate_ssm(ssm_aux, T_ = T_,
                     need_states = TRUE,
                     need_disturbances = need_disturbances)

  # 2. Smooth (Y - Y_sim) under the ORIGINAL ssm. Note: the residual
  #    is computed in OBSERVATION-VALUE space (i.e. including D); the
  #    smoother re-subtracts D internally via Kalman filter, but here
  #    Y_sim was generated from D = 0, so the residual Y - Y_sim still
  #    has the proper D offset that the original ssm's KF will remove.
  smooth <- fast_smoother(Y - sim$Y, ssm,
                          need_states = TRUE)

  # 3. Combine
  out <- list(states = sim$states + smooth$states)
  if (need_disturbances) {
    out$disturbances <- list(
      m_errors = sim$disturbances$m_errors + smooth$disturbances$m_errors,
      shocks   = sim$disturbances$shocks   + smooth$disturbances$shocks
    )
  }
  out
}
