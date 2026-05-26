# R port of MCT-NYFED/functions/update_theta.m
#
# Conjugate update for MA coefficients in:
#
#   x_t = (1 + theta_1 L + ... + theta_q L^q) * sigma_t * eps_t,
#   eps_t ~ N(0, 1)
#
# Given the data x (residuals "whitened" by sigma) and lagged
# residuals as regressors, the posterior of theta is Gaussian under
# a Normal prior with diagonal precision P[k, k] = prec_prior * k^2
# (the k^2 weighting damps higher-order MA terms).
#
# Posterior precision: P + X'X
# Posterior mean: (P + X'X)^{-1} X'y
#
# After sampling theta from this Gaussian, INVERTIBILITY is enforced
# by reflecting any roots of the MA polynomial that lie outside the
# unit circle to their reciprocals. This is a deterministic
# transformation, not strictly a posterior draw, but it keeps the
# subsequent identification clean (without invertibility, the MA
# representation is non-unique — multiple thetas give the same
# autocovariance structure).

#' Update MA coefficients via conjugate Normal posterior + invertibility.
#'
#' @param y          Length-T vector of innovations (the LHS of the
#'                   MA regression: u_t = eps_t + sum_l theta_l eps_{t-l}).
#'                   Effectively: u_t for which we infer the theta that
#'                   generated it.
#' @param x          T x q matrix of lagged eps regressors, where
#'                   x[t, l] = eps_{t - l}.
#' @param prec_prior Positive scalar prior precision on theta_1.
#'                   Higher-lag thetas get tightened by k^2.
#' @return Length-q numeric vector of new theta draws.
update_theta <- function(y, x, prec_prior) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  q_MA <- ncol(x)

  # Prior precision matrix (diagonal, k^2 weighted)
  P <- diag(prec_prior * (seq_len(q_MA))^2, nrow = q_MA)

  # Posterior precision and mean
  Pinv_post <- tryCatch(
    chol2inv(chol(.symmetrize(P + t(x) %*% x))),
    error = function(e) MASS::ginv(.symmetrize(P + t(x) %*% x))
  )
  m_post <- as.numeric(Pinv_post %*% (t(x) %*% y))

  # Sample theta from N(m_post, Pinv_post)
  theta <- as.numeric(MASS::mvrnorm(1, m_post, .symmetrize(Pinv_post)))

  # Enforce invertibility: reflect any MA roots outside unit circle
  # to their reciprocals. The MA polynomial in lag operator L is
  #   1 + theta_1 L + theta_2 L^2 + ... + theta_q L^q
  # ⇔ in standard polynomial form (highest power first):
  #   theta_q * z^q + theta_{q-1} * z^{q-1} + ... + theta_1 * z + 1
  # i.e. c(rev(theta), 1) — matches MATLAB's `roots([flip(theta) 1])`.
  poly_coefs <- c(rev(theta), 1)
  ma_roots <- polyroot(poly_coefs)
  bad <- which(Mod(ma_roots) > 1)
  if (length(bad) > 0L) {
    ma_roots[bad] <- 1 / ma_roots[bad]
    # Rebuild polynomial from reflected roots.
    # R doesn't have a direct equivalent of MATLAB's poly(); construct
    # via cumulative multiplication of factors (z - r_k).
    p <- 1
    for (r in ma_roots) p <- c(p, 0) - c(0, r * p)
    # p is now (highest-power-first) coefficient vector of length q+1.
    # Normalise so the leading coef is 1 (matches MA polynomial form
    # where the constant term is 1).
    p <- p / p[length(p)]
    # The MA coefficients are p[1:q] (excluding the constant 1 at end).
    theta_new <- Re(rev(p[seq_len(q_MA)]))
    # Pad if `roots` dropped leading zeros (matches MATLAB safeguard)
    if (length(theta_new) < q_MA) {
      theta_new <- c(theta_new, rep(0, q_MA - length(theta_new)))
    }
    theta <- theta_new
  }
  theta
}
