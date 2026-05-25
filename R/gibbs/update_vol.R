# R port of MCT-NYFED/functions/update_vol.m
#
# Kim-Shephard-Chib (1998) / Omori-Chib-Shephard-Nakajima (2007)
# 10-component mixture approximation to log(chi-square_1), enabling
# fast Gibbs updates of stochastic volatility processes.
#
# Model:
#   x_t           = sigma_t * eps_t,         eps ~ N(0, 1)
#   ln(sigma_t^2) = ln(sigma_{t-1}^2) + gamma * ups_t,  ups ~ N(0, 1)
#
# Given x (observed) and current sigma + gamma, return ONE draw of
# sigma from p(sigma | x, gamma) via:
#   1. Sample mixture-component indicators z_t at each t, from full
#      conditional given (x_t, sigma_t).
#   2. Conditional on indicators, y_t = log(x_t^2) - mean[z_t] is a
#      Gaussian observation of the RW process ln(sigma_t^2) with known
#      observation variance vars[z_t] and process variance gamma^2.
#   3. Run scalar Kalman filter + backward sample (FFBS) on this
#      Gaussian SSM → draw of ln(sigma_t^2) path.
#   4. Return sigma_t = exp(ln(sigma_t^2) / 2).
#
# Reference: Omori, Chib, Shephard, Nakajima (2007), J. Econometrics 140.

# Mixture-approximation constants — copied verbatim from
# MCT-NYFED/functions/update_vol.m. **Do not edit** — these are
# calibrated to the log(chi-square_1) distribution.
.KSC_PROBS <- c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,
                0.18842, 0.12047, 0.05591, 0.01575,  0.00115)
.KSC_MEANS <- c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173,
                -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
.KSC_VARS  <- c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,
                0.98583, 1.57469, 2.54498, 4.16591,  7.33342)
.KSC_STDVS <- sqrt(.KSC_VARS)
.KSC_BARRIER <- 1e-3   # added to x^2 to handle x == 0 (matches MATLAB)
.KSC_SMALL   <- 1e-6   # numerical guard in backward pass

#' Single Gibbs update of an SV process via the Kim-Shephard mixture.
#'
#' @param x        Length-T vector of observed pre-volatility errors.
#'                 May contain NA for missing data.
#' @param sigma    Length-T vector of CURRENT volatility estimates
#'                 (input to the update).
#' @param gamma    Scalar RW step size on ln(sigma^2).
#' @param var_prior Prior variance on ln(sigma_1^2). Default 1e2.
#' @param mean_prior Prior mean of sigma_1 (used to anchor the filter).
#'                   Default 1.
#' @param AR_prior Persistence in the log-vol AR(1) — for pure RW use 1.
#' @return Length-T vector of NEW sigma values (one draw from the
#'   posterior under the mixture approximation).
update_vol <- function(x, sigma, gamma,
                       var_prior = 1e2,
                       mean_prior = 1,
                       AR_prior = 1) {
  x     <- as.numeric(x)
  sigma <- as.numeric(sigma)
  T_    <- length(x)
  n_comp <- length(.KSC_PROBS)

  # log(eps_t^2) under the current sigma estimate.
  ln_eps <- log(x^2 + .KSC_BARRIER) - log(sigma^2)

  # Posterior probability of each mixture component given the observed
  # log(eps_t^2). Vectorised across t.
  ln_eps_rep <- matrix(ln_eps, T_, n_comp, byrow = FALSE)
  means_rep  <- matrix(.KSC_MEANS, T_, n_comp, byrow = TRUE)
  stdvs_rep  <- matrix(.KSC_STDVS, T_, n_comp, byrow = TRUE)
  probs_rep  <- matrix(.KSC_PROBS, T_, n_comp, byrow = TRUE)

  likelihood <- exp(-0.5 * ((ln_eps_rep - means_rep) / stdvs_rep)^2) /
                stdvs_rep
  pxlikelihood <- likelihood * probs_rep
  xmlikelihood <- rowSums(pxlikelihood)
  posteriors   <- pxlikelihood / xmlikelihood

  # Missing data: fall back to prior weights (matches MATLAB).
  na_rows <- !is.finite(rowSums(posteriors))
  if (any(na_rows)) {
    posteriors[na_rows, ] <- probs_rep[na_rows, ]
  }

  # Draw indicator for each t. rmultinom(n=1, size=1, prob=p) returns
  # a 0/1 column vector with exactly one 1 in the sampled position.
  # We use that column as a row of `weights`.
  weights <- matrix(0, T_, n_comp)
  for (t in seq_len(T_)) {
    weights[t, ] <- rmultinom(1L, 1L, posteriors[t, ])
  }

  # Build the Gaussian SSM for log(sigma_t^2):
  #   y_t           = ln(x_t^2) - mean[z_t]
  #   ln(sigma_t^2) = AR_prior * ln(sigma_{t-1}^2) + gamma * ups_t
  #   observation noise variance = vars[z_t]
  gamsq  <- gamma^2
  ln_x   <- log(x^2 + .KSC_BARRIER)
  mean_t <- as.numeric(weights %*% .KSC_MEANS)
  vars_t <- as.numeric(weights %*% .KSC_VARS)   # variance, not stdev
  y_t    <- ln_x - mean_t

  # Allocate forward-filter storage. Indexing matches MATLAB (T+1
  # entries; t=1 is the prior anchor, t=2..T+1 are the post-update
  # filtered states).
  x1_KF <- numeric(T_ + 1L)
  p1_KF <- numeric(T_ + 1L)
  x2_KF <- numeric(T_ + 1L)
  p2_KF <- numeric(T_ + 1L)

  x1_KF[1L] <- 2 * log(mean_prior)
  p1_KF[1L] <- var_prior
  x1 <- x1_KF[1L]
  p1 <- p1_KF[1L]

  # --- Forward pass (with missing-data fallback) ---------------------
  any_missing <- any(is.na(y_t))
  for (t in seq_len(T_)) {
    # Predict
    x2 <- AR_prior * x1
    p2 <- AR_prior^2 * p1 + gamsq

    # Update (skip if y_t missing)
    if (!any_missing || !is.na(y_t[t])) {
      h <- p2 + vars_t[t]
      k <- p2 / h
      x1 <- x2 + k * (y_t[t] - x2)
      p1 <- p2 - k * p2
    } else {
      x1 <- x2
      p1 <- p2
    }

    x1_KF[t + 1L] <- x1
    p1_KF[t + 1L] <- p1
    x2_KF[t + 1L] <- x2
    p2_KF[t + 1L] <- p2
  }

  # --- Backward sample (FFBS) ----------------------------------------
  utmp <- rnorm(T_ + 1L)
  ln_sigmasq <- numeric(T_ + 1L)

  x3mean <- x1
  p3     <- p1
  x3     <- x3mean + sqrt(p3) * utmp[T_ + 1L]
  ln_sigmasq[T_ + 1L] <- x3

  for (t in T_:1L) {
    x1 <- x1_KF[t]
    p1 <- p1_KF[t]
    x2 <- x2_KF[t + 1L]
    p2 <- p2_KF[t + 1L]
    if (p2 > .KSC_SMALL) {
      p2i    <- 1 / p2
      k      <- AR_prior * p1 * p2i
      x3mean <- x1 + k * (x3 - x2)
      p3     <- p1 - AR_prior * k * p1
    } else {
      x3mean <- x1
      p3     <- p1
    }
    # Numerical guard: p3 can drift slightly negative from floating
    # point. Clip to 0 (matches MATLAB sqrt() behaviour on tiny negs
    # which would NaN; this is safer).
    if (p3 < 0) p3 <- 0
    x3 <- x3mean + sqrt(p3) * utmp[t]
    ln_sigmasq[t] <- x3
  }

  # Return updated sigma_t = exp(ln_sigmasq_t / 2) for t = 1..T.
  # MATLAB skips ln_sigmasq[1] (the prior anchor) and returns 2..T+1.
  exp(ln_sigmasq[2:(T_ + 1L)] / 2)
}
