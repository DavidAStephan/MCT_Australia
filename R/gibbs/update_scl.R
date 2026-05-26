# R port of MCT-NYFED/functions/update_scl.m
#
# Discrete scale-mixture update for the outlier model:
#
#   x_t = s_t * eps_t,  eps_t ~ N(0, 1)
#   s_t ~ Categorical(vals, probs)
#
# Used to allow heavy-tailed measurement noise without changing
# anything else in the SSM. At each t, sample s_t from its full
# conditional posterior P(s_t = vals[k] | x_t) ∝ probs[k] * phi(x_t; vals[k]).
#
# Typical NY Fed setup:
#   vals = c(1, seq(2, 10, length.out = 39))   (40-point grid)
#   probs = c(ps, (1 - ps) / 39, ..., (1 - ps) / 39)
#
# So `vals[1] = 1` is the "normal observation" mass; the larger values
# capture outliers. The probability mass on `vals[1]` is `ps` (close
# to 1); the rest is split equally among the larger scales.

#' Sample scale indicator at each t from its full conditional posterior.
#'
#' @param x     Length-T vector of "whitened" residuals
#'              (residuals / sigma_eps; should be ~N(0, s_t^2)).
#'              NA entries get assigned the prior in lieu of posterior.
#' @param vals  Length-n_s vector of candidate scale values
#'              (smallest typically 1; larger values = outlier scales).
#' @param probs Length-n_s vector of prior mixture probabilities,
#'              summing to 1.
#' @return Length-T vector of sampled scale values (from `vals`).
update_scl <- function(x, vals, probs) {
  x     <- as.numeric(x)
  vals  <- as.numeric(vals)
  probs <- as.numeric(probs)
  T_    <- length(x)
  n_s   <- length(probs)
  stopifnot(length(vals) == n_s,
            abs(sum(probs) - 1) < 1e-6)

  # Gaussian likelihood at each (t, k): phi(x_t / vals[k]) / vals[k]
  x_rep    <- matrix(x, T_, n_s, byrow = FALSE)
  vals_rep <- matrix(vals, T_, n_s, byrow = TRUE)
  probs_rep <- matrix(probs, T_, n_s, byrow = TRUE)

  likelihood   <- exp(-0.5 * (x_rep / vals_rep)^2) / vals_rep
  joint        <- likelihood * probs_rep
  marg         <- rowSums(joint)
  posteriors   <- joint / matrix(marg, T_, n_s)

  # Missing data: use prior weights
  na_rows <- !is.finite(rowSums(posteriors))
  if (any(na_rows)) posteriors[na_rows, ] <- probs_rep[na_rows, ]

  # Sample one indicator per t from the posterior categorical
  s_out <- numeric(T_)
  for (t in seq_len(T_)) {
    k <- sample.int(n_s, size = 1L, prob = posteriors[t, ])
    s_out[t] <- vals[k]
  }
  s_out
}
