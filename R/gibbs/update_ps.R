# R port of MCT-NYFED/functions/update_ps.m
#
# Beta-conjugate update for the probability of "normal" (non-outlier)
# observation in the discrete scale mixture:
#
#   indicator_t = 1 if scale_t == vals[1] (i.e. normal obs)
#               = 0 otherwise (outlier)
#   indicator_t ~ Bernoulli(ps)
#   ps          ~ Beta(a_prior, b_prior)
#
# Posterior: Beta(a_prior + n_normal, b_prior + n_outlier).
#
# Vectorised over N sectors when x is T x N — returns N posterior draws.

#' Draw new ps from Beta posterior.
#'
#' @param x        T x N binary indicator matrix (1 = scale equals
#'                 vals[1] / "normal", 0 = outlier).
#' @param a_prior  Length-N (or scalar) prior alpha.
#' @param b_prior  Length-N (or scalar) prior beta.
#' @return         Length-N vector of new ps draws.
update_ps <- function(x, a_prior, b_prior) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  N <- ncol(x)
  a_prior <- rep_len(as.numeric(a_prior), N)
  b_prior <- rep_len(as.numeric(b_prior), N)
  n_one   <- colSums(x == 1L, na.rm = TRUE)
  n_zero  <- colSums(x != 1L, na.rm = TRUE)
  rbeta(N, a_prior + n_one, b_prior + n_zero)
}
