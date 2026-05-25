# R port of MCT-NYFED/functions/update_gam.m
#
# Inverse-Gamma conjugate update for the scale parameter in:
#
#   x_t = gamma * eps_t,    eps_t ~ N(0, 1).
#
# Prior: gamma^2 ~ Inverse-Gamma(nu_prior/2, nu_prior * s2_prior / 2).
# Posterior given T observations:
#   nu_post = nu_prior + T
#   s2_post = (nu_prior / nu_post) * s2_prior +
#             (1 / nu_post) * sum(x^2)
#   gamma   = 1 / sqrt(Gamma(nu_post/2, scale = 2 / (nu_post * s2_post)))
#
# Used by the main Gibbs loop to update the RW step sizes for each of
# the SV log-volatility processes (sigma_dtau_c, sigma_dtau_i, etc.).
#
# Vectorised over N independent gammas: if x is T x N, returns N draws.

#' Draw new scale(s) gamma from Inverse-Gamma posterior.
#'
#' @param x        T x N matrix (or T vector) of residuals.
#' @param nu_prior Length-N vector (or scalar) of prior degrees of freedom.
#' @param s2_prior Length-N vector (or scalar) of prior scale.
#' @return Length-N vector of new gamma draws.
update_gam <- function(x, nu_prior, s2_prior) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  N        <- ncol(x)
  T_       <- nrow(x)
  nu_prior <- rep_len(as.numeric(nu_prior), N)
  s2_prior <- rep_len(as.numeric(s2_prior), N)

  nu_post <- nu_prior + T_
  ss      <- colSums(x^2)
  s2_post <- (nu_prior / nu_post) * s2_prior + (1 / nu_post) * ss

  # gamma^2 ~ IG(nu_post/2, nu_post*s2_post/2)  ⇔
  # 1/gamma^2 ~ Gamma(shape = nu_post/2, scale = 2/(nu_post*s2_post))
  gam_sq_inv <- rgamma(N, shape = nu_post / 2,
                       scale = 2 / (nu_post * s2_post))
  1 / sqrt(gam_sq_inv)
}
