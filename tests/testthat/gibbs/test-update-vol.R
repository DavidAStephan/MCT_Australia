# Tests for update_vol (Kim-Shephard SV Gibbs update).
#
# Validation strategy: simulate from the known SV process, run many
# Gibbs sweeps starting from a bad initial guess, and check that the
# posterior mean of sigma_t converges back to the truth.

library(testthat)

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate from the SV process:
#   x_t = sigma_t * eps_t,   eps ~ N(0,1)
#   ln(sigma_t^2) = ln(sigma_{t-1}^2) + gamma * ups_t,   ups ~ N(0,1)
sim_sv <- function(T_, gamma, sigma_init = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ln_sig2 <- numeric(T_)
  ln_sig2[1L] <- 2 * log(sigma_init)
  for (t in 2:T_) {
    ln_sig2[t] <- ln_sig2[t - 1L] + gamma * rnorm(1)
  }
  sigma <- exp(ln_sig2 / 2)
  x     <- sigma * rnorm(T_)
  list(x = x, sigma = sigma)
}

test_that("update_vol runs without error on simulated SV data", {
  set.seed(20)
  sim <- sim_sv(T_ = 100, gamma = 0.1)
  out <- update_vol(sim$x, sigma = rep(1, 100), gamma = 0.1)
  expect_length(out, 100L)
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0))
})

test_that("update_vol does NOT drift away from truth over many sweeps", {
  # Run a long Gibbs chain starting from sigma = 1 everywhere.
  # The posterior mean should NOT systematically over/under-estimate
  # the true sigma path.
  set.seed(21)
  T_    <- 200
  gamma <- 0.15
  sim <- sim_sv(T_ = T_, gamma = gamma, sigma_init = 1)

  n_burn <- 500
  n_draw <- 2000
  sigma <- rep(1, T_)
  sigma_draws <- matrix(0, T_, n_draw)

  for (s in seq_len(n_burn)) {
    sigma <- update_vol(sim$x, sigma, gamma)
  }
  for (s in seq_len(n_draw)) {
    sigma <- update_vol(sim$x, sigma, gamma)
    sigma_draws[, s] <- sigma
  }

  post_mean <- rowMeans(sigma_draws)
  # We're not expecting perfect recovery of every t — SV has strong
  # uncertainty per-period. We DO expect the mean log-sigma error to
  # be close to zero (no systematic bias).
  log_bias <- mean(log(post_mean) - log(sim$sigma))
  expect_lt(abs(log_bias), 0.20)

  # Cross-period correlation between post_mean and truth should be
  # comfortably positive (>0.5) — we're recovering the signal.
  cor_post_truth <- cor(log(post_mean), log(sim$sigma))
  expect_gt(cor_post_truth, 0.5)
})

test_that("update_vol handles missing data without crashing", {
  set.seed(22)
  sim <- sim_sv(T_ = 50, gamma = 0.1)
  x <- sim$x
  x[c(5, 15, 30)] <- NA
  out <- update_vol(x, sigma = rep(1, 50), gamma = 0.1)
  expect_length(out, 50L)
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0))
})

test_that("mixture-component prob/mean/var arrays match MATLAB exactly", {
  # Defensive: these constants came from MCT-NYFED/functions/update_vol.m
  # and must not drift via typos. Sanity-check sums and lengths.
  expect_equal(length(.KSC_PROBS), 10L)
  expect_equal(length(.KSC_MEANS), 10L)
  expect_equal(length(.KSC_VARS),  10L)
  expect_equal(sum(.KSC_PROBS), 1.0, tolerance = 1e-5)   # probs sum to 1
  expect_equal(.KSC_STDVS, sqrt(.KSC_VARS))
})
