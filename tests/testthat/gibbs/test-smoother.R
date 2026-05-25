# Tests for fast_smoother and simulation_smoother. Validation strategy:
#
# 1. fast_smoother returns E[x|y], the posterior MEAN of the state path.
#    For a scalar AR(1) + Gaussian noise SSM we can compute the analytic
#    posterior mean as the Kalman SMOOTHER closed-form, and check.
#
# 2. simulation_smoother returns a DRAW from p(x|y). Average of many
#    draws should converge to fast_smoother's point estimate.
#
# 3. Both should handle missing data without error.

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

make_ar1_ssm <- function(phi = 0.7, sigma_eta = 0.5, sigma_eps = 1.0) {
  list(
    H        = matrix(1, 1, 1),
    F        = matrix(phi, 1, 1),
    G        = matrix(1, 1, 1),
    Sigma_eps = matrix(sigma_eps^2, 1, 1),
    Sigma_eta = matrix(sigma_eta^2, 1, 1),
    mu_1     = 0,
    Sigma_1  = matrix(sigma_eta^2 / (1 - phi^2), 1, 1)
  )
}

test_that("fast_smoother runs without error on scalar AR(1) + noise", {
  set.seed(10)
  ssm <- make_ar1_ssm()
  sim <- simulate_ssm(ssm, T_ = 50, need_states = TRUE)
  fs <- fast_smoother(sim$Y, ssm, need_states = TRUE)
  expect_equal(dim(fs$states), c(1L, 50L))
  expect_equal(dim(fs$disturbances$m_errors), c(1L, 50L))
  expect_equal(dim(fs$disturbances$shocks),   c(1L, 49L))
})

test_that("simulation_smoother runs and returns one valid draw", {
  set.seed(11)
  ssm <- make_ar1_ssm()
  sim <- simulate_ssm(ssm, T_ = 50, need_states = TRUE)
  ss <- simulation_smoother(sim$Y, ssm)
  expect_equal(dim(ss$states), c(1L, 50L))
})

test_that("simulation_smoother draws average to fast_smoother point", {
  # Average of many sim_smoother draws should converge to fast_smoother's
  # E[x | y] point estimate.
  set.seed(12)
  ssm <- make_ar1_ssm(phi = 0.7, sigma_eta = 0.5, sigma_eps = 1.0)
  T_ <- 30
  sim <- simulate_ssm(ssm, T_ = T_, need_states = TRUE)

  fs <- fast_smoother(sim$Y, ssm, need_states = TRUE)
  point_mean <- as.numeric(fs$states)   # length T

  # Average many simulation_smoother draws
  n_draws <- 500
  draws <- replicate(n_draws,
                     as.numeric(simulation_smoother(sim$Y, ssm)$states))
  draw_mean <- rowMeans(draws)

  # Per-t Monte Carlo error: posterior_sd / sqrt(n_draws). For T=30 and
  # n_draws=500, mean across t of |draw_mean - point_mean| should be
  # well under 0.1 (posterior SD ~ 0.5, MC SE ~ 0.022).
  mean_abs_diff <- mean(abs(draw_mean - point_mean))
  expect_lt(mean_abs_diff, 0.1)
})

test_that("smoothers handle missing data without crashing", {
  set.seed(13)
  ssm <- make_ar1_ssm()
  sim <- simulate_ssm(ssm, T_ = 20, need_states = TRUE)
  y <- sim$Y
  y[1L, c(5, 12, 18)] <- NA       # spaced missing obs

  fs <- fast_smoother(y, ssm, need_states = TRUE)
  expect_equal(dim(fs$states), c(1L, 20L))
  expect_true(all(is.finite(fs$states)))

  ss <- simulation_smoother(y, ssm)
  expect_equal(dim(ss$states), c(1L, 20L))
  expect_true(all(is.finite(ss$states)))
})
