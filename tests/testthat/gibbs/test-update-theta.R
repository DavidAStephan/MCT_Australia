# Tests for update_theta (MA(q) coefficient conjugate update).

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate from MA(q): u_t = eps_t + sum_{l=1..q} theta_l * eps_{t-l}
sim_ma <- function(T_, theta, sigma = 1, seed = 1L) {
  set.seed(seed)
  q <- length(theta)
  eps <- rnorm(T_ + q, sd = sigma)  # extra at start for lags
  u <- numeric(T_)
  for (t in seq_len(T_)) {
    # eps at position (q + t) is the current eps_t; lags are q+t-1..q+t-q
    lag_terms <- sum(theta * eps[(q + t - 1):(q + t - q)])
    u[t] <- eps[q + t] + lag_terms
  }
  list(u = u, eps = eps, theta = theta)
}

test_that("update_theta runs without error and returns length-q vector", {
  set.seed(100)
  T_ <- 200
  theta_true <- c(0.5, -0.3, 0.1)
  sim <- sim_ma(T_, theta_true, seed = 100)

  # Build regression inputs: y_t = u_t - eps_t = sum theta_l * eps_{t-l}
  # We'd update theta from the conjugate regression on lagged eps's.
  # For this test we feed it the right structure.
  y <- numeric(T_)
  x <- matrix(0, T_, length(theta_true))
  for (t in seq_len(T_)) {
    y[t] <- sim$u[t] - sim$eps[length(theta_true) + t]
    x[t, ] <- sim$eps[length(theta_true) + t - seq_len(length(theta_true))]
  }
  theta_draw <- update_theta(y, x, prec_prior = 0.1)
  expect_length(theta_draw, 3L)
  expect_true(all(is.finite(theta_draw)))
})

test_that("update_theta posterior concentrates around truth as T grows", {
  # With many obs from a known MA(2), repeated draws should have
  # posterior mean close to truth.
  set.seed(101)
  T_ <- 5000
  theta_true <- c(0.4, -0.2)
  sim <- sim_ma(T_, theta_true, seed = 101)

  y <- numeric(T_)
  x <- matrix(0, T_, length(theta_true))
  for (t in seq_len(T_)) {
    y[t] <- sim$u[t] - sim$eps[length(theta_true) + t]
    x[t, ] <- sim$eps[length(theta_true) + t - seq_len(length(theta_true))]
  }
  draws <- replicate(200, update_theta(y, x, prec_prior = 0.1))
  # draws is q x 200; row k is theta_k draws
  for (k in seq_along(theta_true)) {
    med_k <- median(draws[k, ])
    expect_lt(abs(med_k - theta_true[k]), 0.05,
              label = sprintf("theta[%d] median %.3f vs truth %.3f",
                              k, med_k, theta_true[k]))
  }
})

test_that("update_theta enforces invertibility (all MA roots inside unit circle)", {
  # If the unconstrained posterior would put theta in a non-invertible
  # region, the invertibility step should reflect roots so the result
  # IS invertible. Pick a tight prior toward a non-invertible value.
  set.seed(102)
  T_ <- 50
  # Construct y/x so that OLS theta would be ~ (1.5, 0.5) — outside
  # the unit circle. The MA polynomial 1 + 1.5z + 0.5z^2 has roots
  # at z = -1 and z = -2 (the -2 is outside).
  q <- 2
  x <- matrix(rnorm(T_ * q), T_, q)
  theta_target <- c(1.5, 0.5)
  y <- as.numeric(x %*% theta_target) + 0.01 * rnorm(T_)  # very low noise
  theta_draw <- update_theta(y, x, prec_prior = 0.001)   # weak prior

  # The drawn theta's MA polynomial should have all roots inside
  # the unit circle (i.e. |z| > 1 for the polynomial in inverse-lag form
  # — equivalent to |z| < 1 for the polynomial in lag form).
  # Using MATLAB convention: poly is [theta_q, ..., theta_1, 1]
  poly_coefs <- c(rev(theta_draw), 1)
  roots <- polyroot(poly_coefs)
  # All roots should have |z| >= 1 after invertibility enforcement
  # (because they were reflected to reciprocals if originally < 1)
  expect_true(all(Mod(roots) >= 1 - 1e-6),
              label = sprintf("MA roots after reflection: %s",
                              paste(sprintf("%.3f", Mod(roots)),
                                    collapse = ",")))
})
