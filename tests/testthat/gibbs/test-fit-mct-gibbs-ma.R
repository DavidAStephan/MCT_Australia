# Test the MA(q) path of fit_mct_gibbs end-to-end.
# Generate data with KNOWN MA(q) noise structure and check that the
# sampler runs and produces sensible theta posteriors.

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate from MA(q)-augmented MCT model. Uses the NY Fed convention
# where eps state has variance sigma_eps^2.
sim_mct_a_ma <- function(T_, N, rho = 0.7, theta = NULL,
                         sigma_c_true = 0.4, sigma_s_true = 0.25,
                         sigma_eps_true = 0.8, seed = 1L) {
  if (is.null(theta)) theta <- matrix(0, N, 0)
  q <- ncol(theta)
  set.seed(seed)
  c_path <- numeric(T_)
  c_path[1L] <- sigma_c_true / sqrt(1 - rho^2) * rnorm(1)
  for (t in 2:T_) c_path[t] <- rho * c_path[t - 1L] + sigma_c_true * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) s_path[t, ] <- s_path[t - 1L, ] + sigma_s_true * rnorm(N)
  # Generate scaled eps directly: eps[t, i] ~ N(0, sigma_eps_true^2)
  eps <- matrix(sigma_eps_true * rnorm(T_ * N), T_, N)
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    ma_noise <- eps[t, ]
    if (q > 0) {
      for (l in seq_len(q)) {
        if (t - l >= 1) ma_noise <- ma_noise + theta[, l] * eps[t - l, ]
      }
    }
    y[t, ] <- c_path[t] + s_path[t, ] + ma_noise  # lambda = 1 throughout
  }
  list(y = y, c = c_path, s = s_path, eps = eps, theta = theta)
}

test_that("fit_mct_gibbs ma path runs end-to-end with q_MA=2", {
  set.seed(140)
  T_ <- 60; N <- 3
  theta_true <- matrix(c(0.4, -0.2,
                          0.3, -0.1,
                          0.5, 0.0), N, 2, byrow = TRUE)
  sim <- sim_mct_a_ma(T_, N, theta = theta_true, seed = 140)
  fit <- fit_mct_gibbs(sim$y, q_MA = 2L, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  expect_equal(dim(fit$draws$theta), c(N, 2L, 100L))
  expect_equal(dim(fit$draws$c), c(T_, 100L))
  expect_true(all(is.finite(fit$draws$rho)))
})

test_that("MA(2) path recovers theta posterior median direction on sim", {
  # Loose recovery test — short chain, small data, just check
  # we don't recover NEGATIVE theta when truth is POSITIVE.
  set.seed(141)
  T_ <- 150; N <- 3
  theta_true <- matrix(0.4, N, 1)   # MA(1) with positive theta
  sim <- sim_mct_a_ma(T_, N, theta = theta_true, seed = 141)
  fit <- fit_mct_gibbs(sim$y, q_MA = 1L, ref = 1L,
                       n_burn = 200L, n_draw = 400L, verbose = FALSE)
  # Posterior median of theta_1 per sector
  theta_med <- apply(fit$draws$theta, c(1, 2), median)
  # At least 2 of 3 sectors should have positive theta median
  # (matches positive truth; allow 1 sector to fail under MC noise)
  expect_gte(sum(theta_med[, 1] > 0), 2,
             label = sprintf("theta_1 medians: %s",
                             paste(sprintf("%.2f", theta_med[, 1]),
                                   collapse = ", ")))
})

test_that("fit_mct_gibbs rejects combined mixed-freq + MA path", {
  T_ <- 20; N <- 3
  y <- matrix(rnorm(T_ * N), T_, N)
  obs_type <- matrix(1L, T_, N)
  expect_error(
    fit_mct_gibbs(y, obs_type = obs_type, q_MA = 2L, ref = 1L,
                  n_burn = 10L, n_draw = 10L, verbose = FALSE),
    "mixed-freq \\+ MA"
  )
})
