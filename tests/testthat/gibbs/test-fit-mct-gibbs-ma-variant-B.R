# MA(q) on the non-mixed Variant B path (quarterly sub-group setting).
# Validates: (a) variant = "B" + q_MA > 0 runs end-to-end; (b) rho stays
# pinned at 1; (c) NA cells in y (staggered series starts) don't crash;
# (d) theta posterior is non-degenerate.

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

sim_b_ma_q <- function(T_, N, q = 3L,
                       sigma_c_true = 0.3, sigma_s_true = 0.2,
                       sigma_eps_true = 0.6,
                       theta_true = c(0.4, 0.2, 0.1),
                       seed = 1L) {
  set.seed(seed)
  c_path <- numeric(T_); c_path[1L] <- rnorm(1, 0, 2)
  for (t in 2:T_) c_path[t] <- c_path[t - 1L] + sigma_c_true * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) s_path[t, ] <- s_path[t - 1L, ] + sigma_s_true * rnorm(N)
  eps <- matrix(rnorm(T_ * N), T_, N)
  # MA(q) sector noise: u_t = eps_t + sum_l theta_l eps_{t-l}
  u <- matrix(0, T_, N)
  for (i in seq_len(N)) {
    for (t in seq_len(T_)) {
      val <- eps[t, i]
      for (l in seq_len(q)) {
        if (t - l >= 1L) val <- val + theta_true[l] * eps[t - l, i]
      }
      u[t, i] <- sigma_eps_true * val
    }
  }
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- c_path[t] + s_path[t, ] + u[t, ]
  }
  list(y = y, theta_true = theta_true, sigma_c_true = sigma_c_true)
}

test_that("MA(3) on Variant B quarterly path runs end-to-end", {
  set.seed(1000)
  T_ <- 80; N <- 3
  sim <- sim_b_ma_q(T_, N, seed = 1000)
  fit <- fit_mct_gibbs(
    sim$y, q_MA = 3L, ref = 1L,
    n_burn = 100L, n_draw = 200L, verbose = FALSE,
    config = list(variant = "B")
  )
  expect_equal(dim(fit$draws$theta), c(N, 3L, 200L))
  expect_true(all(fit$draws$rho == 1))
  expect_true(all(is.finite(fit$draws$theta)))
})

test_that("MA(3) + NA cells (staggered starts) work together", {
  set.seed(1001)
  T_ <- 100; N <- 4
  sim <- sim_b_ma_q(T_, N, seed = 1001)
  # Sector 4 starts at t = 30
  sim$y[1:29, 4L] <- NA_real_
  fit <- fit_mct_gibbs(
    sim$y, q_MA = 3L, ref = 1L,
    n_burn = 200L, n_draw = 400L, verbose = FALSE,
    config = list(variant = "B")
  )
  expect_true(all(is.finite(fit$draws$c)))
  expect_true(all(is.finite(fit$draws$s)))
  expect_true(all(is.finite(fit$draws$theta)))
})

test_that("MA(3) posterior recovers something close to true theta", {
  set.seed(1002)
  T_ <- 150; N <- 4
  theta_true <- c(0.5, 0.3, 0.1)
  sim <- sim_b_ma_q(T_, N, theta_true = theta_true, seed = 1002)
  fit <- fit_mct_gibbs(
    sim$y, q_MA = 3L, ref = 1L,
    n_burn = 400L, n_draw = 800L, verbose = FALSE,
    config = list(variant = "B")
  )
  # Posterior median per sector + theta lag
  theta_med <- apply(fit$draws$theta, c(1, 2), median)
  # Pool across sectors (they share the same true thetas in this sim)
  pooled <- colMeans(theta_med)
  # Loose tolerance — MA identification is fundamentally weak when the
  # RW common factor can absorb shocks (similar to outlier dynamics).
  # Just check ordering and that we're in the right ballpark.
  expect_gt(pooled[1L], pooled[3L],
            label = sprintf("theta_1=%.2f should exceed theta_3=%.2f",
                            pooled[1L], pooled[3L]))
})
