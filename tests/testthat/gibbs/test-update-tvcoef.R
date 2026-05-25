# Tests for update_tvcoef (Carter-Kohn FFBS for TVP regression).

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate the TVP regression model:
#   y_t          = (alpha_t * x_t) + eps_t,    eps ~ N(0, sigma^2)
#   alpha_t      = alpha_{t-1} + lambda * ups,  ups ~ N(0, 1)
sim_tvp <- function(T_, N, K, lambda, sigma, seed = 1) {
  set.seed(seed)
  alpha <- matrix(0, T_, N * K)  # row t = vec(alpha_t)
  # init at 1 to give the signal something to work with
  alpha[1L, ] <- 1
  for (t in 2:T_) {
    eta <- rep(lambda, K) * rnorm(N * K)
    alpha[t, ] <- alpha[t - 1L, ] + eta
  }
  x <- matrix(rnorm(T_ * K), T_, K)
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    A_t <- matrix(alpha[t, ], N, K)   # vec stacks columns; reshape back
    y[t, ] <- as.numeric(A_t %*% x[t, ]) + sigma * rnorm(N)
  }
  list(y = y, x = x, alpha = alpha)
}

test_that("update_tvcoef runs on small TVP regression", {
  set.seed(40)
  N <- 2; K <- 2; T_ <- 50
  sim <- sim_tvp(T_, N, K, lambda = 0.05, sigma = 0.3)
  sigma_mat <- matrix(0.3, T_, N)
  out <- update_tvcoef(sim$y, sim$x, sigma_mat, lambda = rep(0.05, N))
  expect_equal(dim(out), c(T_, N * K))
  expect_true(all(is.finite(out)))
})

test_that("update_tvcoef draws recover the truth signal direction", {
  # Average many sweeps of FFBS-conditional-on-fixed-hyperparams; the
  # posterior mean should track the true alpha path (modulo MC noise).
  set.seed(41)
  N <- 1; K <- 1; T_ <- 100
  sim <- sim_tvp(T_, N, K, lambda = 0.1, sigma = 0.5)
  sigma_mat <- matrix(0.5, T_, N)

  draws <- replicate(200,
    update_tvcoef(sim$y, sim$x, sigma_mat, lambda = 0.1))
  post_mean <- apply(draws, c(1, 2), mean)

  # Correlation between recovered alpha path and truth — should be high
  # because the obs noise is small and we're using the true hyperparams.
  cor_post_truth <- cor(as.numeric(post_mean), as.numeric(sim$alpha))
  expect_gt(cor_post_truth, 0.8)
})

test_that("update_tvcoef handles missing observations", {
  set.seed(42)
  N <- 2; K <- 2; T_ <- 30
  sim <- sim_tvp(T_, N, K, lambda = 0.05, sigma = 0.3)
  y <- sim$y
  y[c(5, 12, 20), 1] <- NA
  sigma_mat <- matrix(0.3, T_, N)
  out <- update_tvcoef(y, sim$x, sigma_mat, lambda = rep(0.05, N))
  expect_equal(dim(out), c(T_, N * K))
  expect_true(all(is.finite(out)))
})
