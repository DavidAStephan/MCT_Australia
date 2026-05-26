# Tests for update_scl + update_ps (outlier scale mixture).

library(testthat)

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

test_that("update_scl returns values from the candidate grid", {
  set.seed(150)
  vals  <- c(1, 2, 5, 10)
  probs <- c(0.97, 0.01, 0.01, 0.01)
  x <- rnorm(100)
  s <- update_scl(x, vals, probs)
  expect_length(s, 100L)
  expect_true(all(s %in% vals))
})

test_that("update_scl flags large |x| as outliers", {
  # With a clear outlier (x=20), the posterior should heavily favour
  # the largest scale.
  set.seed(151)
  vals  <- c(1, 2, 5, 10)
  probs <- c(0.97, 0.01, 0.01, 0.01)
  draws <- replicate(500, update_scl(20, vals, probs))
  # Median sampled scale should be the largest val (10) — sometimes 5
  # acceptable; but rarely 1 (P(x=20 | s=1) is astronomically small).
  med <- median(draws)
  expect_gte(med, 5)
})

test_that("update_scl returns vals[1] for small |x|", {
  # x close to 0 should always sample the smallest scale.
  set.seed(152)
  vals  <- c(1, 2, 5, 10)
  probs <- c(0.97, 0.01, 0.01, 0.01)
  draws <- replicate(500, update_scl(0.1, vals, probs))
  expect_equal(median(draws), 1)
})

test_that("update_ps recovers truth as T grows", {
  # ps_true = 0.9, T = 10,000. Posterior should concentrate near 0.9.
  set.seed(153)
  ps_true <- 0.9
  T_ <- 10000
  x <- matrix(as.integer(runif(T_) < ps_true), T_, 1)
  draws <- replicate(500, update_ps(x, a_prior = 1, b_prior = 1))
  expect_lt(abs(median(draws) - ps_true), 0.01)
})

test_that("update_ps handles vector N", {
  set.seed(154)
  N <- 3
  ps_true <- c(0.9, 0.5, 0.3)
  T_ <- 5000
  x <- matrix(0, T_, N)
  for (i in seq_len(N)) x[, i] <- as.integer(runif(T_) < ps_true[i])
  draws <- replicate(200, update_ps(x, a_prior = rep(1, N), b_prior = rep(1, N)))
  meds <- rowMeans(draws)
  for (i in seq_len(N)) {
    expect_lt(abs(meds[i] - ps_true[i]), 0.03,
              label = sprintf("ps[%d] %.3f vs truth %.2f",
                              i, meds[i], ps_true[i]))
  }
})
