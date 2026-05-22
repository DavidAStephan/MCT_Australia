# Tests for postprocess helpers. The heavy lifting (extracting from a fit)
# is exercised end-to-end after the real fits complete; here we cover the
# pure-numerics helper `summarise_path` and the array-reshape convention.

library(testthat)
source(file.path(testthat::test_path(), "..", "..", "R", "postprocess.R"))

test_that("summarise_path returns one row per date with monotone bands", {
  set.seed(1)
  D <- 400; T_ <- 20
  truth <- seq(-2, 5, length.out = T_)
  draws <- matrix(rnorm(D * T_, mean = rep(truth, each = D), sd = 1), D, T_)
  dates <- seq(as.Date("2020-01-01"), by = "month", length.out = T_)

  out <- summarise_path(draws, dates)
  expect_equal(nrow(out), T_)
  expect_equal(out$date, dates)
  expect_true(all(out$lower < out$median))
  expect_true(all(out$median < out$upper))
  # Medians should track the truth to within a few SD/sqrt(D) ~ 0.05
  expect_lt(max(abs(out$median - truth)), 0.25)
})

test_that("summarise_path validates inputs", {
  draws <- matrix(rnorm(10 * 5), 10, 5)
  expect_error(summarise_path(draws, dates = seq(as.Date("2020-01-01"),
                                                 by = "month", length.out = 6)))
  expect_error(summarise_path(draws,
                              dates = seq(as.Date("2020-01-01"),
                                          by = "month", length.out = 5),
                              probs = c(0.5, 0.5, 0.5)))  # not increasing
})

test_that("Stan column-major flattening reshapes correctly for (T, N) matrices", {
  # Build a known (D, T, N) array and round-trip it through the same flat
  # ordering the draws_matrix uses: column k is sector_trend[t, i] with
  # t = ((k-1) mod T) + 1, i = ((k-1) %/% T) + 1.
  D <- 3; T_ <- 4; N <- 2
  truth <- array(seq_len(D * T_ * N), dim = c(D, T_, N))
  # Flatten column-major (matching cmdstanr's draws_matrix column order)
  flat <- matrix(as.numeric(truth), nrow = D)
  expect_equal(ncol(flat), T_ * N)
  # Reshape back exactly as extract_matrix_array does
  back <- array(as.numeric(flat), dim = c(D, T_, N))
  expect_equal(back, truth)
})
