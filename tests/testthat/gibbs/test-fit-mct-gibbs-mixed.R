# Step 9a tests: mixed-frequency observation handling.
#
# Strategy:
#   1. fit_mct_gibbs with obs_type = all-1 (monthly everywhere) should
#      run without error and produce sensible output (the augmented-state
#      mixed-freq path subsumes the monthly-only Step 7 case).
#   2. fit_mct_gibbs on truly mixed-freq simulated data (monthly for
#      latest N months, quarterly average pre-cutover) should also run
#      and recover signal (loose thresholds — slow Gibbs mixing).
#   3. The augmented state correctly handles quarterly obs equation
#      (positions [1..3] for c lags, [3i+1..3i+3] for s_i lags).

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Inline simulator from earlier files
sim_mct_a <- function(T_, N, rho = 0.7, lambda = NULL,
                      sigma_c_true = 0.4, sigma_s_true = 0.25,
                      sigma_eps_true = 0.8, ref = 1L, seed = 1L) {
  set.seed(seed)
  if (is.null(lambda)) lambda <- rep(1, N)
  lambda[ref] <- 1
  sigma_c   <- rep(sigma_c_true, T_)
  sigma_s   <- matrix(sigma_s_true, T_, N)
  sigma_eps <- matrix(sigma_eps_true, T_, N)
  c_path <- numeric(T_)
  c_path[1L] <- sigma_c[1L] / sqrt(1 - rho^2) * rnorm(1)
  for (t in 2:T_) c_path[t] <- rho * c_path[t - 1L] + sigma_c[t] * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) {
    s_path[t, ] <- s_path[t - 1L, ] + sigma_s[t, ] * rnorm(N)
  }
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- lambda * c_path[t] + s_path[t, ] + sigma_eps[t, ] * rnorm(N)
  }
  list(y = y, c = c_path, s = s_path, rho = rho)
}

# Convert a monthly y matrix to mixed-frequency: monthly for the last
# `n_monthly_tail` periods, quarterly average (at last month of quarter)
# for everything before that. Returns (y_mixed, obs_type) pair.
to_mixed_freq <- function(y, n_monthly_tail) {
  T_ <- nrow(y); N <- ncol(y)
  y_mixed <- matrix(NA_real_, T_, N)
  obs_type <- matrix(0L, T_, N)
  cutover <- T_ - n_monthly_tail
  # Pre-cutover: quarterly average at last month of each quarter
  for (t in seq_len(cutover)) {
    if (t >= 3 && t %% 3 == 0) {
      y_mixed[t, ] <- colMeans(y[(t - 2L):t, , drop = FALSE])
      obs_type[t, ] <- 2L
    }
  }
  # Post-cutover: monthly
  for (t in (cutover + 1L):T_) {
    y_mixed[t, ] <- y[t, ]
    obs_type[t, ] <- 1L
  }
  list(y = y_mixed, obs_type = obs_type)
}

test_that("fit_mct_gibbs with obs_type=all-1 runs end-to-end", {
  set.seed(70)
  sim <- sim_mct_a(T_ = 50, N = 3, seed = 70)
  obs_type <- matrix(1L, 50, 3)
  fit <- fit_mct_gibbs(sim$y, obs_type = obs_type, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  expect_equal(dim(fit$draws$c), c(50L, 100L))
  expect_equal(dim(fit$draws$s), c(50L, 3L, 100L))
  expect_length(fit$draws$rho, 100L)
})

test_that("mixed-freq augmented-state c extraction matches monthly path", {
  # On identical monthly data, the mixed-freq path with all obs_type=1
  # uses an augmented 3(N+1)-dim state but should extract c and s with
  # values comparable in scale to the monthly-only Step 7 path. We
  # don't expect EXACT match (different RNG sequence in the augmented
  # SSM), but moments should be in the same range.
  set.seed(71)
  sim <- sim_mct_a(T_ = 80, N = 3, seed = 71)

  fit_mo <- fit_mct_gibbs(sim$y, ref = 1L,
                          n_burn = 100L, n_draw = 200L, verbose = FALSE)
  obs_type <- matrix(1L, 80, 3)
  fit_mx <- fit_mct_gibbs(sim$y, obs_type = obs_type, ref = 1L,
                          n_burn = 100L, n_draw = 200L, verbose = FALSE)

  # Compare posterior medians of c at each t
  c_mo <- apply(fit_mo$draws$c, 1, median)
  c_mx <- apply(fit_mx$draws$c, 1, median)
  # Should correlate strongly (both estimating the same posterior)
  expect_gt(cor(c_mo, c_mx), 0.85)
})

test_that("fit_mct_gibbs handles mixed-freq simulated data without error", {
  set.seed(72)
  T_ <- 120
  sim <- sim_mct_a(T_ = T_, N = 4, seed = 72)
  mixed <- to_mixed_freq(sim$y, n_monthly_tail = 12)
  expect_equal(sum(mixed$obs_type == 0L) + sum(mixed$obs_type == 1L) +
                 sum(mixed$obs_type == 2L), T_ * 4L)

  fit <- fit_mct_gibbs(mixed$y, obs_type = mixed$obs_type, ref = 1L,
                       n_burn = 100L, n_draw = 200L, verbose = FALSE)
  expect_equal(dim(fit$draws$c), c(T_, 200L))
  # Finite outputs
  expect_true(all(is.finite(fit$draws$c)))
  expect_true(all(is.finite(fit$draws$s)))
})

test_that("mixed-freq Gibbs recovers some c signal on simulated data", {
  set.seed(73)
  T_ <- 120
  sim <- sim_mct_a(T_ = T_, N = 5, seed = 73)
  mixed <- to_mixed_freq(sim$y, n_monthly_tail = 24)
  fit <- fit_mct_gibbs(mixed$y, obs_type = mixed$obs_type, ref = 1L,
                       n_burn = 300L, n_draw = 600L, verbose = FALSE)
  c_post_mean <- rowMeans(fit$draws$c)
  # Loose threshold — mixed-freq + short chain + only 24/120 monthly
  # obs gives weaker identification than the all-monthly case.
  cor_c <- cor(c_post_mean, sim$c)
  expect_gt(cor_c, 0.3)
})
