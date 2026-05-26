# Tests for the marginal MH update for rho (R/gibbs/update_rho_mh.R).
#
# Implementation notes (corrected understanding 2026-05-27):
#   - The marginal MH for rho was implemented under the suspicion that
#     the conditional Gibbs update had a downward bias on rho.
#   - The on-real-data validation in Step 11 showed both updates give
#     ESSENTIALLY THE SAME posterior — and that posterior matches Stan
#     (e.g. on real ABS data both give rho ≈ 0.36, matching Stan's 0.37).
#   - The "rho bias on sim data" diagnosis was actually a finite-sample
#     identification issue: with sigma_eps=0.8 and sigma_c=0.4, the
#     marginal log-lik of rho really IS centred near 0.45, not at the
#     truth 0.7. Both sampling strategies correctly recover this.
#   - The MH is therefore not faster or more accurate. It's available as
#     `use_marginal_mh_rho = TRUE` (default for safety; set FALSE to
#     skip the extra KF pass per iter and save ~50% wall time).
#
# These tests verify that MH and conditional Gibbs give CONSISTENT
# posteriors (within MC noise) — confirming both are sampling the true
# posterior, not that either recovers ground truth (which the data
# doesn't support at this SNR).

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

sim_mct_a <- function(T_, N, rho = 0.7,
                      sigma_c_true = 0.4, sigma_s_true = 0.25,
                      sigma_eps_true = 0.8, seed = 1L) {
  set.seed(seed)
  sigma_c <- rep(sigma_c_true, T_)
  sigma_s <- matrix(sigma_s_true, T_, N)
  sigma_eps <- matrix(sigma_eps_true, T_, N)
  c_path <- numeric(T_)
  c_path[1L] <- sigma_c[1L] / sqrt(1 - rho^2) * rnorm(1)
  for (t in 2:T_) c_path[t] <- rho * c_path[t - 1L] + sigma_c[t] * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) s_path[t, ] <- s_path[t - 1L, ] + sigma_s[t, ] * rnorm(N)
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- c_path[t] + s_path[t, ] + sigma_eps[t, ] * rnorm(N)
  }
  list(y = y, c = c_path, s = s_path, rho = rho)
}

test_that("marginal MH runs and emits acceptance rate", {
  set.seed(80)
  sim <- sim_mct_a(T_ = 80, N = 4, seed = 80)
  obs_type <- matrix(1L, 80, 4)
  fit <- fit_mct_gibbs(sim$y, obs_type = obs_type, ref = 1L,
                       n_burn = 100L, n_draw = 200L, verbose = FALSE,
                       config = list(use_marginal_mh_rho = TRUE,
                                     rho_prop_sd = 0.05))
  expect_true(!is.null(fit$rho_accept_rate))
  expect_gt(fit$rho_accept_rate, 0.30)   # healthy acceptance
  expect_lt(fit$rho_accept_rate, 1.00)
  expect_length(fit$draws$rho, 200L)
})

test_that("MH and conditional Gibbs give consistent rho posteriors", {
  # On the same data, the two update strategies should land at the
  # SAME posterior (both are valid samplers of the same target).
  # Tight test: posterior medians within 0.1 of each other.
  set.seed(81)
  sim <- sim_mct_a(T_ = 150, N = 4, rho = 0.7, seed = 81)
  obs_type <- matrix(1L, 150, 4)

  fit_mh <- fit_mct_gibbs(
    sim$y, obs_type = obs_type, ref = 1L,
    n_burn = 500L, n_draw = 1000L, verbose = FALSE,
    config = list(use_marginal_mh_rho = TRUE, rho_prop_sd = 0.05)
  )
  fit_gibbs <- fit_mct_gibbs(
    sim$y, obs_type = obs_type, ref = 1L,
    n_burn = 500L, n_draw = 1000L, verbose = FALSE,
    config = list(use_marginal_mh_rho = FALSE)
  )

  rho_mh    <- median(fit_mh$draws$rho)
  rho_gibbs <- median(fit_gibbs$draws$rho)
  # MC tolerance on single chains of this length
  expect_lt(abs(rho_mh - rho_gibbs), 0.1,
            label = sprintf("MH median %.3f vs Gibbs median %.3f",
                            rho_mh, rho_gibbs))
})

test_that("acceptance rate scales sensibly with prop_sd", {
  # Smaller proposal sd ⇒ higher acceptance (proposals stay close).
  set.seed(82)
  sim <- sim_mct_a(T_ = 80, N = 4, seed = 82)
  obs_type <- matrix(1L, 80, 4)

  acc_small <- fit_mct_gibbs(
    sim$y, obs_type = obs_type, ref = 1L,
    n_burn = 100L, n_draw = 300L, verbose = FALSE,
    config = list(use_marginal_mh_rho = TRUE, rho_prop_sd = 0.03)
  )$rho_accept_rate
  acc_large <- fit_mct_gibbs(
    sim$y, obs_type = obs_type, ref = 1L,
    n_burn = 100L, n_draw = 300L, verbose = FALSE,
    config = list(use_marginal_mh_rho = TRUE, rho_prop_sd = 0.30)
  )$rho_accept_rate
  expect_gt(acc_small, acc_large)
})
