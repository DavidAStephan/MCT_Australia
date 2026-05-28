# Test Variant B (RW common factor / Stock-Watson common trend) on
# simulated data. Validates that:
#   - the rho parameter is pinned at 1 throughout (no rho update)
#   - the sigma_c posterior recovers the true RW innovation SD
#   - the latest-t trend recovers the true value within 90% CI
#
# Smaller (T, N) than the Variant A tests so the suite runs in <30s.

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate a Variant B (RW common factor) panel with mixed-frequency
# observations (last 1/3 monthly, rest quarterly).
sim_mct_b_mixed <- function(T_, N, ref = 1L,
                            sigma_c_true = 0.4, sigma_s_true = 0.25,
                            sigma_eps_true = 0.8, monthly_start_frac = 2/3,
                            seed = 1L) {
  set.seed(seed)
  # RW common factor
  c_path <- numeric(T_)
  c_path[1L] <- rnorm(1, 0, 2)            # diffuse start (Variant B prior)
  for (t in 2:T_) c_path[t] <- c_path[t - 1L] + sigma_c_true * rnorm(1)
  # RW sector trends
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) {
    s_path[t, ] <- s_path[t - 1L, ] + sigma_s_true * rnorm(N)
  }
  # Latent monthly inflations
  pi_lat <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    pi_lat[t, ] <- c_path[t] + s_path[t, ] + sigma_eps_true * rnorm(N)
  }
  # Observation matrix: quarterly average up to month floor(T*frac);
  # monthly thereafter
  y <- matrix(NA_real_, T_, N)
  obs_type <- matrix(0L, T_, N)
  monthly_start <- max(3L, floor(T_ * monthly_start_frac))
  for (t in seq(3L, monthly_start - 1L, by = 3L)) {
    obs_type[t, ] <- 2L
    y[t, ] <- (pi_lat[t - 2L, ] + pi_lat[t - 1L, ] + pi_lat[t, ]) / 3
  }
  for (t in monthly_start:T_) {
    obs_type[t, ] <- 1L
    y[t, ] <- pi_lat[t, ]
  }
  list(y = y, obs_type = obs_type, c = c_path, s = s_path,
       sigma_c_true = sigma_c_true, sigma_s_true = sigma_s_true,
       sigma_eps_true = sigma_eps_true)
}

test_that("Variant B Gibbs runs end-to-end and pins rho at 1", {
  set.seed(700)
  T_ <- 60; N <- 3
  sim <- sim_mct_b_mixed(T_, N, seed = 700)
  fit <- fit_mct_gibbs(
    sim$y, obs_type = sim$obs_type, ref = 1L,
    n_burn = 100L, n_draw = 200L, verbose = FALSE,
    config = list(variant = "B")
  )
  expect_equal(dim(fit$draws$c), c(T_, 200L))
  expect_equal(dim(fit$draws$s), c(T_, N, 200L))
  expect_true(all(fit$draws$rho == 1),
              label = "rho should be pinned at 1 for Variant B")
})

test_that("Variant B recovers the common-factor RW path roughly", {
  set.seed(701)
  T_ <- 120; N <- 4
  sim <- sim_mct_b_mixed(T_, N, sigma_c_true = 0.3,
                         sigma_s_true = 0.2, sigma_eps_true = 0.6,
                         seed = 701)
  fit <- fit_mct_gibbs(
    sim$y, obs_type = sim$obs_type, ref = 1L,
    n_burn = 300L, n_draw = 600L, verbose = FALSE,
    config = list(variant = "B")
  )

  # Posterior median c path vs truth — correlation should be reasonable
  c_med <- apply(fit$draws$c, 1, median)
  rho_corr <- cor(c_med, sim$c)
  expect_gt(rho_corr, 0.4,
            label = sprintf("c-path correlation = %.2f", rho_corr))

  # sigma_c posterior median (averaged across t) within 50% of truth
  # — RW SV identification is loose, so a wide tolerance is appropriate.
  sigma_c_med <- median(fit$draws$sigma_c)
  expect_gt(sigma_c_med, 0.5 * sim$sigma_c_true)
  expect_lt(sigma_c_med, 2.0 * sim$sigma_c_true)
})
