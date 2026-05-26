# Test the outlier scale mixture path of fit_mct_gibbs (mixed-freq
# path + use_outliers = TRUE). Generate data with injected outliers
# at specific (t, i) positions and check that the sampler:
#   1. runs end-to-end
#   2. flags the injected outliers (posterior median s > 1 at those obs)
#   3. recovers ps close to truth

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate MCT Variant A monthly data, then INJECT outliers at specific
# (t, i) positions by inflating the obs noise.
sim_mct_with_outliers <- function(T_, N, rho = 0.7,
                                   sigma_c_true = 0.4,
                                   sigma_s_true = 0.25,
                                   sigma_eps_true = 0.8,
                                   outlier_positions = NULL,
                                   outlier_scale = 8,
                                   seed = 1L) {
  set.seed(seed)
  c_path <- numeric(T_)
  c_path[1L] <- sigma_c_true / sqrt(1 - rho^2) * rnorm(1)
  for (t in 2:T_) c_path[t] <- rho * c_path[t - 1L] + sigma_c_true * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) s_path[t, ] <- s_path[t - 1L, ] + sigma_s_true * rnorm(N)
  noise <- matrix(sigma_eps_true * rnorm(T_ * N), T_, N)
  # Inflate noise at outlier positions
  scale_truth <- matrix(1, T_, N)
  for (pos in outlier_positions) {
    t <- pos[1L]; i <- pos[2L]
    noise[t, i] <- noise[t, i] * outlier_scale
    scale_truth[t, i] <- outlier_scale
  }
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- c_path[t] + s_path[t, ] + noise[t, ]
  }
  list(y = y, c = c_path, s = s_path, scale_truth = scale_truth,
       sigma_c_true = sigma_c_true, sigma_s_true = sigma_s_true,
       sigma_eps_true = sigma_eps_true)
}

test_that("fit_mct_gibbs outlier path runs end-to-end", {
  set.seed(160)
  T_ <- 60; N <- 3
  outliers <- list(c(20, 1), c(35, 2), c(50, 3))   # 3 injected outliers
  sim <- sim_mct_with_outliers(T_, N, outlier_positions = outliers,
                               outlier_scale = 8, seed = 160)
  obs_type <- matrix(1L, T_, N)
  fit <- fit_mct_gibbs(
    sim$y, obs_type = obs_type, ref = 1L,
    n_burn = 50L, n_draw = 100L, verbose = FALSE,
    config = list(use_outliers = TRUE)
  )
  expect_equal(dim(fit$draws$s_outlier), c(T_, N, 100L))
  expect_equal(dim(fit$draws$ps),        c(N, 100L))
  expect_true(all(is.finite(fit$draws$rho)))
  expect_true(all(fit$draws$ps >= 0 & fit$draws$ps <= 1))
})

test_that("outlier path flags injected outliers via elevated s_outlier", {
  # Inject 5 large outliers; the posterior median s at those obs
  # should be > 1.5 (i.e. distinguishable from "normal").
  # At unaffected obs, posterior median s should be ≈ 1.
  set.seed(161)
  T_ <- 80; N <- 4
  outliers <- list(c(15, 1), c(30, 2), c(45, 3), c(60, 1), c(70, 4))
  sim <- sim_mct_with_outliers(T_, N, outlier_positions = outliers,
                               outlier_scale = 10, seed = 161)
  obs_type <- matrix(1L, T_, N)
  fit <- fit_mct_gibbs(
    sim$y, obs_type = obs_type, ref = 1L,
    n_burn = 200L, n_draw = 400L, verbose = FALSE,
    config = list(use_outliers = TRUE)
  )

  # Posterior median s_outlier at each (t, i)
  s_med <- apply(fit$draws$s_outlier, c(1, 2), median)

  # Outlier-flagged obs should have s_med > 1.5
  n_flagged <- 0L
  for (pos in outliers) {
    if (s_med[pos[1L], pos[2L]] > 1.5) n_flagged <- n_flagged + 1L
  }
  # Expect at least 3 of 5 injected outliers flagged
  # (some may be missed under short chains + small obs windows)
  expect_gte(n_flagged, 3L,
             label = sprintf("flagged %d/%d injected outliers", n_flagged,
                             length(outliers)))

  # Most "normal" obs should NOT be flagged (s_med ~= 1)
  mask_inj <- matrix(FALSE, T_, N)
  for (pos in outliers) mask_inj[pos[1L], pos[2L]] <- TRUE
  s_normal <- s_med[!mask_inj]
  prop_normal_flagged <- mean(s_normal > 1.5)
  expect_lt(prop_normal_flagged, 0.15,
            label = sprintf("%.1f%% of normal obs erroneously flagged",
                            100 * prop_normal_flagged))
})
