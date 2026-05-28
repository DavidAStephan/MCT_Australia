# Tests for the outlier scale mixture on the non-mixed-frequency Gibbs
# path (quarterly-clocked panels). Mirrors the existing
# test-fit-mct-gibbs-outliers.R (mixed-freq monthly-only path) but
# exercises the new code in R/gibbs/gibbs_sweep.R.

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate Variant B on a quarterly clock with INJECTED outliers at
# specific (t, i) positions by inflating the obs noise.
sim_b_q_with_outliers <- function(T_, N, sigma_c_true = 0.4,
                                  sigma_s_true = 0.25,
                                  sigma_eps_true = 0.6,
                                  outlier_positions = NULL,
                                  outlier_scale = 8, seed = 1L) {
  set.seed(seed)
  c_path <- numeric(T_); c_path[1L] <- rnorm(1, 0, 2)
  for (t in 2:T_) c_path[t] <- c_path[t - 1L] + sigma_c_true * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) s_path[t, ] <- s_path[t - 1L, ] + sigma_s_true * rnorm(N)
  noise <- matrix(sigma_eps_true * rnorm(T_ * N), T_, N)
  scale_truth <- matrix(1, T_, N)
  for (pos in outlier_positions) {
    noise[pos[1L], pos[2L]] <- noise[pos[1L], pos[2L]] * outlier_scale
    scale_truth[pos[1L], pos[2L]] <- outlier_scale
  }
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- c_path[t] + s_path[t, ] + noise[t, ]
  }
  list(y = y, scale_truth = scale_truth)
}

test_that("quarterly Variant B with outliers runs end-to-end", {
  set.seed(900)
  T_ <- 80; N <- 4
  outliers <- list(c(20, 1), c(45, 3))
  sim <- sim_b_q_with_outliers(T_, N, outlier_positions = outliers,
                               outlier_scale = 8, seed = 900)
  fit <- fit_mct_gibbs(
    sim$y, ref = 1L,
    n_burn = 100L, n_draw = 200L, verbose = FALSE,
    config = list(variant = "B", use_outliers = TRUE)
  )
  expect_equal(dim(fit$draws$s_outlier), c(T_, N, 200L))
  expect_equal(dim(fit$draws$ps),         c(N, 200L))
  expect_true(all(fit$draws$rho == 1))
  expect_true(all(fit$draws$ps >= 0 & fit$draws$ps <= 1))
})

test_that("outlier path flags injected outliers at quarterly clock", {
  set.seed(901)
  T_ <- 100; N <- 5
  outliers <- list(c(25, 1), c(40, 2), c(55, 3), c(70, 4), c(85, 5))
  sim <- sim_b_q_with_outliers(T_, N, outlier_positions = outliers,
                               outlier_scale = 10, seed = 901)
  fit <- fit_mct_gibbs(
    sim$y, ref = 1L,
    n_burn = 300L, n_draw = 600L, verbose = FALSE,
    config = list(variant = "B", use_outliers = TRUE)
  )

  s_med <- apply(fit$draws$s_outlier, c(1, 2), median)

  # Expect at least 2 of 5 injected outliers flagged. Identification
  # is intentionally looser here than in the AR(1) common-factor test
  # because Variant B's RW common factor can absorb point spikes into
  # a single-period trend movement, leaving less posterior mass on
  # large s_outlier values. This is the same dynamic at play with
  # COVID-period AU data — we expect outliers to help but not flag
  # everything we'd expect a priori.
  n_flagged <- sum(vapply(outliers,
                          function(pos) s_med[pos[1L], pos[2L]] > 1.5,
                          logical(1)))
  expect_gte(n_flagged, 2L,
             label = sprintf("flagged %d/%d injected", n_flagged,
                             length(outliers)))

  # Most "normal" obs should NOT be flagged
  mask_inj <- matrix(FALSE, T_, N)
  for (pos in outliers) mask_inj[pos[1L], pos[2L]] <- TRUE
  prop_normal_flagged <- mean(s_med[!mask_inj] > 1.5)
  expect_lt(prop_normal_flagged, 0.15,
            label = sprintf("%.1f%% normal obs flagged",
                            100 * prop_normal_flagged))
})

test_that("outliers + NA handling work together (staggered sectors)", {
  set.seed(902)
  T_ <- 80; N <- 4
  outliers <- list(c(40, 1))
  sim <- sim_b_q_with_outliers(T_, N, outlier_positions = outliers,
                               seed = 902)
  # Sector 4 starts at t = 30
  sim$y[1:29, 4L] <- NA_real_
  fit <- fit_mct_gibbs(
    sim$y, ref = 1L,
    n_burn = 200L, n_draw = 400L, verbose = FALSE,
    config = list(variant = "B", use_outliers = TRUE)
  )
  expect_true(all(is.finite(fit$draws$c)))
  expect_true(all(is.finite(fit$draws$s)))
  expect_true(all(is.finite(fit$draws$s_outlier)))
})
