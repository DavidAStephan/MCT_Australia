# Tests for Variant B on the non-mixed-frequency Gibbs path
# (quarterly-clocked sub-group setting). Validates:
#   - the non-mixed gibbs_sweep accepts variant = "B"
#   - rho stays pinned at 1 throughout
#   - sigma_c posterior recovers the true RW innovation SD
#   - NA cells in y are handled (staggered-start sectors)

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate Variant B on a single (quarterly) clock — no mixed-freq.
# Some sectors have NA observations at the start (staggered series).
sim_b_quarterly <- function(T_, N, sigma_c_true = 0.4, sigma_s_true = 0.25,
                            sigma_eps_true = 0.6, ref = 1L,
                            staggered_start = NULL, seed = 1L) {
  set.seed(seed)
  c_path <- numeric(T_)
  c_path[1L] <- rnorm(1, 0, 2)            # diffuse Variant B prior
  for (t in 2:T_) c_path[t] <- c_path[t - 1L] + sigma_c_true * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) {
    s_path[t, ] <- s_path[t - 1L, ] + sigma_s_true * rnorm(N)
  }
  y <- matrix(NA_real_, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- c_path[t] + s_path[t, ] + sigma_eps_true * rnorm(N)
  }
  # Apply staggered starts: NA out the first `staggered_start[i]` obs of sector i
  if (!is.null(staggered_start)) {
    for (i in seq_len(N)) {
      k <- staggered_start[i]
      if (k > 0) y[seq_len(k), i] <- NA_real_
    }
  }
  list(y = y, c = c_path, s = s_path,
       sigma_c_true = sigma_c_true)
}

test_that("Variant B non-mixed path runs and pins rho at 1", {
  set.seed(800)
  T_ <- 80; N <- 3
  sim <- sim_b_quarterly(T_, N, seed = 800)
  fit <- fit_mct_gibbs(
    sim$y, ref = 1L,
    n_burn = 100L, n_draw = 200L, verbose = FALSE,
    config = list(variant = "B")
  )
  expect_equal(dim(fit$draws$c), c(T_, 200L))
  expect_true(all(fit$draws$rho == 1),
              label = "Variant B rho pinned at 1")
})

test_that("Variant B quarterly NA-handling works (staggered starts)", {
  set.seed(801)
  T_ <- 100; N <- 4
  # Sector 4 starts 30 quarters in (mimicking 'Financial services')
  sim <- sim_b_quarterly(T_, N, sigma_c_true = 0.3,
                         sigma_s_true = 0.2, sigma_eps_true = 0.5,
                         staggered_start = c(0, 0, 0, 30), seed = 801)
  fit <- fit_mct_gibbs(
    sim$y, ref = 1L,
    n_burn = 200L, n_draw = 400L, verbose = FALSE,
    config = list(variant = "B")
  )
  expect_true(all(is.finite(fit$draws$c)))
  expect_true(all(is.finite(fit$draws$s)),
              label = "no NaN/Inf in sector-trend draws despite staggered NAs")
  # sigma_c posterior median within 50% of truth
  sigma_c_med <- median(fit$draws$sigma_c)
  expect_gt(sigma_c_med, 0.5 * sim$sigma_c_true)
  expect_lt(sigma_c_med, 2.0 * sim$sigma_c_true)
})
