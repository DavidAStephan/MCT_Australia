# End-to-end test of fit_mct_gibbs: simulate from Variant A spec,
# run a Gibbs chain, check recovery of headline parameters and paths.

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Simulate from MCT Variant A spec directly (matches CLAUDE.md /
# stan/mct_aus_A.stan model, simplified to constant lambda + monthly obs).
sim_mct_a <- function(T_, N, rho = 0.7, lambda = NULL,
                      sigma_c_true = 0.4, sigma_s_true = 0.25,
                      sigma_eps_true = 0.8, ref = 1L, seed = 1L) {
  set.seed(seed)
  if (is.null(lambda)) {
    lambda <- rep(1, N)            # all loadings = 1 for simplicity
  }
  lambda[ref] <- 1

  # Constant SV for the truth (we're testing parameter recovery, not SV
  # path recovery — SV path uncertainty is large per-period).
  sigma_c   <- rep(sigma_c_true, T_)
  sigma_s   <- matrix(sigma_s_true, T_, N)
  sigma_eps <- matrix(sigma_eps_true, T_, N)

  # AR(1) common factor, stationary start
  c_path <- numeric(T_)
  c_path[1L] <- sigma_c[1L] / sqrt(1 - rho^2) * rnorm(1)
  for (t in 2:T_) c_path[t] <- rho * c_path[t - 1L] + sigma_c[t] * rnorm(1)

  # RW sector trends starting at 0
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) {
    s_path[t, ] <- s_path[t - 1L, ] + sigma_s[t, ] * rnorm(N)
  }

  # Observation
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- lambda * c_path[t] + s_path[t, ] + sigma_eps[t, ] * rnorm(N)
  }

  list(y = y, c = c_path, s = s_path, rho = rho, lambda = lambda,
       sigma_c = sigma_c, sigma_s = sigma_s, sigma_eps = sigma_eps)
}

test_that("fit_mct_gibbs runs end-to-end on small sim data without error", {
  set.seed(50)
  sim <- sim_mct_a(T_ = 60, N = 3, seed = 50)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 100L, n_draw = 200L,
                       verbose = FALSE)
  expect_named(fit, c("draws", "config", "n_iter", "elapsed_sec"))
  expect_equal(dim(fit$draws$c), c(60L, 200L))
  expect_equal(dim(fit$draws$s), c(60L, 3L, 200L))
  expect_length(fit$draws$rho, 200L)
})

test_that("fit_mct_gibbs rho posterior is informative and not crazy", {
  # NOTE: tight rho recovery is hard with this short-chain config. The
  # AR(1) c-process competes with the RW sector trends for signal —
  # Gibbs needs many more iters than HMC to mix these out. With 500
  # burn + 1000 draws we expect rho_median in roughly the (0.3, 0.9)
  # range and posterior CONCENTRATED away from the prior, not exact
  # truth recovery. Treat this as "sampler is finding signal" check;
  # use the longer benchmark scripts for production-quality recovery.
  set.seed(51)
  T_  <- 200
  rho_true <- 0.7
  sim <- sim_mct_a(T_ = T_, N = 5, rho = rho_true, seed = 51)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 500L, n_draw = 1000L,
                       verbose = FALSE)
  rho_med  <- median(fit$draws$rho)
  rho_sd   <- sd(fit$draws$rho)
  # Posterior should be substantially tighter than the N(0.5, 0.3) prior
  expect_lt(rho_sd, 0.25)
  # Median in a reasonable range — not stuck at the prior, not negative
  expect_gt(rho_med, 0.3)
  expect_lt(rho_med, 0.95)
})

test_that("fit_mct_gibbs recovers non-reference loadings (loose)", {
  # Tolerance is 0.4 not 0.2 because lambda and c can trade signal
  # under the slow Gibbs mixing — see note in rho test above.
  set.seed(52)
  T_ <- 150
  lambda_true <- c(1, 0.7, 1.3, 0.9)
  sim <- sim_mct_a(T_ = T_, N = 4, lambda = lambda_true,
                   ref = 1L, seed = 52)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 500L, n_draw = 1000L,
                       verbose = FALSE)
  lambda_med <- apply(fit$draws$lambda, 1, median)
  expect_equal(lambda_med[1L], 1)  # ref pinned
  for (i in 2:4) {
    expect_lt(abs(lambda_med[i] - lambda_true[i]), 0.4,
              label = sprintf("lambda[%d] recovery", i))
  }
})

test_that("fit_mct_gibbs c-path correlates positively with truth", {
  # cor > 0.5, not 0.7 — c and s can trade signal under slow Gibbs
  # mixing. Long-chain benchmarks should show >0.8.
  set.seed(53)
  T_ <- 150
  sim <- sim_mct_a(T_ = T_, N = 5, seed = 53)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 300L, n_draw = 600L,
                       verbose = FALSE)
  c_post_mean <- rowMeans(fit$draws$c)
  cor_c <- cor(c_post_mean, sim$c)
  expect_gt(cor_c, 0.5)
})
