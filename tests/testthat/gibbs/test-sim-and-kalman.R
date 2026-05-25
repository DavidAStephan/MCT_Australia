# Sanity tests for the ported simulator + Kalman filter. Validates
# the two together on a tiny scalar AR(1) state-space model where the
# log-likelihood and innovations can be checked analytically.

library(testthat)
suppressMessages(library(MASS))

# Source the gibbs ports
for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

# Manual multivariate-normal log-density (avoids the mvtnorm dependency).
.dmvnorm_log <- function(x, mean, sigma) {
  d  <- length(x)
  L  <- chol(sigma)                       # upper triangular
  z  <- backsolve(L, x - mean, transpose = TRUE)
  -0.5 * (d * log(2 * pi) +
          2 * sum(log(diag(L))) +         # log det Sigma = 2 * sum log diag(L)
          sum(z * z))
}

# ---------------------------------------------------------------------
# Build a scalar AR(1) + Gaussian-noise observation model:
#   x_{t+1} = phi * x_t + eta_t,  eta ~ N(0, sigma_eta^2)
#   y_t     =        x_t + eps_t,  eps ~ N(0, sigma_eps^2)
#   x_1     ~ N(0, sigma_eta^2 / (1 - phi^2))   (stationary)
# ---------------------------------------------------------------------
make_ar1_ssm <- function(phi = 0.7, sigma_eta = 0.5, sigma_eps = 1.0) {
  list(
    H        = matrix(1, 1, 1),
    F        = matrix(phi, 1, 1),
    G        = matrix(1, 1, 1),
    Sigma_eps = matrix(sigma_eps^2, 1, 1),
    Sigma_eta = matrix(sigma_eta^2, 1, 1),
    mu_1     = 0,
    Sigma_1  = matrix(sigma_eta^2 / (1 - phi^2), 1, 1)
  )
}

test_that("simulate_ssm produces output of the right shape", {
  set.seed(1)
  ssm <- make_ar1_ssm()
  out <- simulate_ssm(ssm, T_ = 50,
                     need_states = TRUE,
                     need_disturbances = TRUE)
  expect_equal(dim(out$Y), c(1L, 50L))
  expect_equal(dim(out$states), c(1L, 50L))
  expect_equal(dim(out$disturbances$m_errors), c(1L, 50L))
  expect_equal(dim(out$disturbances$shocks),   c(1L, 49L))
})

test_that("simulated AR(1) data has approximately the right marginal", {
  # Stationary variance: Var(x) = sigma_eta^2 / (1 - phi^2);
  # Var(y) = Var(x) + sigma_eps^2.
  set.seed(2)
  phi <- 0.7; se <- 0.5; sm <- 1.0
  ssm <- make_ar1_ssm(phi, se, sm)
  # Many reps to average over
  ys <- replicate(50, simulate_ssm(ssm, T_ = 500)$Y)
  # ys is 1 x 500 x 50; flatten each replicate's series before var()
  var_y_sim <- mean(apply(ys, 3, function(m) var(as.numeric(m))))
  var_x_theory <- se^2 / (1 - phi^2)
  var_y_theory <- var_x_theory + sm^2
  # Allow ~10% Monte Carlo tolerance
  expect_lt(abs(var_y_sim / var_y_theory - 1), 0.10)
})

test_that("kalman_filter log-likelihood matches dmvnorm of marginal y", {
  # For a fixed simulated y, the SSM defines a joint Gaussian on
  # (y_1, ..., y_T). The Kalman filter log-likelihood should equal
  # the joint log-density.
  set.seed(3)
  ssm <- make_ar1_ssm(phi = 0.5, sigma_eta = 0.4, sigma_eps = 0.8)
  T_ <- 10
  sim <- simulate_ssm(ssm, T_ = T_)
  y <- as.numeric(sim$Y)

  # Build the joint covariance of (y_1, ..., y_T) by hand for a scalar
  # AR(1) + iid noise: y_t = x_t + eps_t, x_t = sum_{k>=0} phi^k eta_{t-k}.
  phi <- 0.5; se2 <- 0.4^2; sm2 <- 0.8^2
  var_x <- se2 / (1 - phi^2)
  cov_x <- outer(seq_len(T_), seq_len(T_),
                 Vectorize(function(s, t) phi^abs(s - t) * var_x))
  cov_y <- cov_x + sm2 * diag(T_)
  ll_truth <- .dmvnorm_log(y, mean = rep(0, T_), sigma = cov_y)

  ll_kf <- kalman_filter(matrix(y, 1, T_), ssm)$log_likelihood
  expect_equal(ll_kf, ll_truth, tolerance = 1e-8)
})

test_that("kalman_filter handles missing data by skipping the update", {
  # If we set y_3 to NA, the filter should produce the same log-lik
  # as if we'd analytically integrated y_3 out (i.e. it doesn't enter
  # the joint density at all).
  set.seed(4)
  ssm <- make_ar1_ssm(phi = 0.6, sigma_eta = 0.5, sigma_eps = 1.0)
  T_ <- 8
  y_full <- as.numeric(simulate_ssm(ssm, T_ = T_)$Y)
  y_miss <- y_full
  y_miss[3] <- NA

  # Drop position 3 from the analytic marginal
  phi <- 0.6; se2 <- 0.5^2; sm2 <- 1.0^2
  var_x <- se2 / (1 - phi^2)
  keep <- c(1, 2, 4, 5, 6, 7, 8)
  cov_x <- outer(keep, keep, Vectorize(function(s, t) phi^abs(s - t) * var_x))
  cov_y_keep <- cov_x + sm2 * diag(length(keep))
  ll_truth <- .dmvnorm_log(y_full[keep], rep(0, length(keep)), cov_y_keep)

  ll_kf <- kalman_filter(matrix(y_miss, 1, T_), ssm)$log_likelihood
  expect_equal(ll_kf, ll_truth, tolerance = 1e-8)
})
