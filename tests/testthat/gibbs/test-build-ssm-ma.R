# Test the MA-augmented SSM constructor by simulating from it and
# checking that:
#   - The state has the expected dim (1 + N + (q+1)*N)
#   - Empirical c-path has the right AR(1) persistence
#   - Empirical residuals exhibit MA(q) autocorrelation structure

library(testthat)
suppressMessages(library(MASS))

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

test_that("build_mct_ssm_ma constructs correct-dim SSM", {
  T_ <- 50; N <- 3; q_MA <- 3
  theta <- matrix(c(0.5, -0.3, 0.1,
                    0.4, -0.2, 0.05,
                    0.3,  0.1, 0.0), N, q_MA, byrow = TRUE)
  ssm <- build_mct_ssm_ma(
    rho = 0.7, lambda = rep(1, N),
    sigma_c = rep(0.4, T_),
    sigma_s = matrix(0.25, T_, N),
    sigma_eps = matrix(0.8, T_, N),
    theta = theta, T_ = T_, N = N, q_MA = q_MA
  )
  expected_D <- 1L + N + (q_MA + 1L) * N
  expect_equal(nrow(ssm$F), expected_D)
  expect_equal(ncol(ssm$F), expected_D)
  expect_equal(dim(ssm$H), c(N, expected_D, T_))
})

test_that("MA-augmented simulator gives AR(1) c-path with right rho", {
  set.seed(120)
  T_ <- 1000; N <- 3; q_MA <- 3
  rho_true <- 0.7
  theta <- matrix(0, N, q_MA)
  ssm <- build_mct_ssm_ma(
    rho = rho_true, lambda = rep(1, N),
    sigma_c = rep(0.4, T_),
    sigma_s = matrix(0.25, T_, N),
    sigma_eps = matrix(0.8, T_, N),
    theta = theta, T_ = T_, N = N, q_MA = q_MA
  )
  sim <- simulate_ssm(ssm, T_ = T_, need_states = TRUE)
  c_path <- as.numeric(sim$states[ssm$.meta$c_pos, ])
  rho_emp <- cor(c_path[-1], c_path[-T_])
  expect_lt(abs(rho_emp - rho_true), 0.1)
})

test_that("MA-augmented obs residuals exhibit the expected MA(q) ACF", {
  # With theta_1 = 0.5, the residuals should have lag-1 ACF ≈
  # theta_1 / (1 + theta_1^2 + theta_2^2 + theta_3^2). For pure MA(1)
  # with theta_1=0.5: lag-1 ACF = 0.5 / 1.25 = 0.4.
  set.seed(121)
  T_ <- 2000; N <- 1; q_MA <- 1
  theta <- matrix(0.5, N, q_MA)
  ssm <- build_mct_ssm_ma(
    rho = 0, lambda = rep(1, N),
    sigma_c = rep(1e-3, T_),   # near-zero common factor
    sigma_s = matrix(1e-3, T_, N),   # near-zero sector trend
    sigma_eps = matrix(1, T_, N),
    theta = theta, T_ = T_, N = N, q_MA = q_MA,
    var_s_init = 1e-3
  )
  sim <- simulate_ssm(ssm, T_ = T_, need_states = TRUE)
  y <- as.numeric(sim$Y)
  # Theoretical lag-1 ACF for MA(1) with theta=0.5: 0.5 / 1.25 = 0.4
  lag1_acf <- acf(y, lag.max = 1, plot = FALSE)$acf[2]
  expect_lt(abs(lag1_acf - 0.4), 0.05,
            label = sprintf("MA(1) theta=0.5 ⇒ lag1 ACF empirical %.3f vs theory 0.4",
                            lag1_acf))
})
