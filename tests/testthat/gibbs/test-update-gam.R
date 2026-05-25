# Tests for update_gam (Inverse-Gamma conjugate update for SV step sizes).

library(testthat)

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}

test_that("update_gam runs on vector input and returns positive scalar", {
  set.seed(30)
  x <- rnorm(100, sd = 0.5)
  out <- update_gam(x, nu_prior = 10, s2_prior = 0.25)
  expect_length(out, 1L)
  expect_true(is.finite(out))
  expect_gt(out, 0)
})

test_that("update_gam returns N draws for N-column input", {
  set.seed(31)
  N <- 4
  x <- matrix(rnorm(50 * N, sd = c(0.2, 0.5, 1.0, 2.0)),
              nrow = 50, ncol = N, byrow = TRUE)
  out <- update_gam(x, nu_prior = rep(10, N), s2_prior = rep(0.5, N))
  expect_length(out, N)
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0))
})

test_that("update_gam posterior concentrates around truth as T grows", {
  # With T=10,000 obs from N(0, true_gam^2), the posterior mean of gamma
  # should be very close to true_gam.
  set.seed(32)
  true_gam <- 0.7
  T_ <- 10000
  x <- rnorm(T_, sd = true_gam)
  draws <- replicate(500, update_gam(x, nu_prior = 5, s2_prior = 1))
  post_mean <- mean(draws)
  # Posterior CV at T=10K is ~1/sqrt(2T) ~ 0.7%. Easy 2% tolerance.
  expect_lt(abs(post_mean / true_gam - 1), 0.02)
})

test_that("update_gam respects prior when data is uninformative (T = 0)", {
  # With T=0, posterior should equal prior. Sample many draws and check
  # the implied gamma distribution looks like the prior.
  # Prior: gamma^2 ~ IG(nu_prior/2, nu_prior*s2_prior/2)
  # E[gamma^2] = nu_prior*s2_prior/(nu_prior - 2)  (for nu_prior > 2)
  set.seed(33)
  nu  <- 10
  s2  <- 0.5
  draws <- replicate(2000, update_gam(matrix(0, 0, 1),
                                      nu_prior = nu, s2_prior = s2))
  emp_mean_g2 <- mean(draws^2)
  prior_mean_g2 <- nu * s2 / (nu - 2)
  expect_lt(abs(emp_mean_g2 / prior_mean_g2 - 1), 0.10)
})
