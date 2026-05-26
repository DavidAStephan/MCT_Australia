# Tests for the mct_gibbs_fit S3 wrapper. Validates that:
#   - $draws(variables, format) returns the expected shapes + names
#   - $diagnostic_summary() and $save_object() work
#   - the wrapper composes with extract_path_matrix() and trend_path()
#     from R/postprocess.R (i.e. downstream code is unchanged)

library(testthat)
suppressMessages({
  library(MASS)
  library(posterior)
})

for (f in list.files("../../../R/gibbs", "\\.R$", full.names = TRUE)) {
  source(f)
}
# Also source postprocess.R so we can test composition
source("../../../R/postprocess.R")

# Inline simulator (same as in test-fit-mct-gibbs.R). Inlining avoids
# the test_that-bodies-fire-on-source problem in cross-file sourcing.
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
  list(y = y, c = c_path, s = s_path, rho = rho, lambda = lambda,
       sigma_c = sigma_c, sigma_s = sigma_s, sigma_eps = sigma_eps)
}

# Helper: build a minimal stan_data list compatible with the wrapper.
.make_stan_data <- function(y, dates = NULL, groups = NULL, w = NULL) {
  T_ <- nrow(y); N <- ncol(y)
  if (is.null(dates))  dates  <- seq(as.Date("2000-01-01"),
                                     by = "month", length.out = T_)
  if (is.null(groups)) groups <- paste0("g", seq_len(N))
  if (is.null(w))      w      <- rep(1 / N, N)
  # Monthly-only obs: every (t, i) observed
  t_m <- rep(seq_len(T_), each = N)
  i_m <- rep(seq_len(N), times = T_)
  y_m <- as.numeric(t(y))                  # row-major over (t, i)
  list(
    T = T_, N = N, dates = dates, groups = groups, w = w,
    n_obs_m = length(y_m),
    t_m = as.integer(t_m), i_m = as.integer(i_m), y_m = y_m,
    n_obs_q = 0L,
    t_q = integer(0), i_q = integer(0), y_q = numeric(0)
  )
}

test_that("mct_gibbs_fit builds and exposes correct methods", {
  set.seed(60)
  sim <- sim_mct_a(T_ = 30, N = 3, seed = 60)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  sd_obj <- .make_stan_data(sim$y)
  wrapped <- mct_gibbs_fit(fit, sd_obj)

  expect_s3_class(wrapped, "mct_gibbs_fit")
  expect_true(is.function(wrapped$draws))
  expect_true(is.function(wrapped$diagnostic_summary))
  expect_true(is.function(wrapped$save_object))

  diag <- wrapped$diagnostic_summary()
  expect_equal(diag$num_divergent, 0L)
  expect_equal(diag$ebfmi, 1.0)
})

test_that("$draws('rho') returns correct shape + variable name", {
  set.seed(61)
  sim <- sim_mct_a(T_ = 30, N = 3, seed = 61)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  wrapped <- mct_gibbs_fit(fit, .make_stan_data(sim$y))
  m <- wrapped$draws(variables = "rho", format = "draws_matrix")
  expect_equal(nrow(m), 100L)
  expect_equal(ncol(m), 1L)
  expect_equal(colnames(m), "rho")
})

test_that("$draws('trend') returns T columns named trend[1]..trend[T]", {
  set.seed(62)
  T_ <- 30
  sim <- sim_mct_a(T_ = T_, N = 3, seed = 62)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  wrapped <- mct_gibbs_fit(fit, .make_stan_data(sim$y))
  m <- wrapped$draws(variables = "trend")
  expect_equal(ncol(m), T_)
  expect_equal(colnames(m)[1L], "trend[1]")
  expect_equal(colnames(m)[T_], sprintf("trend[%d]", T_))
})

test_that("$draws('sector_trend') returns T*N columns col-major", {
  set.seed(63)
  T_ <- 20; N <- 3
  sim <- sim_mct_a(T_ = T_, N = N, seed = 63)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  wrapped <- mct_gibbs_fit(fit, .make_stan_data(sim$y))
  m <- wrapped$draws(variables = "sector_trend")
  expect_equal(ncol(m), T_ * N)
  # Column-major: first T columns are "[t,1]" for t=1..T, then "[t,2]" etc.
  expect_equal(colnames(m)[1L], "sector_trend[1,1]")
  expect_equal(colnames(m)[T_], sprintf("sector_trend[%d,1]", T_))
  expect_equal(colnames(m)[T_ + 1L], "sector_trend[1,2]")
})

test_that("composes with extract_path_matrix from postprocess.R", {
  set.seed(64)
  sim <- sim_mct_a(T_ = 30, N = 3, seed = 64)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  wrapped <- mct_gibbs_fit(fit, .make_stan_data(sim$y))
  m <- extract_path_matrix(wrapped, "trend")
  expect_equal(dim(m), c(100L, 30L))
})

test_that("composes with trend_path from postprocess.R", {
  set.seed(65)
  T_ <- 30
  sim <- sim_mct_a(T_ = T_, N = 3, seed = 65)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  sd_obj <- .make_stan_data(sim$y)
  wrapped <- mct_gibbs_fit(fit, sd_obj)
  tr <- trend_path(wrapped, sd_obj$dates)
  expect_equal(nrow(tr), T_)
  expect_named(tr, c("date", "lower", "median", "upper"))
  expect_true(all(is.finite(tr$median)))
})

test_that("$save_object roundtrip preserves the wrapper", {
  set.seed(66)
  sim <- sim_mct_a(T_ = 20, N = 3, seed = 66)
  fit <- fit_mct_gibbs(sim$y, ref = 1L,
                       n_burn = 50L, n_draw = 100L, verbose = FALSE)
  wrapped <- mct_gibbs_fit(fit, .make_stan_data(sim$y))
  tmp <- tempfile(fileext = ".rds")
  wrapped$save_object(tmp)
  reloaded <- readRDS(tmp)
  expect_s3_class(reloaded, "mct_gibbs_fit")
  expect_equal(reloaded$draws("rho"), wrapped$draws("rho"))
  unlink(tmp)
})
