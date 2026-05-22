# Unit tests for data_prep.R. The most important one is the encoding sanity
# check from the brief: when we have known monthly inflations and aggregate
# them to a quarterly observation, the encoding (`obs_type=2`, placed at the
# last month of the quarter) must round-trip — i.e. when the Stan model
# averages three consecutive latent months it should recover the quarterly
# value we put in.

library(testthat)
source(file.path(testthat::test_path(), "..", "..", "R", "data_fetch.R"))
source(file.path(testthat::test_path(), "..", "..", "R", "data_prep.R"))

# Stub weights aligned to abs_cpi_groups() for tests that don't care about the
# actual weight values — only the build_stan_data plumbing.
stub_weights <- function() {
  g <- abs_cpi_groups()
  tibble::tibble(group = factor(g, levels = g),
                 weight = rep(1 / length(g), length(g)))
}

test_that("build_inflation_series annualises correctly", {
  df <- tibble::tibble(
    date = seq(as.Date("2020-01-01"), as.Date("2020-04-01"), by = "month"),
    group = factor(rep("Housing", 4)),
    index = c(100, 101, 102.01, 103.0301),  # +1% per month
    source = "monthly"
  )
  out <- build_inflation_series(df)
  # log(1.01) * 12 * 100 = ~11.94 (annualised monthly inflation, percent)
  expect_equal(unique(round(out$inflation, 4)), round(log(1.01) * 12 * 100, 4))
})

test_that("demean_series subtracts the 2000-2019 quarterly mean", {
  # Two quarterly obs in the window, one outside
  df <- tibble::tibble(
    date = c(as.Date("2010-03-01"), as.Date("2010-06-01"), as.Date("2025-03-01")),
    group = factor(rep("Housing", 3)),
    inflation = c(2, 4, 10),
    source = rep("quarterly", 3)
  )
  out <- demean_series(df)
  # In-window mean = (2+4)/2 = 3; demeaned values: -1, 1, 7
  expect_equal(out$inflation, c(-1, 1, 7))
  expect_equal(unique(out$series_mean), 3)
})

test_that("build_stan_data places quarterly obs at the last month of the quarter, and monthly obs from the cutover", {
  groups <- abs_cpi_groups()
  # Two quarters of quarterly data (Mar 2024, Jun 2024) for Housing only, then
  # three months of monthly data (Apr-Jun 2024) — at and after the cutover the
  # monthly observations should take precedence.
  df_q <- tibble::tibble(
    date = c(as.Date("2024-03-01"), as.Date("2024-06-01")),
    group = factor("Housing", levels = groups),
    inflation = c(3.0, 4.0),
    source = "quarterly",
    series_mean = 0
  )
  df_m <- tibble::tibble(
    date = seq(as.Date("2024-04-01"), as.Date("2024-06-01"), by = "month"),
    group = factor("Housing", levels = groups),
    inflation = c(2.5, 5.0, 4.5),
    source = "monthly",
    series_mean = 0
  )
  # Add zeros for the other 10 groups at the same dates so build_stan_data has
  # 11 group levels to align to (factor preserves levels even if no rows).
  infl <- dplyr::bind_rows(df_q, df_m)

  sd_ <- build_stan_data(infl,
                         weights = stub_weights(),
                         start = as.Date("2024-01-01"),
                         monthly_cutover = as.Date("2024-04-01"))

  # Index of Housing
  h <- match("Housing", sd_$groups)

  # Mar 2024 (pre-cutover) -> quarterly obs, value 3.0
  t_mar <- match(as.Date("2024-03-01"), sd_$dates)
  expect_equal(unname(sd_$obs_type[t_mar, h]), 2L)
  expect_equal(unname(sd_$y[t_mar, h]), 3.0)

  # Apr 2024 -> monthly obs (cutover), value 2.5
  t_apr <- match(as.Date("2024-04-01"), sd_$dates)
  expect_equal(unname(sd_$obs_type[t_apr, h]), 1L)
  expect_equal(unname(sd_$y[t_apr, h]), 2.5)

  # Jun 2024 -> monthly obs (NOT the quarterly 4.0). The quarterly obs at
  # Jun 2024 should be ignored because monthly takes precedence post-cutover.
  t_jun <- match(as.Date("2024-06-01"), sd_$dates)
  expect_equal(unname(sd_$obs_type[t_jun, h]), 1L)
  expect_equal(unname(sd_$y[t_jun, h]), 4.5)

  # Jan, Feb, May 2024 for Housing -> unobserved (no monthly data for Jan/Feb,
  # no quarterly placement; May only has monthly which IS present).
  t_jan <- match(as.Date("2024-01-01"), sd_$dates)
  t_feb <- match(as.Date("2024-02-01"), sd_$dates)
  expect_equal(unname(sd_$obs_type[t_jan, h]), 0L)
  expect_equal(unname(sd_$obs_type[t_feb, h]), 0L)
})

test_that("build_stan_data emits correctly-shaped flat arrays for Stan", {
  groups <- abs_cpi_groups()
  df_q <- tibble::tibble(
    date = c(as.Date("2024-03-01"), as.Date("2023-12-01")),
    group = factor("Housing", levels = groups),
    inflation = c(3.0, 2.0),
    source = "quarterly",
    series_mean = 0
  )
  df_m <- tibble::tibble(
    date = as.Date("2024-04-01"),
    group = factor("Housing", levels = groups),
    inflation = 2.5,
    source = "monthly",
    series_mean = 0
  )
  infl <- dplyr::bind_rows(df_q, df_m)
  sd_ <- build_stan_data(infl,
                         weights = stub_weights(),
                         start = as.Date("2024-01-01"),
                         monthly_cutover = as.Date("2024-04-01"))

  # T = 4 months (Jan-Apr 2024); the quarterly obs at 2023-12 falls before
  # `start` and should be dropped.
  expect_equal(sd_$T, 4L)
  expect_equal(sd_$N, 11L)
  expect_equal(sd_$ref, match("Housing", sd_$groups))
  expect_equal(length(sd_$w), 11L)
  expect_equal(sum(sd_$w), 1, tolerance = 1e-12)

  # Monthly observation: 1 obs (Housing at t=4)
  expect_equal(sd_$n_obs_m, 1L)
  expect_equal(sd_$t_m, 4L)
  expect_equal(sd_$i_m, match("Housing", sd_$groups))
  expect_equal(sd_$y_m, 2.5)

  # Quarterly observation: would be at t=3 (Mar 2024); needs t >= 3 to be kept.
  expect_equal(sd_$n_obs_q, 1L)
  expect_equal(sd_$t_q, 3L)
  expect_equal(sd_$i_q, match("Housing", sd_$groups))
  expect_equal(sd_$y_q, 3.0)
})

test_that("3-month average of latent monthly inflations recovers the quarterly observation by construction (encoding sanity)", {
  # Simulate the Stan model's view: invent N=11 series of latent monthly
  # annualised inflations on the 1990-2024 clock, aggregate to quarterly via
  # the simple average (the model's measurement equation), feed the quarterly
  # series into build_stan_data, and confirm that for every observed
  # quarterly cell the y value equals (pi[t-2] + pi[t-1] + pi[t]) / 3.

  set.seed(42)
  groups <- abs_cpi_groups()
  dates_m <- seq(as.Date("2020-01-01"), as.Date("2023-12-01"), by = "month")
  pi_latent <- matrix(
    rnorm(length(dates_m) * length(groups), mean = 0, sd = 2),
    nrow = length(dates_m), ncol = length(groups),
    dimnames = list(NULL, groups)
  )

  # Build the quarterly observation as the average of 3 consecutive months,
  # placed at the last month of each quarter (Mar, Jun, Sep, Dec).
  is_qend <- lubridate::month(dates_m) %in% c(3, 6, 9, 12)
  q_idx <- which(is_qend)
  q_idx <- q_idx[q_idx >= 3]  # need a full 3-month window
  q_avg <- t(vapply(q_idx, function(t) colMeans(pi_latent[(t - 2):t, ]), numeric(length(groups))))

  infl_q <- tibble::tibble(
    date = rep(dates_m[q_idx], length(groups)),
    group = factor(rep(groups, each = length(q_idx)), levels = groups),
    inflation = as.vector(q_avg),
    source = "quarterly",
    series_mean = 0
  )

  sd_ <- build_stan_data(infl_q,
                         weights = stub_weights(),
                         start = as.Date("2020-01-01"),
                         monthly_cutover = as.Date("2099-01-01"))  # disable monthly path

  # For each quarter-end cell with obs_type == 2, confirm y matches the avg of
  # latent months t-2, t-1, t.
  for (t in q_idx) {
    expected <- colMeans(pi_latent[(t - 2):t, , drop = FALSE])
    for (n in seq_along(groups)) {
      expect_equal(unname(sd_$obs_type[t, n]), 2L)
      expect_equal(unname(sd_$y[t, n]), unname(expected[n]),
                   tolerance = 1e-12,
                   info = paste0("t=", dates_m[t], " group=", groups[n]))
    }
  }
})
