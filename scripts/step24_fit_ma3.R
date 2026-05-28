# Step 24 — Fit Variant B with MA(3) sector noise on the quarterly
# sub-group panel. Replaces (or stacks on top of) the outlier scale
# mixture as a second pass at smoothing the COVID-period trend.
#
# State dim at q=3, N=33: 1 + 33 + 4*33 = 166 (vs the no-MA model's 34).
# Per-iter cost scales roughly cubically with state dim → ~120x slower
# than no-MA, so ~10 hours wall time at n_burn=1000, n_draw=2000.
# Smoke test first to get an accurate budget.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

sd_q <- readRDS("outputs/draws/stan_data.rds")
y    <- sd_q$y

cat(sprintf(
  "[step24] sd_q: T=%d N=%d  n_obs_q=%d  ref=%s\n",
  sd_q$T, sd_q$N, sd_q$n_obs_q, sd_q$groups[sd_q$ref]
))

# ----- Smoke test: 10 iters to get per-iter time ----------------------
cat("\n[step24] Smoke test (10 iters)...\n")
set.seed(99L)
t0 <- Sys.time()
fit_smoke <- fit_mct_gibbs(
  y, q_MA = 3L, ref = sd_q$ref,
  n_burn = 5L, n_draw = 5L, verbose = FALSE,
  config = list(variant = "B",
                s2_gam_prior = 1 / (60 * 4))
)
sec_per_iter <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 10
cat(sprintf(
  "[step24] Smoke per-iter: %.2f sec; projected 3000-iter = %.0f min (%.1f h)\n",
  sec_per_iter, sec_per_iter * 3000 / 60, sec_per_iter * 3000 / 3600
))

# ----- Full fit -------------------------------------------------------
out_path <- "outputs/draws/fit_B_ma3.rds"
if (file.exists(out_path)) {
  cat("\n[step24] cached at ", out_path, " — loading.\n", sep = "")
  fit <- readRDS(out_path)
} else {
  cat("\n[step24] Fitting Variant B + MA(3) ...\n")
  set.seed(20260528L)
  t0 <- Sys.time()
  fit <- fit_mct_gibbs(
    y, q_MA = 3L, ref = sd_q$ref,
    n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
    config = list(variant       = "B",
                  s2_gam_prior  = 1 / (60 * 4))
  )
  cat(sprintf("[step24] done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit, out_path)
  cat("[step24] saved ", out_path, "\n", sep = "")
}

# ----- Comparison vs no-MA (current production) ----------------------
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

fit_curr <- readRDS("outputs/draws/fit_B.rds")
fit_curr_wrap <- mct_gibbs_fit(fit_curr, sd_q)
fit_ma_wrap   <- mct_gibbs_fit(fit,      sd_q)

trend_curr <- trend_path(fit_curr_wrap, sd_q$dates)
trend_ma   <- trend_path(fit_ma_wrap,   sd_q$dates)

cat("\n=== No-MA (current production w/ outliers) vs MA(3) ===\n")
cat(sprintf("Latest-t trend:\n"))
cat(sprintf("  No-MA: median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(trend_curr$median, 1L),
            tail(trend_curr$lower, 1L), tail(trend_curr$upper, 1L)))
cat(sprintf("  MA(3): median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(trend_ma$median, 1L),
            tail(trend_ma$lower, 1L), tail(trend_ma$upper, 1L)))

covid_dates <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",
                         "2020-12-01", "2021-03-01"))
cat("\nTrend around COVID:\n")
cat(sprintf("%-12s  %-22s  %-22s  diff\n", "date",
            "No-MA (median 68%CI)", "MA(3) (median 68%CI)"))
for (d in covid_dates) {
  i_c <- which.min(abs(as.numeric(trend_curr$date - d)))
  i_m <- which.min(abs(as.numeric(trend_ma$date   - d)))
  cat(sprintf(
    "%s   %6.3f [%6.3f,%6.3f]   %6.3f [%6.3f,%6.3f]   %+.3f\n",
    format(trend_curr$date[i_c], "%Y-%m"),
    trend_curr$median[i_c], trend_curr$lower[i_c], trend_curr$upper[i_c],
    trend_ma$median[i_m],   trend_ma$lower[i_m],   trend_ma$upper[i_m],
    trend_ma$median[i_m] - trend_curr$median[i_c]
  ))
}

# Posterior theta summary
theta_med <- apply(fit$draws$theta, c(1, 2), median)
cat("\nPosterior median theta per sub-group (quarterly MA(3)):\n")
cat(sprintf("%-50s  %6s  %6s  %6s\n", "sub-group", "theta1", "theta2", "theta3"))
for (i in seq_len(sd_q$N)) {
  cat(sprintf("  %-50s  %+6.3f %+6.3f %+6.3f\n",
              sd_q$groups[i],
              theta_med[i, 1L], theta_med[i, 2L], theta_med[i, 3L]))
}

cat("\n[step24] done.\n")
