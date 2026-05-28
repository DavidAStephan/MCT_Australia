# Step 25 — MA(1) with tighter theta prior on the quarterly sub-group
# panel. MA(3) at the default prior blew up (theta values ranging
# +33 to -155, structurally non-invertible). Reducing q and tightening
# the prior should constrain identification and avoid the blow-up.
#
# Setup:
#   - q_MA = 1 (state aug 1 + N + 2*N = 67 vs MA(3)'s 166)
#   - theta_prec_prior = 10 (prior SD on theta ≈ 0.32; the default 0.1
#     allowed SD ≈ 3.2 which is way too loose).
#   - Variant B, quarterly cadence, no outliers (test MA alone first).
# Projected fit time: ~15-20 min.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

sd_q <- readRDS("outputs/draws/stan_data.rds")
y    <- sd_q$y

cat(sprintf("[step25] sd_q: T=%d N=%d  ref=%s\n",
            sd_q$T, sd_q$N, sd_q$groups[sd_q$ref]))

# Smoke test
set.seed(99L)
t0 <- Sys.time()
fit_smoke <- fit_mct_gibbs(
  y, q_MA = 1L, ref = sd_q$ref,
  n_burn = 3L, n_draw = 3L, verbose = FALSE,
  config = list(variant            = "B",
                s2_gam_prior       = 1 / (60 * 4),
                theta_prec_prior   = 10)
)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[step25] Smoke %.1fs (6 iters); projected 3000-iter = %.0f min\n",
            secs, secs / 6 * 3000 / 60))

# Full fit
out_path <- "outputs/draws/fit_B_ma1.rds"
if (file.exists(out_path)) {
  cat("\n[step25] cached at ", out_path, "\n", sep = "")
  fit <- readRDS(out_path)
} else {
  cat("\n[step25] Fitting Variant B + MA(1) ...\n")
  set.seed(20260528L)
  t0 <- Sys.time()
  fit <- fit_mct_gibbs(
    y, q_MA = 1L, ref = sd_q$ref,
    n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
    config = list(variant          = "B",
                  s2_gam_prior     = 1 / (60 * 4),
                  theta_prec_prior = 10)
  )
  cat(sprintf("[step25] done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit, out_path)
}

# Theta diagnostic — should be in (-1, 1) for invertibility
theta_med <- apply(fit$draws$theta, c(1, 2), median)
cat("\nPosterior median theta_1 per sub-group (should be in (-1, 1)):\n")
in_range <- 0L; out_range <- 0L
for (i in seq_len(sd_q$N)) {
  ok <- abs(theta_med[i, 1L]) < 1
  if (ok) in_range <- in_range + 1L else out_range <- out_range + 1L
  cat(sprintf("  %-50s theta_1=%+6.3f  %s\n",
              sd_q$groups[i], theta_med[i, 1L],
              if (ok) "" else "  ** out of (-1,1) **"))
}
cat(sprintf("\nIn range (-1,1): %d / %d sub-groups\n",
            in_range, sd_q$N))

# Headline comparison vs outliers-only production
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)
fit_prod <- readRDS("outputs/draws/fit_B.rds")
fit_prod_wrap <- mct_gibbs_fit(fit_prod, sd_q)
fit_ma_wrap   <- mct_gibbs_fit(fit,     sd_q)

trend_prod <- trend_path(fit_prod_wrap, sd_q$dates)
trend_ma   <- trend_path(fit_ma_wrap,   sd_q$dates)
cs_prod    <- common_share_path(fit_prod_wrap, sd_q$dates)
cs_ma      <- common_share_path(fit_ma_wrap,   sd_q$dates)

cat("\n=== Outliers-only (production) vs MA(1) ===\n")
cat(sprintf("Latest-t trend:\n"))
cat(sprintf("  Outliers: median=%.3f  68%%=[%.3f, %.3f]  width=%.2fpp\n",
            tail(trend_prod$median, 1L),
            tail(trend_prod$lower, 1L), tail(trend_prod$upper, 1L),
            tail(trend_prod$upper, 1L) - tail(trend_prod$lower, 1L)))
cat(sprintf("  MA(1):    median=%.3f  68%%=[%.3f, %.3f]  width=%.2fpp\n",
            tail(trend_ma$median, 1L),
            tail(trend_ma$lower, 1L), tail(trend_ma$upper, 1L),
            tail(trend_ma$upper, 1L) - tail(trend_ma$lower, 1L)))
cat(sprintf("\nLatest-t common-share:\n"))
cat(sprintf("  Outliers: %.3f\n", tail(cs_prod$median, 1L)))
cat(sprintf("  MA(1):    %.3f\n", tail(cs_ma$median, 1L)))

covid_dates <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",
                         "2020-12-01"))
cat("\nTrend around COVID:\n")
cat(sprintf("%-12s  %-22s  %-22s  diff\n", "date",
            "Outliers (median 68%CI)", "MA(1)   (median 68%CI)"))
for (d in covid_dates) {
  i_p <- which.min(abs(as.numeric(trend_prod$date - d)))
  i_m <- which.min(abs(as.numeric(trend_ma$date   - d)))
  cat(sprintf(
    "%s   %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]  %+.3f\n",
    format(trend_prod$date[i_p], "%Y-%m"),
    trend_prod$median[i_p], trend_prod$lower[i_p], trend_prod$upper[i_p],
    trend_ma$median[i_m],   trend_ma$lower[i_m],   trend_ma$upper[i_m],
    trend_ma$median[i_m] - trend_prod$median[i_p]
  ))
}

cat("\n[step25] done.\n")
