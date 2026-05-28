# Step 23 — Re-fit the production model (Variant B, quarterly sub-groups)
# with the outlier scale mixture turned ON. Compare to the no-outliers
# fit to see whether the COVID-period trend dip is dampened.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("R",       "\\.R$", full.names = TRUE)) source(f)

sd_q <- readRDS("outputs/draws/stan_data.rds")
y    <- sd_q$y

cat(sprintf(
  "[step23] sd_q: T=%d N=%d  n_obs_q=%d  ref=%s\n",
  sd_q$T, sd_q$N, sd_q$n_obs_q, sd_q$groups[sd_q$ref]
))

out_path <- "outputs/draws/fit_B_outliers.rds"
if (file.exists(out_path)) {
  cat("\n[step23] cached at ", out_path, " — loading.\n", sep = "")
  fit_out <- readRDS(out_path)
} else {
  cat("\n[step23] Fitting Variant B with use_outliers = TRUE ...\n")
  set.seed(20260528L)
  t0 <- Sys.time()
  fit_out <- fit_mct_gibbs(
    y, ref = sd_q$ref,
    n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
    config = list(variant = "B",
                  s2_gam_prior = 1 / (60 * 4),
                  use_outliers = TRUE)
  )
  cat(sprintf("[step23] done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit_out, out_path)
  cat("[step23] saved ", out_path, "\n", sep = "")
}

# Headline comparison: outliers ON vs OFF
fit_base <- readRDS("outputs/draws/fit_B.rds")
fit_base_wrap <- mct_gibbs_fit(fit_base, sd_q)
fit_out_wrap  <- mct_gibbs_fit(fit_out,  sd_q)

trend_base <- trend_path(fit_base_wrap, sd_q$dates)
trend_out  <- trend_path(fit_out_wrap,  sd_q$dates)
cs_base    <- common_share_path(fit_base_wrap, sd_q$dates)
cs_out     <- common_share_path(fit_out_wrap,  sd_q$dates)

cat("\n=== Outliers OFF vs ON: latest-t trend ===\n")
cat(sprintf("  OFF: median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(trend_base$median, 1L),
            tail(trend_base$lower, 1L), tail(trend_base$upper, 1L)))
cat(sprintf("  ON:  median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(trend_out$median, 1L),
            tail(trend_out$lower, 1L), tail(trend_out$upper, 1L)))

cat("\n=== Outliers OFF vs ON: common-share ===\n")
cat(sprintf("  OFF: latest median=%.3f\n", tail(cs_base$median, 1L)))
cat(sprintf("  ON:  latest median=%.3f\n", tail(cs_out$median, 1L)))

# Trend at COVID period: 2020-Q2 (the smoking gun)
covid_dates <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",
                         "2020-12-01", "2021-03-01"))
cat("\n=== Trend trajectory around COVID (the test) ===\n")
cat(sprintf("%-12s  %-22s  %-22s  diff\n", "date", "OFF (median 68%CI)",
            "ON  (median 68%CI)"))
for (d in covid_dates) {
  i_b <- which.min(abs(as.numeric(trend_base$date - d)))
  i_o <- which.min(abs(as.numeric(trend_out$date  - d)))
  cat(sprintf(
    "%s   %6.3f [%6.3f,%6.3f]   %6.3f [%6.3f,%6.3f]   %+.3f\n",
    format(trend_base$date[i_b], "%Y-%m"),
    trend_base$median[i_b], trend_base$lower[i_b], trend_base$upper[i_b],
    trend_out$median[i_o],  trend_out$lower[i_o],  trend_out$upper[i_o],
    trend_out$median[i_o] - trend_base$median[i_b]
  ))
}

# Outlier flagging: which (date, sub-group) cells were flagged?
s_med <- apply(fit_out$draws$s_outlier, c(1, 2), median)
n_obs <- sum(!is.na(sd_q$y))
n_flagged <- sum(s_med > 1.5 & !is.na(sd_q$y))
cat(sprintf("\nObservations flagged as outliers (s_med > 1.5): %d / %d (%.1f%%)\n",
            n_flagged, n_obs, 100 * n_flagged / n_obs))

cat("\nMean ps per sub-group (probability of normal obs):\n")
ps_mean <- rowMeans(fit_out$draws$ps)
for (i in seq_len(sd_q$N)) {
  cat(sprintf("  %-50s  ps=%.3f\n", sd_q$groups[i], ps_mean[i]))
}

# Top-5 flagged (date, sub-group) cells
cat("\nTop-10 flagged (date, sub-group) by posterior median s_outlier:\n")
flat_idx <- order(-s_med * !is.na(sd_q$y))[1:10]
for (k in flat_idx) {
  t_i <- ((k - 1L) %% sd_q$T) + 1L
  i_i <- ((k - 1L) %/% sd_q$T) + 1L
  cat(sprintf("  %s  %-45s  s_med=%.2f\n",
              format(sd_q$dates[t_i], "%Y-Q%q"),
              sd_q$groups[i_i], s_med[t_i, i_i]))
}

saveRDS(list(trend_base = trend_base, trend_out = trend_out,
             cs_base = cs_base, cs_out = cs_out, s_med = s_med),
        "outputs/draws/outliers_comparison.rds")
cat("\n[step23] Saved comparison to outputs/draws/outliers_comparison.rds\n")
