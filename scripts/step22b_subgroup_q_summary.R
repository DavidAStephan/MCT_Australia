# Step 22b — Re-build the quarterly sub-group stan_data (with the
# updated builder that surfaces obs in *_m slots so the log_lik
# computation works) and produce the headline summary against the
# already-fit fit_gibbs_B_subgroup_q.rds.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("R",       "\\.R$", full.names = TRUE)) source(f)

cpi_sub  <- fetch_cpi_quarterly_subgroup()
w_sub    <- fetch_weights_subgroup()
infl_sub <- build_inflation_series(cpi_sub) |>
  demean_series(window = c(1993, 2019))

sd_q <- build_stan_data_quarterly(
  infl_sub, weights = w_sub,
  ref = "Rents", start = as.Date("1993-01-01")
)
saveRDS(sd_q, "outputs/draws/stan_data_subgroup_q.rds")

fit <- readRDS("outputs/draws/fit_gibbs_B_subgroup_q.rds")
fit_wrap <- mct_gibbs_fit(fit, sd_q)
trend    <- trend_path(fit_wrap, sd_q$dates)
cs       <- common_share_path(fit_wrap, sd_q$dates)

cat("=== Quarterly sub-group Variant B headline ===\n")
cat(sprintf("T=%d quarters  N=%d  n_obs=%d  reference=%s\n",
            sd_q$T, sd_q$N, sd_q$n_obs_m, sd_q$groups[sd_q$ref]))
cat(sprintf("\nMedian posterior SVs:\n"))
cat(sprintf("  sigma_c (RW common-trend SD per quarter): %.3f\n",
            median(fit$draws$sigma_c)))
cat(sprintf("  sigma_s (pooled):                          %.3f\n",
            median(fit$draws$sigma_s)))
cat(sprintf("  sigma_eps (pooled):                        %.3f\n",
            median(fit$draws$sigma_eps)))

cat(sprintf("\nLatest-t (%s) trend: median=%.3f  68%%=[%.3f, %.3f]  width=%.2fpp\n",
            format(tail(sd_q$dates, 1L), "%Y-%m"),
            tail(trend$median, 1L),
            tail(trend$lower, 1L), tail(trend$upper, 1L),
            tail(trend$upper, 1L) - tail(trend$lower, 1L)))
cat(sprintf("Latest-t common-share: median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(cs$median, 1L), tail(cs$lower, 1L), tail(cs$upper, 1L)))

# Compare to production 11-group Variant A (monthly mixed-freq Stan)
trend_prod <- suppressMessages(targets::tar_read(trend_A))
cs_prod    <- suppressMessages(targets::tar_read(common_share_A))
cat("\n=== Comparison to current production (11-group monthly mixed-freq Variant A) ===\n")
cat(sprintf("Production trend latest-t: median=%.3f  68%%=[%.3f, %.3f]  width=%.2fpp\n",
            tail(trend_prod$median, 1L),
            tail(trend_prod$lower, 1L), tail(trend_prod$upper, 1L),
            tail(trend_prod$upper, 1L) - tail(trend_prod$lower, 1L)))
cat(sprintf("Production common-share latest-t: median=%.3f\n",
            tail(cs_prod$median, 1L)))

# Show the trend trajectory at several key dates
cat("\nQuarterly sub-group B trend trajectory:\n")
keys <- as.Date(c("1995-03-01", "2000-03-01", "2005-03-01", "2010-03-01",
                  "2015-03-01", "2020-03-01", "2024-03-01", "2026-03-01"))
for (d in keys) {
  i <- which.min(abs(as.numeric(trend$date - d)))
  cat(sprintf("  %s  trend=%.3f  68%%=[%.3f, %.3f]\n",
              format(trend$date[i], "%Y-%m"),
              trend$median[i], trend$lower[i], trend$upper[i]))
}

# Save processed paths for the dashboard wiring later
saveRDS(list(stan_data = sd_q, fit = fit, fit_wrap = fit_wrap,
             trend = trend, common_share = cs),
        "outputs/draws/subgroup_q_B_summary.rds")
cat("\n[step22b] Saved outputs/draws/subgroup_q_B_summary.rds\n")
