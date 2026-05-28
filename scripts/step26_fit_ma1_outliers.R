# Step 26 — MA(1) + outliers (stacked) on the quarterly sub-group panel.
# Tests whether the two complementary mechanisms (MA bleeds out shocks
# over time, outliers absorb anomalous-volatility periods) improve over
# MA(1) alone (step25) and outliers alone (current production).

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

sd_q <- readRDS("outputs/draws/stan_data.rds")
y    <- sd_q$y

out_path <- "outputs/draws/fit_B_ma1_outliers.rds"
if (file.exists(out_path)) {
  cat("[step26] cached at ", out_path, "\n", sep = "")
  fit <- readRDS(out_path)
} else {
  cat("[step26] Fitting Variant B + MA(1) + outliers ...\n")
  set.seed(20260528L)
  t0 <- Sys.time()
  fit <- fit_mct_gibbs(
    y, q_MA = 1L, ref = sd_q$ref,
    n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
    config = list(variant          = "B",
                  s2_gam_prior     = 1 / (60 * 4),
                  theta_prec_prior = 10,
                  use_outliers     = TRUE)
  )
  cat(sprintf("[step26] done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit, out_path)
}

# Diagnostics
theta_med <- apply(fit$draws$theta, c(1, 2), median)
in_range <- sum(abs(theta_med[, 1L]) < 1)
cat(sprintf("\nThetas in (-1,1): %d / %d\n", in_range, sd_q$N))

s_med <- apply(fit$draws$s_outlier, c(1, 2), median)
n_obs <- sum(!is.na(sd_q$y))
n_flag <- sum(s_med > 1.5 & !is.na(sd_q$y))
cat(sprintf("Obs flagged (s_med > 1.5): %d / %d (%.1f%%)\n",
            n_flag, n_obs, 100 * n_flag / n_obs))

# Three-way comparison
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)
fit_out_only <- readRDS("outputs/draws/fit_B.rds")            # outliers, no MA
fit_ma_only  <- readRDS("outputs/draws/fit_B_ma1.rds")         # MA(1), no outliers

wraps <- list(
  out_only = mct_gibbs_fit(fit_out_only, sd_q),
  ma_only  = mct_gibbs_fit(fit_ma_only,  sd_q),
  both     = mct_gibbs_fit(fit,          sd_q)
)
trends <- lapply(wraps, function(w) trend_path(w, sd_q$dates))
css    <- lapply(wraps, function(w) common_share_path(w, sd_q$dates))

cat("\n=== Three-way comparison (latest-t) ===\n")
cat(sprintf("%-12s  %-22s  %-22s  %-8s\n",
            "config", "trend median 68%CI", "common-share 68%CI", "width(pp)"))
for (nm in names(wraps)) {
  tr <- trends[[nm]]; cs <- css[[nm]]
  cat(sprintf("%-12s  %6.3f [%6.3f,%6.3f]  %5.3f [%5.3f,%5.3f]  %.2f\n",
              nm,
              tail(tr$median, 1L), tail(tr$lower, 1L), tail(tr$upper, 1L),
              tail(cs$median, 1L), tail(cs$lower, 1L), tail(cs$upper, 1L),
              tail(tr$upper, 1L) - tail(tr$lower, 1L)))
}

covid_dates <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",
                         "2020-12-01"))
cat("\n=== Trend around COVID ===\n")
cat(sprintf("%-12s  %-18s  %-18s  %-18s\n", "date",
            "outliers-only", "MA(1)-only", "MA(1)+outliers"))
for (d in covid_dates) {
  row <- format(d, "%Y-%m")
  for (nm in names(wraps)) {
    tr <- trends[[nm]]
    i <- which.min(abs(as.numeric(tr$date - d)))
    row <- paste0(row, sprintf("  %6.3f [%5.2f,%5.2f]",
                               tr$median[i], tr$lower[i], tr$upper[i]))
  }
  cat(row, "\n")
}
