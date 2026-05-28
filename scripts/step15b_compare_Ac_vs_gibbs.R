# Step 15b — Post-fit apples-to-apples comparison: Stan Ac vs Gibbs baseline.
# Both are constant-lambda AR(1) common, no outliers — the closest
# possible comparison. Re-runnable without re-fitting.

suppressMessages({
  library(cmdstanr)
  library(posterior)
})
for (f in list.files("R",       "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

stan_data   <- readRDS("outputs/draws/stan_data.rds")
fit_Ac      <- readRDS("outputs/draws/fit_Ac.rds")
fit_gib_raw <- readRDS("outputs/draws/fit_gibbs_baseline.rds")

rho_Ac  <- as.numeric(fit_Ac$draws("rho", format = "matrix"))
rho_gib <- fit_gib_raw$draws$rho
cat("\n[step15b] rho posterior (Stan Ac vs Gibbs baseline):\n")
cat(sprintf("  Stan Ac:   median=%.3f  68%%=[%.3f, %.3f]\n",
            median(rho_Ac),
            quantile(rho_Ac,  0.16), quantile(rho_Ac,  0.84)))
cat(sprintf("  Gibbs:     median=%.3f  68%%=[%.3f, %.3f]\n",
            median(rho_gib),
            quantile(rho_gib, 0.16), quantile(rho_gib, 0.84)))

trend_Ac     <- trend_path(fit_Ac, stan_data$dates)
fit_gib_wrap <- mct_gibbs_fit(fit_gib_raw, stan_data)
trend_gib    <- trend_path(fit_gib_wrap, stan_data$dates)

cs_Ac  <- common_share_path(fit_Ac,    stan_data$dates)
cs_gib <- common_share_path(fit_gib_wrap, stan_data$dates)

cat(sprintf(
  "\n[step15b] latest-t trend (median):  Stan Ac=%.3f  Gibbs=%.3f  (diff %+.3f)\n",
  tail(trend_Ac$median, 1L), tail(trend_gib$median, 1L),
  tail(trend_gib$median, 1L) - tail(trend_Ac$median, 1L)
))
cat(sprintf(
  "[step15b] latest-t trend 68%%:  Stan Ac=[%.3f, %.3f]  Gibbs=[%.3f, %.3f]\n",
  tail(trend_Ac$lower, 1L), tail(trend_Ac$upper, 1L),
  tail(trend_gib$lower, 1L), tail(trend_gib$upper, 1L)))
cat(sprintf(
  "[step15b] latest-t common-share (median):  Stan Ac=%.3f  Gibbs=%.3f\n",
  tail(cs_Ac$median, 1L), tail(cs_gib$median, 1L)))

tr_diff <- trend_gib$median - trend_Ac$median
cat(sprintf(
  "\n[step15b] RMS difference of median trends (all t=1..%d):  %.3f\n",
  nrow(trend_Ac), sqrt(mean(tr_diff^2))))
cat(sprintf("[step15b] max |diff| of medians:  %.3f at t=%d\n",
            max(abs(tr_diff)), which.max(abs(tr_diff))))
