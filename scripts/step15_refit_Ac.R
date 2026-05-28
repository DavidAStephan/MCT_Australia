# Step 15 — Refit Stan Variant Ac (constant-lambda AR(1)) and save the
# fit object so it can be used for the Gibbs-vs-Stan apples-to-apples
# comparison.
#
# Background: Per next_session.md, the Step 11 trend-level gap between
# Stan Variant A (TVP λ) and Gibbs (constant λ) is a model-spec
# difference, not a sampler bug. The fair Stan comparison is against
# Variant Ac (constant λ), but the previous Ac fit was deleted in the
# May speedup exploration. This script regenerates it.
#
# Difference from benchmark_Ac.R: save the fit object via save_dir.

suppressMessages({
  library(cmdstanr)
  library(posterior)
})
for (f in list.files("R",       "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

stan_data <- readRDS("outputs/draws/stan_data.rds")
cat(sprintf(
  "[step15] stan_data: T=%d N=%d  n_obs_m=%d n_obs_q=%d\n",
  stan_data$T, stan_data$N, stan_data$n_obs_m, stan_data$n_obs_q
))

t0 <- Sys.time()
fit_Ac <- fit_mct(
  stan_data,
  variant         = "Ac",
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1500,
  iter_sampling   = 1500,
  adapt_delta     = 0.97,
  max_treedepth   = 12,
  seed            = 20260526,
  save_dir        = "outputs/draws"   # writes outputs/draws/fit_Ac.rds
)
elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

cat(sprintf("\n[step15] Ac fit complete in %.1f min\n", elapsed_min))
cat("[step15] diagnostics:\n")
print(fit_Ac$diagnostic_summary())

# Gibbs apples-to-apples comparison: trend + rho posterior
gibbs_path <- "outputs/draws/fit_gibbs_baseline.rds"
if (file.exists(gibbs_path)) {
  fit_gib_raw <- readRDS(gibbs_path)
  rho_Ac  <- as.numeric(fit_Ac$draws("rho", format = "matrix"))
  rho_gib <- fit_gib_raw$draws$rho
  cat("\n[step15] rho posterior (Stan Ac vs Gibbs baseline):\n")
  cat(sprintf("  Stan Ac:   median=%.3f  68%%=[%.3f, %.3f]\n",
              median(rho_Ac),
              quantile(rho_Ac,  0.16), quantile(rho_Ac,  0.84)))
  cat(sprintf("  Gibbs:     median=%.3f  68%%=[%.3f, %.3f]\n",
              median(rho_gib),
              quantile(rho_gib, 0.16), quantile(rho_gib, 0.84)))

  # Latest-t trend
  trend_Ac <- trend_path(fit_Ac, stan_data$dates)
  fit_gib_wrap <- mct_gibbs_fit(fit_gib_raw, stan_data)
  trend_gib    <- trend_path(fit_gib_wrap, stan_data$dates)
  cat(sprintf(
    "\n[step15] latest-t trend: Stan Ac = %.3f  |  Gibbs = %.3f  (diff %+.3f)\n",
    tail(trend_Ac$median, 1L), tail(trend_gib$median, 1L),
    tail(trend_gib$median, 1L) - tail(trend_Ac$median, 1L)
  ))
} else {
  cat("\n[step15] (Gibbs baseline fit not found; skipped comparison.)\n")
}

cat("\n[step15] Saved: outputs/draws/fit_Ac.rds\n")
