# Step 17 — Fit Variant A on the 1993-onward sample (post-RBA-inflation-
# targeting cut) with the 1993-2019 demeaning window.
#
# Doesn't touch the production _targets pipeline. Saves the new
# stan_data and the new Variant A fit alongside the existing 1990-based
# files so we can compare directly.

suppressMessages({
  library(targets)
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(lubridate)
})
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

cpi_q_clean <- tar_read(cpi_q_clean)
cpi_m_clean <- tar_read(cpi_m_clean)
weights     <- tar_read(weights)

combined <- dplyr::bind_rows(cpi_q_clean, cpi_m_clean)

# New demeaning: 1993-2019 (full pre-pandemic in-sample period).
infl_demeaned <- demean_series(combined, window = c(1993, 2019))

stan_data_1993 <- build_stan_data(
  infl_demeaned, weights = weights,
  start = as.Date("1993-01-01")
)
saveRDS(stan_data_1993, "outputs/draws/stan_data_1993.rds")
cat(sprintf(
  "[step17] stan_data_1993: T=%d  N=%d  n_obs_m=%d  n_obs_q=%d  dates: %s -> %s\n",
  stan_data_1993$T, stan_data_1993$N,
  stan_data_1993$n_obs_m, stan_data_1993$n_obs_q,
  min(stan_data_1993$dates), max(stan_data_1993$dates)
))

# Stash the original for comparison
saveRDS(combined, "outputs/draws/cpi_combined_cached.rds")

cat("[step17] Fitting Variant A (~90 min) ...\n")
t0 <- Sys.time()

# Save into a separate file name so the 1990-onward fit_A.rds is untouched.
# We sidestep fit_mct's hard-coded save_dir/filename by saving manually.
fit_A_1993 <- fit_mct(
  stan_data_1993,
  variant         = "A",
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1500,
  iter_sampling   = 1500,
  adapt_delta     = 0.95,
  max_treedepth   = 12,
  seed            = 19930101,
  save_dir        = NULL   # don't auto-save with the canonical name
)
fit_A_1993$save_object("outputs/draws/fit_A_1993.rds")

elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\n[step17] Variant A (1993) fit complete in %.1f min\n",
            elapsed_min))
cat("[step17] diagnostics:\n")
print(fit_A_1993$diagnostic_summary())

# Headline comparison vs the original 1990-onward fit
fit_A_1990 <- readRDS("outputs/draws/fit_A.rds")
stan_data_1990 <- readRDS("outputs/draws/stan_data.rds")

rho_1990 <- as.numeric(fit_A_1990$draws("rho", format = "matrix"))
rho_1993 <- as.numeric(fit_A_1993$draws("rho", format = "matrix"))

cat("\n[step17] rho posterior:\n")
cat(sprintf("  1990 sample:  median=%.3f  68%%=[%.3f, %.3f]\n",
            median(rho_1990),
            quantile(rho_1990, 0.16), quantile(rho_1990, 0.84)))
cat(sprintf("  1993 sample:  median=%.3f  68%%=[%.3f, %.3f]\n",
            median(rho_1993),
            quantile(rho_1993, 0.16), quantile(rho_1993, 0.84)))

# Latest-t trend + common-share (overlap: both samples end at same date)
trend_1990 <- trend_path(fit_A_1990, stan_data_1990$dates)
trend_1993 <- trend_path(fit_A_1993, stan_data_1993$dates)
cs_1990    <- common_share_path(fit_A_1990, stan_data_1990$dates)
cs_1993    <- common_share_path(fit_A_1993, stan_data_1993$dates)

cat(sprintf(
  "\n[step17] latest-t trend median:  1990=%.3f  1993=%.3f  (diff %+.3f)\n",
  tail(trend_1990$median, 1L), tail(trend_1993$median, 1L),
  tail(trend_1993$median, 1L) - tail(trend_1990$median, 1L)))
cat(sprintf(
  "[step17] latest-t trend 68%%:     1990=[%.3f, %.3f]  1993=[%.3f, %.3f]\n",
  tail(trend_1990$lower, 1L), tail(trend_1990$upper, 1L),
  tail(trend_1993$lower, 1L), tail(trend_1993$upper, 1L)))
cat(sprintf(
  "[step17] latest-t common-share:   1990=%.3f  1993=%.3f\n",
  tail(cs_1990$median, 1L), tail(cs_1993$median, 1L)))

cat("\n[step17] Saved: outputs/draws/fit_A_1993.rds\n")
cat("[step17] done.\n")
