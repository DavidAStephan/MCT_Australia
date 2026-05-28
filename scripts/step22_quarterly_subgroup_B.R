# Step 22 — Stock-Watson MCT (Variant B, RW common trend) on the
# QUARTERLY-clocked sub-group panel (N=33, T=133). This is the
# "match NY Fed in spirit, use AU's quarterly data honestly" model.
#
# Differences from the failed Step 21:
#   - Latent clock is quarterly, not monthly. Avoids the 3:1
#     states-to-obs problem that blew up sigma_s.
#   - Uses the non-mixed-frequency Gibbs path (`gibbs_sweep` +
#     `build_mct_ssm`), now Variant-B-aware.
#   - SV prior s2_gam_prior calibrated for quarterly (1/(60*4))
#     not monthly (1/(60*12)).

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("R",       "\\.R$", full.names = TRUE)) source(f)

# Build a fresh quarterly stan_data.
cpi_sub  <- fetch_cpi_quarterly_subgroup()
w_sub    <- fetch_weights_subgroup()
infl_sub <- build_inflation_series(cpi_sub) |>
  demean_series(window = c(1993, 2019))

sd_q <- build_stan_data_quarterly(
  infl_sub, weights = w_sub,
  ref = "Rents", start = as.Date("1993-01-01")
)
saveRDS(sd_q, "outputs/draws/stan_data_subgroup_q.rds")

cat(sprintf(
  "[step22] Quarterly sub-group stan_data: T=%d quarters  N=%d  n_obs_q=%d\n",
  sd_q$T, sd_q$N, sd_q$n_obs_q
))
cat("Date range:", as.character(min(sd_q$dates)),
    "->", as.character(max(sd_q$dates)), "\n")
cat("Reference sub-group (ref=", sd_q$ref, "): ", sd_q$groups[sd_q$ref],
    "\n", sep = "")

# Use the dense y matrix directly (sd_q$y has NA for missing cells —
# the non-mixed Gibbs path handles those).
y <- sd_q$y

# Smoke test first
t0 <- Sys.time()
fit_smoke <- fit_mct_gibbs(
  y, ref = sd_q$ref,
  n_burn = 5L, n_draw = 5L, verbose = FALSE,
  config = list(variant = "B",
                s2_gam_prior = 1 / (60 * 4))   # quarterly cadence
)
sec_per_iter <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 10
cat(sprintf("\n[step22] Smoke per-iter: %.2f sec; projected 3000-iter = %.0f min\n",
            sec_per_iter, sec_per_iter * 3000 / 60))

# Full fit
out_path <- "outputs/draws/fit_gibbs_B_subgroup_q.rds"
if (file.exists(out_path)) {
  cat("\n[step22] cached at ", out_path, " — loading.\n", sep = "")
  fit <- readRDS(out_path)
} else {
  cat("\n[step22] Fitting Variant B (quarterly clock, sub-groups) ...\n")
  set.seed(1L)
  t0 <- Sys.time()
  fit <- fit_mct_gibbs(
    y, ref = sd_q$ref,
    n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
    config = list(variant = "B",
                  s2_gam_prior = 1 / (60 * 4))
  )
  cat(sprintf("[step22] done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit, out_path)
  cat("[step22] saved ", out_path, "\n", sep = "")
}

# Headline summary
fit_wrap <- mct_gibbs_fit(fit, sd_q)
trend    <- trend_path(fit_wrap, sd_q$dates)
cs       <- common_share_path(fit_wrap, sd_q$dates)

cat("\n=== Quarterly sub-group Variant B headline ===\n")
cat(sprintf("Median posterior sigma_c (RW common-trend SD per quarter): %.3f\n",
            median(fit$draws$sigma_c)))
cat(sprintf("Median posterior sigma_s (pooled): %.3f\n",
            median(fit$draws$sigma_s)))
cat(sprintf("Median posterior sigma_eps (pooled): %.3f\n",
            median(fit$draws$sigma_eps)))

cat(sprintf("\nLatest-t (%s) trend: median=%.3f  68%%=[%.3f, %.3f]  width=%.2fpp\n",
            format(tail(sd_q$dates, 1L), "%Y-%m"),
            tail(trend$median, 1L),
            tail(trend$lower, 1L), tail(trend$upper, 1L),
            tail(trend$upper, 1L) - tail(trend$lower, 1L)))
cat(sprintf("Latest-t common-share: median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(cs$median, 1L),
            tail(cs$lower, 1L), tail(cs$upper, 1L)))

# Comparison to production 11-group (Stan Variant A on monthly mixed-freq)
cat("\n=== Comparison to current production (11-group mixed-freq Variant A) ===\n")
tar_load_silent <- function(x) suppressMessages(targets::tar_read(x))
trend_prod <- tar_load_silent(trend_A)
cs_prod    <- tar_load_silent(common_share_A)
cat(sprintf("Production 11-group A trend latest-t: median=%.3f  68%%=[%.3f, %.3f]  width=%.2fpp\n",
            tail(trend_prod$median, 1L),
            tail(trend_prod$lower, 1L), tail(trend_prod$upper, 1L),
            tail(trend_prod$upper, 1L) - tail(trend_prod$lower, 1L)))
cat(sprintf("Production 11-group A common-share latest-t: median=%.3f\n",
            tail(cs_prod$median, 1L)))
