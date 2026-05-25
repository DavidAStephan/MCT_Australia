# Benchmark: Variant Ac (constant-lambda A) vs the cached Variant A baseline.
# Same production config as the original fit so the wall-time delta is a
# clean read on Phase 1 #2 (constant-lambda). Phase 1 #5 (tighter adapt
# buffers) is also enabled via init_buffer/term_buffer/window.
#
# Run with: Rscript scripts/benchmark_Ac.R 2>&1 | tee outputs/benchmark_Ac.log
#
# Compare-to-baseline: README states fit_A took ~93 min on the M1.

suppressMessages({
  library(cmdstanr)
  library(posterior)
})
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

stan_data <- readRDS("outputs/draws/stan_data.rds")
cat(sprintf(
  "[benchmark] stan_data: T=%d N=%d  n_obs_m=%d n_obs_q=%d\n",
  stan_data$T, stan_data$N, stan_data$n_obs_m, stan_data$n_obs_q
))

t0 <- Sys.time()
fit_Ac <- fit_mct(
  stan_data,
  variant = "Ac",
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 1500,
  adapt_delta = 0.97,
  max_treedepth = 12,
  # Stan defaults for adapt buffers — isolate Ac-vs-A delta cleanly.
  # The wider-buffer config (#5) will be measured separately in a
  # follow-up run after we have a clean Ac baseline.
  seed = 20260524,
  # Benchmark only — don't write a ~2 GB rds. Comparison happens in-memory.
  # If Ac validates, future runs will save as the canonical fit_A.rds and
  # the original Variant A model can be retired.
  save_dir = NULL
)
elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

cat(sprintf(
  "\n[benchmark] Ac fit complete in %.1f min (baseline A: ~93 min)\n",
  elapsed_min
))

cat("\n[benchmark] Per-chain diagnostics:\n")
print(fit_Ac$diagnostic_summary())

cat("\n[benchmark] Stan profile output (where the time went):\n")
prof <- tryCatch(fit_Ac$profiles(), error = function(e) NULL)
if (!is.null(prof) && length(prof) > 0) {
  # Aggregate across chains: sum total_time per profile name.
  prof_df <- do.call(rbind, prof)
  agg <- aggregate(prof_df$total_time, by = list(name = prof_df$name), FUN = sum)
  agg <- agg[order(-agg$x), ]
  for (i in seq_len(nrow(agg))) {
    cat(sprintf("  %-22s  %7.1f s  (sum across chains)\n",
                agg$name[i], agg$x[i]))
  }
} else {
  cat("  (no profile output captured)\n")
}

cat("\n[benchmark] Quick consistency check vs baseline fit_A:\n")
baseline_path <- "outputs/draws/fit_A.rds"
if (file.exists(baseline_path)) {
  fit_A <- readRDS(baseline_path)
  trend_A  <- trend_path(fit_A, stan_data$dates)
  trend_Ac <- trend_path(fit_Ac, stan_data$dates)
  cs_A     <- common_share_path(fit_A, stan_data$dates)
  cs_Ac    <- common_share_path(fit_Ac, stan_data$dates)
  cat(sprintf(
    "  latest trend     median:   A=%6.3f   Ac=%6.3f   diff=%+.3f\n",
    tail(trend_A$median, 1), tail(trend_Ac$median, 1),
    tail(trend_Ac$median, 1) - tail(trend_A$median, 1)
  ))
  cat(sprintf(
    "  latest trend     lower:    A=%6.3f   Ac=%6.3f\n",
    tail(trend_A$lower, 1), tail(trend_Ac$lower, 1)
  ))
  cat(sprintf(
    "  latest trend     upper:    A=%6.3f   Ac=%6.3f\n",
    tail(trend_A$upper, 1), tail(trend_Ac$upper, 1)
  ))
  cat(sprintf(
    "  latest common-share med:   A=%6.3f   Ac=%6.3f   diff=%+.3f\n",
    tail(cs_A$median, 1), tail(cs_Ac$median, 1),
    tail(cs_Ac$median, 1) - tail(cs_A$median, 1)
  ))
} else {
  cat("  (baseline fit_A.rds not found; skipped)\n")
}

cat("\n[benchmark] done.\n")
