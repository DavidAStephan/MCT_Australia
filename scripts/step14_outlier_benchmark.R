# Step 14 of gibbs_port_plan.md — Real-data outlier benchmark.
#
# Runs two mixed-frequency Gibbs fits on the cached AU CPI data:
#   (a) baseline (use_outliers = FALSE)
#   (b) outlier scale-mixture (use_outliers = TRUE)
# Then reports the posterior median s_outlier flags for sanity check,
# and runs a LOO comparison between the two fits.
#
# Usage:
#   Rscript scripts/step14_outlier_benchmark.R [n_burn] [n_draw]
# Defaults: n_burn = 1000, n_draw = 2000.

suppressMessages({
  library(MASS)
})

args <- commandArgs(trailingOnly = TRUE)
n_burn <- if (length(args) >= 1L) as.integer(args[[1L]]) else 1000L
n_draw <- if (length(args) >= 2L) as.integer(args[[2L]]) else 2000L

# Load all gibbs sources
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

sd_real <- readRDS("outputs/draws/stan_data.rds")
T_ <- sd_real$T; N <- sd_real$N

# Re-materialise the dense y / obs_type matrices from sparse stan_data.
# stan_data was produced with the sparse-list layout the Stan model
# needs; the Gibbs mixed-freq sweep wants dense (T x N).
y        <- matrix(NA_real_, T_, N)
obs_type <- matrix(0L,      T_, N)
for (k in seq_along(sd_real$y_m)) {
  y[sd_real$t_m[k], sd_real$i_m[k]] <- sd_real$y_m[k]
  obs_type[sd_real$t_m[k], sd_real$i_m[k]] <- 1L
}
for (k in seq_along(sd_real$y_q)) {
  y[sd_real$t_q[k], sd_real$i_q[k]] <- sd_real$y_q[k]
  obs_type[sd_real$t_q[k], sd_real$i_q[k]] <- 2L
}

dir.create("outputs/draws", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

cat("=== Step 14 benchmark === T=", T_, " N=", N,
    " n_burn=", n_burn, " n_draw=", n_draw, "\n", sep = "")
cat("Monthly obs: ", length(sd_real$y_m),
    "  Quarterly obs: ", length(sd_real$y_q), "\n", sep = "")

# ------------------------------------------------------------------
# (a) Baseline mixed-freq Gibbs (no outliers)
# ------------------------------------------------------------------
baseline_path <- "outputs/draws/fit_gibbs_baseline.rds"
if (file.exists(baseline_path)) {
  cat("\n[baseline] Cached fit found at ", baseline_path, " — loading.\n",
      sep = "")
  fit_base <- readRDS(baseline_path)
} else {
  cat("\n[baseline] Fitting Gibbs (use_outliers = FALSE) ...\n")
  set.seed(1L)
  t0 <- Sys.time()
  fit_base <- fit_mct_gibbs(
    y, obs_type = obs_type, ref = sd_real$ref,
    n_burn = n_burn, n_draw = n_draw, verbose = TRUE,
    config = list(use_outliers        = FALSE,
                  use_marginal_mh_rho = FALSE)
  )
  cat(sprintf("[baseline] elapsed = %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit_base, baseline_path)
  cat("[baseline] saved to ", baseline_path, "\n", sep = "")
}

# ------------------------------------------------------------------
# (b) Outlier-enabled mixed-freq Gibbs
# ------------------------------------------------------------------
outlier_path <- "outputs/draws/fit_gibbs_outliers.rds"
if (file.exists(outlier_path)) {
  cat("\n[outliers] Cached fit found at ", outlier_path, " — loading.\n",
      sep = "")
  fit_out <- readRDS(outlier_path)
} else {
  cat("\n[outliers] Fitting Gibbs (use_outliers = TRUE) ...\n")
  set.seed(2L)
  t0 <- Sys.time()
  fit_out <- fit_mct_gibbs(
    y, obs_type = obs_type, ref = sd_real$ref,
    n_burn = n_burn, n_draw = n_draw, verbose = TRUE,
    config = list(use_outliers        = TRUE,
                  use_marginal_mh_rho = FALSE)
  )
  cat(sprintf("[outliers] elapsed = %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit_out, outlier_path)
  cat("[outliers] saved to ", outlier_path, "\n", sep = "")
}

# ------------------------------------------------------------------
# Outlier flag inspection
# ------------------------------------------------------------------
cat("\n=== Outlier flags: per-sector top-5 dates (posterior median s) ===\n")
s_med <- apply(fit_out$draws$s_outlier, c(1, 2), median)

flag_rows <- vector("list", N)
for (i in seq_len(N)) {
  top <- order(-s_med[, i])[1:5]
  flag_rows[[i]] <- data.frame(
    sector = sd_real$groups[i],
    rank   = 1:5,
    date   = format(sd_real$dates[top], "%Y-%m"),
    s_med  = round(s_med[top, i], 2),
    obs_type = obs_type[cbind(top, rep(i, 5))]
  )
  cat(sprintf("%-50s ", sd_real$groups[i]))
  cat(paste(format(sd_real$dates[top], "%Y-%m"),
            sprintf("(s=%.1f)", s_med[top, i])),
      sep = "  ")
  cat("\n")
}
flag_df <- do.call(rbind, flag_rows)
write.csv(flag_df, "outputs/tables/outlier_top5_per_sector.csv",
          row.names = FALSE)

# Fraction of obs flagged as outliers (s_med > 1.5 as a rough cut)
n_obs_total <- sum(obs_type > 0L)
n_flagged   <- sum(s_med > 1.5 & obs_type > 0L)
cat(sprintf("\nObservations flagged (s_med > 1.5): %d / %d (%.1f%%)\n",
            n_flagged, n_obs_total, 100 * n_flagged / n_obs_total))

# Sector posterior-mean of ps (probability of normal obs)
cat("\nPosterior mean ps (probability of normal obs) per sector:\n")
ps_mean <- rowMeans(fit_out$draws$ps)
for (i in seq_len(N)) {
  cat(sprintf("  %-50s %.3f\n", sd_real$groups[i], ps_mean[i]))
}

cat("\n=== Step 14 complete ===\n")
cat("Files written:\n")
cat("  outputs/draws/fit_gibbs_baseline.rds\n")
cat("  outputs/draws/fit_gibbs_outliers.rds\n")
cat("  outputs/tables/outlier_top5_per_sector.csv\n")
