# Step 10 of gibbs_port_plan.md: smoke-test fit_mct_gibbs on REAL
# ABS data (T=435, N=11, mixed-freq).
#
# Goal: confirm the Gibbs sampler runs end-to-end on production-scale
# data and produces sensible posterior medians comparable to the Stan
# Variant A baseline (cached fit at outputs/draws/fit_A.rds).
#
# Short chain (100 burn + 100 draws) first to time it. If reasonable,
# extend to full production.

suppressMessages({
  library(MASS); library(posterior); library(dplyr)
})
for (f in list.files("R", "\\.R$", full.names = TRUE, recursive = TRUE)) {
  source(f)
}

cat("[step10] loading real stan_data ...\n")
sd_real <- readRDS("outputs/draws/stan_data.rds")
T_ <- sd_real$T; N <- sd_real$N

# Build y matrix + obs_type matrix from sparse stan_data lists.
y <- matrix(NA_real_, T_, N)
obs_type <- matrix(0L, T_, N)
for (k in seq_along(sd_real$y_m)) {
  y[sd_real$t_m[k], sd_real$i_m[k]] <- sd_real$y_m[k]
  obs_type[sd_real$t_m[k], sd_real$i_m[k]] <- 1L
}
for (k in seq_along(sd_real$y_q)) {
  y[sd_real$t_q[k], sd_real$i_q[k]] <- sd_real$y_q[k]
  obs_type[sd_real$t_q[k], sd_real$i_q[k]] <- 2L
}
cat(sprintf("[step10] y/obs_type built: %d monthly, %d quarterly, %d missing\n",
            sum(obs_type == 1L), sum(obs_type == 2L),
            sum(obs_type == 0L)))

# Smoke test showed ~0.09 sec/iter on this hardware. Run a fuller
# chain (1000 burn + 2000 draws ≈ 4.5 min) so the comparison to the
# Stan baseline is meaningful.
n_burn  <- 1000L
n_draw  <- 2000L
cat(sprintf("\n[step10] fit_mct_gibbs (mixed-freq) — %d burn + %d draws ...\n",
            n_burn, n_draw))
t0 <- Sys.time()
fit <- fit_mct_gibbs(y, obs_type = obs_type, ref = sd_real$ref,
                     n_burn = n_burn, n_draw = n_draw, verbose = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\n[step10] done in %.1f sec (%.2f sec/iter)\n",
            elapsed, elapsed / (n_burn + n_draw)))

# Inspect posterior medians for headline parameters
cat("\n[step10] posterior medians:\n")
cat(sprintf("  rho                  : %.3f  [%.3f, %.3f]\n",
            median(fit$draws$rho),
            quantile(fit$draws$rho, 0.05),
            quantile(fit$draws$rho, 0.95)))
lambda_med <- apply(fit$draws$lambda, 1, median)
cat(sprintf("  lambda[1..N] medians : %s\n",
            paste(sprintf("%.2f", lambda_med), collapse = ", ")))

if (file.exists("outputs/draws/fit_A.rds")) {
  cat("\n[step10] comparing to Stan Variant A baseline:\n")
  fit_A <- readRDS("outputs/draws/fit_A.rds")

  # Stan latest-t trend draws (column T_ of the trend matrix)
  trend_stan <- as.numeric(
    fit_A$draws(variables = "trend", format = "draws_matrix")[, T_])

  # Gibbs latest-t trend: trend[t, draw] = sum_i w_i * s[t, i, draw]
  # draws$s is (T_, N, n_keep). Take latest t → (N x n_keep), transpose
  # → (n_keep x N), multiply by w → length n_keep.
  s_T <- fit$draws$s[T_, , ]                       # N x n_keep
  trend_gibbs <- as.numeric(t(s_T) %*% sd_real$w)  # length n_keep

  cat(sprintf("  latest-t trend median  Stan: %.3f   Gibbs: %.3f\n",
              median(trend_stan), median(trend_gibbs)))
  cat(sprintf("  latest-t trend 68%%      Stan: [%.2f, %.2f]   Gibbs: [%.2f, %.2f]\n",
              quantile(trend_stan, 0.16), quantile(trend_stan, 0.84),
              quantile(trend_gibbs, 0.16), quantile(trend_gibbs, 0.84)))

  # Compare trend posterior medians across the WHOLE history
  trend_stan_med  <- apply(
    fit_A$draws(variables = "trend", format = "draws_matrix"),
    2, median
  )
  trend_gibbs_med <- numeric(T_)
  for (t in seq_len(T_)) {
    trend_gibbs_med[t] <- median(as.numeric(
      t(fit$draws$s[t, , ]) %*% sd_real$w))
  }
  cor_full <- cor(trend_stan_med, trend_gibbs_med)
  cat(sprintf("  trend posterior median cor(Stan, Gibbs) across t: %.3f\n",
              cor_full))
  mae <- mean(abs(trend_stan_med - trend_gibbs_med))
  cat(sprintf("  trend posterior median MAE Stan vs Gibbs: %.3f\n", mae))
}

# Extrapolate full-run cost: standard production = 1500 burn + 1500 draws.
proj_full <- elapsed / (n_burn + n_draw) * (1500 + 1500)
cat(sprintf("\n[step10] extrapolated full-run (1500+1500): %.0f sec = %.1f min = %.1f hr\n",
            proj_full, proj_full / 60, proj_full / 3600))
