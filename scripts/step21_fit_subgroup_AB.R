# Step 21 — Fit Variants A and B via Gibbs on the sub-group panel.
#
# Sub-group panel: T=399, N=33, all quarterly observations. Per-iter
# cost scales roughly as N^3 for the Kalman filter / smoother — at
# N=33 vs the production N=11, that's ~27x. Budget: ~2-3 hours per
# fit at n_burn=1000, n_draw=2000.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

stan_data_sg <- readRDS("outputs/draws/stan_data_subgroup.rds")
T_ <- stan_data_sg$T; N <- stan_data_sg$N

# Re-materialise dense y / obs_type from the sparse stan_data lists.
y        <- matrix(NA_real_, T_, N)
obs_type <- matrix(0L,      T_, N)
for (k in seq_along(stan_data_sg$y_q)) {
  y[stan_data_sg$t_q[k], stan_data_sg$i_q[k]] <- stan_data_sg$y_q[k]
  obs_type[stan_data_sg$t_q[k], stan_data_sg$i_q[k]] <- 2L
}
# (No monthly obs in sub-group panel — n_obs_m = 0.)

cat(sprintf("[step21] Sub-group fit setup: T=%d N=%d  n_obs_q=%d\n",
            T_, N, stan_data_sg$n_obs_q))
cat("Reference sub-group (ref=", stan_data_sg$ref, "): ",
    stan_data_sg$groups[stan_data_sg$ref], "\n", sep = "")

# ----- Variant B (RW common trend — primary model) -----
fit_B_path <- "outputs/draws/fit_gibbs_B_subgroup.rds"
if (file.exists(fit_B_path)) {
  cat("\n[step21] B cached at ", fit_B_path, " — loading.\n", sep = "")
  fit_B <- readRDS(fit_B_path)
} else {
  cat("\n[step21] Fitting Variant B on sub-group panel...\n")
  set.seed(1L)
  t0 <- Sys.time()
  fit_B <- fit_mct_gibbs(
    y, obs_type = obs_type, ref = stan_data_sg$ref,
    n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
    config = list(variant = "B")
  )
  cat(sprintf("[step21] B done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit_B, fit_B_path)
  cat("[step21] saved ", fit_B_path, "\n", sep = "")
}

# ----- Variant A (AR(1) common cycle — comparison) -----
fit_A_path <- "outputs/draws/fit_gibbs_A_subgroup.rds"
if (file.exists(fit_A_path)) {
  cat("\n[step21] A cached at ", fit_A_path, " — loading.\n", sep = "")
  fit_A <- readRDS(fit_A_path)
} else {
  cat("\n[step21] Fitting Variant A on sub-group panel...\n")
  set.seed(2L)
  t0 <- Sys.time()
  fit_A <- fit_mct_gibbs(
    y, obs_type = obs_type, ref = stan_data_sg$ref,
    n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
    config = list(variant = "A",
                  use_marginal_mh_rho = FALSE)
  )
  cat(sprintf("[step21] A done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(fit_A, fit_A_path)
  cat("[step21] saved ", fit_A_path, "\n", sep = "")
}

# ----- Headline comparison -----
cat("\n=== Sub-group fit results ===\n")
cat(sprintf("Variant A rho (median, 68%%): %.3f  [%.3f, %.3f]\n",
            median(fit_A$draws$rho),
            quantile(fit_A$draws$rho, 0.16),
            quantile(fit_A$draws$rho, 0.84)))

# Wrap and compute headline metrics
fit_A_wrap <- mct_gibbs_fit(fit_A, stan_data_sg)
fit_B_wrap <- mct_gibbs_fit(fit_B, stan_data_sg)
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

trend_A <- trend_path(fit_A_wrap, stan_data_sg$dates)
trend_B <- trend_path(fit_B_wrap, stan_data_sg$dates)
cs_A    <- common_share_path(fit_A_wrap, stan_data_sg$dates)
cs_B    <- common_share_path(fit_B_wrap, stan_data_sg$dates)

cat(sprintf("\nLatest-t (%s) trend:\n",
            format(tail(stan_data_sg$dates, 1L), "%Y-%m")))
cat(sprintf("  Variant A:  median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(trend_A$median, 1L),
            tail(trend_A$lower, 1L), tail(trend_A$upper, 1L)))
cat(sprintf("  Variant B:  median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(trend_B$median, 1L),
            tail(trend_B$lower, 1L), tail(trend_B$upper, 1L)))

cat("\nLatest-t common-share:\n")
cat(sprintf("  Variant A:  median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(cs_A$median, 1L),
            tail(cs_A$lower, 1L), tail(cs_A$upper, 1L)))
cat(sprintf("  Variant B:  median=%.3f  68%%=[%.3f, %.3f]\n",
            tail(cs_B$median, 1L),
            tail(cs_B$lower, 1L), tail(cs_B$upper, 1L)))

# 2-way LOO comparison
cat("\n=== LOO comparison ===\n")
ll_A <- as.matrix(fit_A_wrap$draws(variables = "log_lik"))
ll_B <- as.matrix(fit_B_wrap$draws(variables = "log_lik"))
keep <- (colSums(abs(ll_A)) > 0) & (colSums(abs(ll_B)) > 0)
ll_A <- ll_A[, keep, drop = FALSE]
ll_B <- ll_B[, keep, drop = FALSE]

r_eff_A <- loo::relative_eff(exp(ll_A), chain_id = rep(1L, nrow(ll_A)))
r_eff_B <- loo::relative_eff(exp(ll_B), chain_id = rep(1L, nrow(ll_B)))
loo_A <- loo::loo(ll_A, r_eff = r_eff_A)
loo_B <- loo::loo(ll_B, r_eff = r_eff_B)
print(loo::loo_compare(list(A = loo_A, B = loo_B)))

saveRDS(list(A = loo_A, B = loo_B),
        "outputs/draws/loo_subgroup_AB.rds")
cat("\n[step21] done.\n")
