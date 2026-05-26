# Tight Gibbs vs Stan comparison on real ABS data, using a LONGER
# Gibbs chain (with the marginal-MH rho update for defensiveness)
# and reporting actual overlap of posterior CIs.

suppressMessages({
  library(MASS); library(posterior); library(dplyr)
})
for (f in list.files("R", "\\.R$", full.names = TRUE, recursive = TRUE)) {
  source(f)
}

cat("[step10b] loading real stan_data + cached Stan fit_A ...\n")
sd_real <- readRDS("outputs/draws/stan_data.rds")
fit_A   <- readRDS("outputs/draws/fit_A.rds")
T_ <- sd_real$T; N <- sd_real$N

# Build y/obs_type for the gibbs path
y        <- matrix(NA_real_, T_, N)
obs_type <- matrix(0L, T_, N)
for (k in seq_along(sd_real$y_m)) {
  y[sd_real$t_m[k], sd_real$i_m[k]] <- sd_real$y_m[k]
  obs_type[sd_real$t_m[k], sd_real$i_m[k]] <- 1L
}
for (k in seq_along(sd_real$y_q)) {
  y[sd_real$t_q[k], sd_real$i_q[k]] <- sd_real$y_q[k]
  obs_type[sd_real$t_q[k], sd_real$i_q[k]] <- 2L
}

# Run Gibbs with marginal-MH rho and a longer chain
n_burn <- 2000L; n_draw <- 5000L
cat(sprintf("\n[step10b] fit_mct_gibbs (marginal-MH rho, %d burn + %d draws)...\n",
            n_burn, n_draw))
t0 <- Sys.time()
fit_g <- fit_mct_gibbs(
  y, obs_type = obs_type, ref = sd_real$ref,
  n_burn = n_burn, n_draw = n_draw, verbose = FALSE,
  config = list(use_marginal_mh_rho = TRUE, rho_prop_sd = 0.05)
)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("[step10b] done in %.1f min (acceptance %.1f%%)\n",
            el, 100 * (fit_g$rho_accept_rate %||% NA)))

# ---- Compare rho, lambda, latest-t trend ----
cat("\n[step10b] === rho posterior ===\n")
rho_stan <- as.numeric(fit_A$draws(variables = "rho", format = "draws_matrix"))
rho_g    <- fit_g$draws$rho
report <- function(label, x) {
  cat(sprintf("  %-7s  median %.3f   68%% [%.3f, %.3f]   90%% [%.3f, %.3f]\n",
              label, median(x),
              quantile(x, 0.16), quantile(x, 0.84),
              quantile(x, 0.05), quantile(x, 0.95)))
}
report("Stan",  rho_stan)
report("Gibbs", rho_g)

cat("\n[step10b] === latest-t trend posterior (March 2026) ===\n")
trend_stan_T <- as.numeric(
  fit_A$draws(variables = "trend", format = "draws_matrix")[, T_])
trend_g_T <- as.numeric(t(fit_g$draws$s[T_, , ]) %*% sd_real$w)
report("Stan",  trend_stan_T)
report("Gibbs", trend_g_T)

cat("\n[step10b] === trend posterior MEDIAN across all t ===\n")
trend_stan_full <- apply(
  fit_A$draws(variables = "trend", format = "draws_matrix"),
  2, median
)
trend_g_full <- numeric(T_)
for (t in seq_len(T_)) {
  trend_g_full[t] <- median(as.numeric(t(fit_g$draws$s[t, , ]) %*% sd_real$w))
}
cat(sprintf("  cor(Stan, Gibbs) = %.3f    MAE = %.3f    range Stan: [%.2f, %.2f]\n",
            cor(trend_stan_full, trend_g_full),
            mean(abs(trend_stan_full - trend_g_full)),
            min(trend_stan_full), max(trend_stan_full)))
