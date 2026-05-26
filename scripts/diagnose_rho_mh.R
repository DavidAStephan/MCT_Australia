# Diagnose: rho-MH acceptance + recovery across different prop_sd values
suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

set.seed(51); T_ <- 200; N <- 5; rho_true <- 0.7
sc <- 0.4; ss <- 0.25; se <- 0.8
cp <- numeric(T_); cp[1] <- sc / sqrt(1 - rho_true^2) * rnorm(1)
for (t in 2:T_) cp[t] <- rho_true * cp[t - 1] + sc * rnorm(1)
sp <- matrix(0, T_, N)
for (t in 2:T_) sp[t, ] <- sp[t - 1, ] + ss * rnorm(N)
y <- matrix(0, T_, N)
for (t in 1:T_) y[t, ] <- cp[t] + sp[t, ] + se * rnorm(N)

obs_type <- matrix(1L, T_, N)

for (psd in c(0.05, 0.10, 0.20, 0.30)) {
  cat(sprintf("\n=== prop_sd = %.2f ===\n", psd))
  t0 <- Sys.time()
  fit <- fit_mct_gibbs(
    y, obs_type = obs_type, ref = 1L,
    n_burn = 500L, n_draw = 1000L, verbose = FALSE,
    config = list(use_marginal_mh_rho = TRUE, rho_prop_sd = psd)
  )
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  rho median = %.3f  [%.3f, %.3f]   truth = %.2f\n",
              median(fit$draws$rho),
              quantile(fit$draws$rho, 0.05),
              quantile(fit$draws$rho, 0.95), rho_true))
  cat(sprintf("  acceptance = %.1f%%   elapsed = %.1f s\n",
              100 * fit$rho_accept_rate, el))
}
