# Definitive diagnostic: pin all params at TRUTH and run only the
# rho MH update. If rho recovers, the issue is the joint Gibbs chain
# drifting. If rho doesn't recover even with truth, there's a bug
# in the marginal MH itself or in the SSM/KF.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

# Simulate
set.seed(51); T_ <- 200; N <- 5; rho_true <- 0.7
sc <- 0.4; ss <- 0.25; se <- 0.8
cp <- numeric(T_); cp[1] <- sc / sqrt(1 - rho_true^2) * rnorm(1)
for (t in 2:T_) cp[t] <- rho_true * cp[t - 1] + sc * rnorm(1)
sp <- matrix(0, T_, N)
for (t in 2:T_) sp[t, ] <- sp[t - 1, ] + ss * rnorm(N)
y <- matrix(0, T_, N)
for (t in 1:T_) y[t, ] <- cp[t] + sp[t, ] + se * rnorm(N)

obs_type <- matrix(1L, T_, N)
state_truth <- list(
  rho       = 0.5,                          # initial; MH will move it
  lambda    = rep(1, N),                    # truth
  sigma_c   = rep(sc, T_),                  # truth (constant)
  sigma_s   = matrix(ss, T_, N),            # truth (constant)
  sigma_eps = matrix(se, T_, N)             # truth (constant)
)
cfg <- list(T_ = T_, N = N, var_s_init = 0.25)

# Run pure MH for rho with everything else pinned at truth
cat("=== rho MH only, everything else PINNED AT TRUTH ===\n")
for (psd in c(0.05, 0.10)) {
  cat(sprintf("\n  prop_sd = %.2f\n", psd))
  rho_draws <- numeric(2000)
  rho_cur <- 0.5
  n_acc <- 0
  for (i in seq_along(rho_draws)) {
    step <- update_rho_marginal_mh(
      rho_cur, y, obs_type, state_truth, cfg,
      prop_sd = psd, prior_mean = 0.5, prior_sd = 0.3
    )
    rho_cur <- step$rho
    if (step$accepted) n_acc <- n_acc + 1
    rho_draws[i] <- rho_cur
  }
  # Drop first 500 as burn-in
  rho_post <- rho_draws[501:2000]
  cat(sprintf("    rho median = %.3f  [%.3f, %.3f]   truth = %.2f\n",
              median(rho_post),
              quantile(rho_post, 0.05),
              quantile(rho_post, 0.95), rho_true))
  cat(sprintf("    acceptance = %.1f%%\n", 100 * n_acc / 2000))
}
