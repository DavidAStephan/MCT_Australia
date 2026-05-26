# Compute log p(y | rho) for a grid of rho values with all other
# params at truth. If the log-lik peaks NEAR rho=0.7, my MH is fine
# (must be inferring posterior properly). If it peaks at ~0.45, there's
# a bug in the SSM or KF.

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
state_truth <- list(
  lambda = rep(1, N), sigma_c = rep(sc, T_),
  sigma_s = matrix(ss, T_, N), sigma_eps = matrix(se, T_, N)
)
cfg <- list(T_ = T_, N = N, var_s_init = 0.25)

cat("Log-lik scan over rho (other params at truth):\n")
rhos <- seq(0.05, 0.95, by = 0.05)
ll <- numeric(length(rhos))
for (k in seq_along(rhos)) {
  ssm <- build_mct_ssm_mixed(
    rho = rhos[k], lambda = state_truth$lambda,
    sigma_c = state_truth$sigma_c, sigma_s = state_truth$sigma_s,
    sigma_eps = state_truth$sigma_eps,
    obs_type = obs_type, T_ = T_, N = N,
    var_s_init = cfg$var_s_init
  )
  ll[k] <- kalman_filter(t(y), ssm)$log_likelihood
}
for (k in seq_along(rhos)) {
  cat(sprintf("  rho = %.2f  log_lik = %.2f\n", rhos[k], ll[k]))
}
cat(sprintf("\nMax log-lik at rho = %.2f\n", rhos[which.max(ll)]))
