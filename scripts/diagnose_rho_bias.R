# Probe whether rho bias is in:
#   A. .update_rho (the math) — test by feeding TRUTH c + sigma_c
#   B. upstream c-draws — test by inspecting posterior sigma_c against truth
#   C. somewhere else

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

sim_mct_a <- function(T_, N, rho = 0.7,
                      sigma_c_true = 0.4, sigma_s_true = 0.25,
                      sigma_eps_true = 0.8, seed = 1L) {
  set.seed(seed)
  sigma_c <- rep(sigma_c_true, T_)
  sigma_s <- matrix(sigma_s_true, T_, N)
  sigma_eps <- matrix(sigma_eps_true, T_, N)
  c_path <- numeric(T_)
  c_path[1L] <- sigma_c[1L] / sqrt(1 - rho^2) * rnorm(1)
  for (t in 2:T_) c_path[t] <- rho * c_path[t - 1L] + sigma_c[t] * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) s_path[t, ] <- s_path[t - 1L, ] + sigma_s[t, ] * rnorm(N)
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- c_path[t] + s_path[t, ] + sigma_eps[t, ] * rnorm(N)
  }
  list(y = y, c = c_path, s = s_path, rho = rho,
       sigma_c = sigma_c, sigma_s = sigma_s, sigma_eps = sigma_eps)
}

set.seed(51); T_ <- 200; N <- 5
sim <- sim_mct_a(T_ = T_, N = N, rho = 0.7, seed = 51)

# Direct probe A: feed TRUTH c and sigma_c into .update_rho. Should
# return rho draws concentrated near 0.7.
cat("=== A) .update_rho on TRUTH c + sigma_c ===\n")
rho_draws_truth <- replicate(1000, .update_rho(sim$c, sim$sigma_c))
cat(sprintf("  rho draws: median=%.3f  [q05=%.3f, q95=%.3f]\n",
            median(rho_draws_truth),
            quantile(rho_draws_truth, 0.05),
            quantile(rho_draws_truth, 0.95)))

# Probe B: run Gibbs and inspect sigma_c posterior vs truth (0.4).
cat("\n=== B) Gibbs run + sigma_c posterior median ===\n")
fit <- fit_mct_gibbs(sim$y, ref = 1L,
                    n_burn = 1000L, n_draw = 2000L, verbose = FALSE)
sigma_c_med <- apply(fit$draws$sigma_c, 1, median)
cat(sprintf("  truth sigma_c = %.3f (constant)\n", 0.4))
cat(sprintf("  posterior sigma_c median across t: mean=%.3f  min=%.3f  max=%.3f\n",
            mean(sigma_c_med), min(sigma_c_med), max(sigma_c_med)))

# Probe C: do the same with TRUTH c but POSTERIOR sigma_c — does
# rho update still recover 0.7?
sigma_c_chain <- fit$draws$sigma_c   # T_ x n_keep
cat("\n=== C) .update_rho on TRUTH c + each posterior sigma_c draw ===\n")
rho_with_truth_c <- numeric(ncol(sigma_c_chain))
for (k in seq_len(ncol(sigma_c_chain))) {
  rho_with_truth_c[k] <- .update_rho(sim$c, sigma_c_chain[, k])
}
cat(sprintf("  rho draws (truth-c, post-sigma): median=%.3f  q05=%.3f q95=%.3f\n",
            median(rho_with_truth_c),
            quantile(rho_with_truth_c, 0.05),
            quantile(rho_with_truth_c, 0.95)))

# Probe D: take each posterior c draw, run .update_rho with TRUTH
# sigma_c. Isolates whether c-draws have the right autocorrelation.
cat("\n=== D) .update_rho on each posterior c + TRUTH sigma_c ===\n")
c_chain <- fit$draws$c
rho_with_post_c <- numeric(ncol(c_chain))
for (k in seq_len(ncol(c_chain))) {
  rho_with_post_c[k] <- .update_rho(c_chain[, k], sim$sigma_c)
}
cat(sprintf("  rho draws (post-c, truth-sigma): median=%.3f  q05=%.3f q95=%.3f\n",
            median(rho_with_post_c),
            quantile(rho_with_post_c, 0.05),
            quantile(rho_with_post_c, 0.95)))
