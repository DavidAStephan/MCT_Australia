# Quick diagnostic: does smart initialisation + a longer chain fix the
# rho-mixing bias? Three variants on the same sim data:
#
#   1. default init (sigma_c = sigma_s = sigma_eps = 0.5), short chain
#   2. default init, LONG chain (n_burn=3000, n_draw=5000)
#   3. smart init (sigma_eps from sd(y); sigma_c smaller; sigma_s smaller)
#      with default-length chain
#
# Compare rho posterior medians and c-path correlation with truth.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

sim_mct_a <- function(T_, N, rho = 0.7, lambda = NULL,
                      sigma_c_true = 0.4, sigma_s_true = 0.25,
                      sigma_eps_true = 0.8, ref = 1L, seed = 1L) {
  set.seed(seed)
  if (is.null(lambda)) lambda <- rep(1, N)
  lambda[ref] <- 1
  sigma_c <- rep(sigma_c_true, T_)
  sigma_s <- matrix(sigma_s_true, T_, N)
  sigma_eps <- matrix(sigma_eps_true, T_, N)
  c_path <- numeric(T_)
  c_path[1L] <- sigma_c[1L] / sqrt(1 - rho^2) * rnorm(1)
  for (t in 2:T_) c_path[t] <- rho * c_path[t - 1L] + sigma_c[t] * rnorm(1)
  s_path <- matrix(0, T_, N)
  for (t in 2:T_) {
    s_path[t, ] <- s_path[t - 1L, ] + sigma_s[t, ] * rnorm(N)
  }
  y <- matrix(0, T_, N)
  for (t in seq_len(T_)) {
    y[t, ] <- lambda * c_path[t] + s_path[t, ] + sigma_eps[t, ] * rnorm(N)
  }
  list(y = y, c = c_path, s = s_path, rho = rho)
}

set.seed(51)
T_ <- 200; N <- 5
sim <- sim_mct_a(T_ = T_, N = N, rho = 0.7, seed = 51)

# Smart init: c0 = mean across sectors of demeaned y; sigma_eps from
# OLS residuals; sigma_c and sigma_s from path innovations.
make_smart_init <- function(y, ref = 1L) {
  T_ <- nrow(y); N <- ncol(y)
  c0 <- rowMeans(y)                    # rough common factor
  s0 <- y - matrix(c0, T_, N)           # rough sector trends
  init_sigma_eps <- 0.7 * matrix(apply(y, 2, sd), T_, N, byrow = TRUE)
  init_sigma_c <- rep(sd(diff(c0)), T_)
  init_sigma_s <- matrix(sd(diff(s0)), T_, N)
  list(
    c = c0, s = s0,
    sigma_c = init_sigma_c,
    sigma_s = init_sigma_s,
    sigma_eps = init_sigma_eps,
    gamma_c = 0.05,
    gamma_s = rep(0.05, N),
    gamma_eps = rep(0.05, N),
    rho = 0.5,
    lambda = rep(1, N)
  )
}

cat("=== Variant 1: default init, short chain (500+1000) ===\n")
t0 <- Sys.time()
fit1 <- fit_mct_gibbs(sim$y, ref = 1L,
                     n_burn = 500L, n_draw = 1000L, verbose = FALSE)
cat(sprintf("  time = %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  rho median = %.3f  (truth 0.7)\n", median(fit1$draws$rho)))
cat(sprintf("  rho q05/q95 = %.3f / %.3f\n",
            quantile(fit1$draws$rho, 0.05),
            quantile(fit1$draws$rho, 0.95)))
cat(sprintf("  cor(c_post_mean, truth) = %.3f\n",
            cor(rowMeans(fit1$draws$c), sim$c)))

cat("\n=== Variant 2: default init, LONG chain (3000+5000) ===\n")
t0 <- Sys.time()
fit2 <- fit_mct_gibbs(sim$y, ref = 1L,
                     n_burn = 3000L, n_draw = 5000L, verbose = FALSE)
cat(sprintf("  time = %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  rho median = %.3f  (truth 0.7)\n", median(fit2$draws$rho)))
cat(sprintf("  rho q05/q95 = %.3f / %.3f\n",
            quantile(fit2$draws$rho, 0.05),
            quantile(fit2$draws$rho, 0.95)))
cat(sprintf("  cor(c_post_mean, truth) = %.3f\n",
            cor(rowMeans(fit2$draws$c), sim$c)))

cat("\n=== Variant 3: smart init, short chain (500+1000) ===\n")
t0 <- Sys.time()
fit3 <- fit_mct_gibbs(sim$y, ref = 1L,
                     init = make_smart_init(sim$y),
                     n_burn = 500L, n_draw = 1000L, verbose = FALSE)
cat(sprintf("  time = %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  rho median = %.3f  (truth 0.7)\n", median(fit3$draws$rho)))
cat(sprintf("  rho q05/q95 = %.3f / %.3f\n",
            quantile(fit3$draws$rho, 0.05),
            quantile(fit3$draws$rho, 0.95)))
cat(sprintf("  cor(c_post_mean, truth) = %.3f\n",
            cor(rowMeans(fit3$draws$c), sim$c)))
