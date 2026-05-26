# Simulate from my SSM via simulate_ssm; check whether the implied
# c-path has the expected AR(1) coefficient. If not, there's a bug
# in the SSM construction.

suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

T_ <- 500; N <- 5; rho_true <- 0.7
sc <- 0.4; ss <- 0.25; se <- 0.8

obs_type <- matrix(1L, T_, N)
ssm <- build_mct_ssm_mixed(
  rho = rho_true, lambda = rep(1, N),
  sigma_c = rep(sc, T_), sigma_s = matrix(ss, T_, N),
  sigma_eps = matrix(se, T_, N),
  obs_type = obs_type, T_ = T_, N = N, var_s_init = 0.25
)

cat("SSM check:\n")
cat("  State dim D =", nrow(ssm$F), "\n")
cat("  F[1, 1:3] =", ssm$F[1, 1:3], "  (should be rho=0.7, 0, 0)\n")
cat("  F[2, 1:3] =", ssm$F[2, 1:3], "  (should be 1, 0, 0 — c lag shift)\n")
cat("  F[3, 1:3] =", ssm$F[3, 1:3], "  (should be 0, 1, 0)\n")
cat("  F[4, 4]   =", ssm$F[4, 4], "  (should be 1 — s_1 carries)\n")
cat("  Sigma_eta[1,1,1] =", ssm$Sigma_eta[1,1,1], "  (should be sc^2=0.16)\n")
cat("  G[1,1] =", ssm$G[1,1], "  (should be 1)\n")
cat("  G[4,2] =", ssm$G[4,2], "  (should be 1 — s_1 gets innovation 2)\n")
cat("  H[1, 1] =", ssm$H[1, 1, 1], "  (should be lambda[1]=1)\n")
cat("  H[1, 4] =", ssm$H[1, 4, 1], "  (should be 1 — s_1 enters obs 1)\n")

cat("\nSimulate from SSM, extract c-path, check empirical rho:\n")
set.seed(99)
sim <- simulate_ssm(ssm, T_ = T_, need_states = TRUE)
c_path <- as.numeric(sim$states[1L, ])
# Empirical AR(1) coef
rho_emp <- cor(c_path[-1], c_path[-T_])
cat(sprintf("  Empirical rho(c) from sim = %.3f (truth %.2f)\n",
            rho_emp, rho_true))
cat(sprintf("  Var(c) = %.3f  (theory sigma_c^2/(1-rho^2) = %.3f)\n",
            var(c_path), sc^2 / (1 - rho_true^2)))

# Also: extract s_1 path, check it's RW (rho(diff(s_1)) ~ 0)
s1_path <- as.numeric(sim$states[4L, ])
cat(sprintf("  Var(diff(s_1)) = %.3f  (theory sigma_s^2 = %.3f)\n",
            var(diff(s1_path)), ss^2))

# Sanity: does the SSM simulator give the same y distribution as the
# direct DGP?
set.seed(100)
direct_c <- numeric(T_)
direct_c[1] <- sc / sqrt(1 - rho_true^2) * rnorm(1)
for (t in 2:T_) direct_c[t] <- rho_true * direct_c[t-1] + sc * rnorm(1)
cat(sprintf("  Direct DGP: empirical rho(c) = %.3f\n",
            cor(direct_c[-1], direct_c[-T_])))
