library(rstan)
m_kalman_non_centered <- stan_model("code/eg2-three-indicator-ar1/dsem_kalman_non_centered.stan")
saveRDS(m_kalman_non_centered, "code/eg2-three-indicator-ar1/dsem_kalman_non_centered.rds")
