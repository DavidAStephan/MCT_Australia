# Marginal Metropolis-Hastings update for the AR(1) coefficient rho.
#
# The conditional Gibbs update (.update_rho in gibbs_sweep.R) is
# systematically biased: each FFBS draw of c adds posterior noise that
# is approximately iid relative to the AR(1) signal, deflating the
# lag-1 autocorrelation of any individual c-draw. Updating rho from the
# regression of c[t] on c[t-1] then biases rho downward — confirmed in
# scripts/diagnose_rho_bias.R (truth 0.7 ⇒ posterior median ~0.42).
#
# The fix: propose rho ~ N(current, prop_sd), accept based on the
# MARGINAL likelihood p(y | rho, other_params), which is the Kalman
# filter log-likelihood with c integrated out. This is exact (no
# conditional-c approximation) and recovers the correct rho posterior.
#
# Cost: one additional KF forward pass per Gibbs iteration. On real
# AU data at T=435 N=11, that's ~30-50 ms — Gibbs goes from ~100 ms
# to ~150 ms per iter, still 10x+ faster than Stan.

#' Marginal MH update for rho. Returns updated rho + acceptance flag.
#'
#' @param rho_cur   Current rho value.
#' @param y         T x N observation matrix (NA for missing).
#' @param obs_type  T x N integer indicator (0/1/2) — used to rebuild
#'                  the SSM at the proposed rho.
#' @param state     Current state list (sigma_c, sigma_s, sigma_eps,
#'                  lambda).
#' @param config    Config list (T_, N, var_s_init, etc.).
#' @param prop_sd   Proposal standard deviation for rho ~ N(cur, sd).
#' @param prior_mean,prior_sd Prior parameters (truncated to (-1, 1)).
#' @return List with `rho` (new value) and `accepted` (bool).
update_rho_marginal_mh <- function(rho_cur, y, obs_type, state, config,
                                   prop_sd = 0.05,
                                   prior_mean = 0.5, prior_sd = 0.3) {
  # Propose
  rho_prop <- rnorm(1, mean = rho_cur, sd = prop_sd)

  # Reject immediately if outside (-1, 1) (proposal density is symmetric
  # so this is just a hard reject — no MH correction needed)
  if (rho_prop <= -1 || rho_prop >= 1) {
    return(list(rho = rho_cur, accepted = FALSE))
  }

  # Build SSMs at current and proposed rho
  ssm_cur <- build_mct_ssm_mixed(
    rho = rho_cur, lambda = state$lambda,
    sigma_c = state$sigma_c, sigma_s = state$sigma_s,
    sigma_eps = state$sigma_eps,
    obs_type = obs_type, T_ = config$T_, N = config$N,
    var_s_init = config$var_s_init
  )
  ssm_prop <- build_mct_ssm_mixed(
    rho = rho_prop, lambda = state$lambda,
    sigma_c = state$sigma_c, sigma_s = state$sigma_s,
    sigma_eps = state$sigma_eps,
    obs_type = obs_type, T_ = config$T_, N = config$N,
    var_s_init = config$var_s_init
  )

  # Marginal log-likelihoods via KF (no smoother needed; just need ll)
  ll_cur  <- kalman_filter(t(y), ssm_cur)$log_likelihood
  ll_prop <- kalman_filter(t(y), ssm_prop)$log_likelihood

  # Prior contributions — Gaussian, then truncation to (-1, 1)
  # (truncation normaliser cancels because both are inside the support)
  log_prior_cur  <- dnorm(rho_cur,  prior_mean, prior_sd, log = TRUE)
  log_prior_prop <- dnorm(rho_prop, prior_mean, prior_sd, log = TRUE)

  log_ratio <- (ll_prop + log_prior_prop) - (ll_cur + log_prior_cur)

  if (log(runif(1)) < log_ratio) {
    return(list(rho = rho_prop, accepted = TRUE))
  }
  list(rho = rho_cur, accepted = FALSE)
}
