# Gibbs sweep with MA(q) measurement errors (monthly-only mode).
#
# Step 12 of gibbs_port_plan.md. Adds:
#   - Augmented state with (q+1) eps lags per sector
#   - update_theta for the MA coefficients
#   - eps-derived residuals for SV update of sigma_eps
#
# State convention (NY Fed): eps state's variance = sigma_eps^2 (the
# sigma is folded into state innovation noise, NOT the obs equation H).
# This makes the SV update on sigma_eps directly applicable: the
# smoothed eps draws ARE the scaled noise x_t in x_t = sigma_t * unit_eps.
#
# Limitations of v1:
#   - Monthly observations only (no quarterly mixing). The state dim
#     would be 1 + N + (q+1)*N + extra for quarterly lags — leave that
#     to a follow-up.
#   - Constant lambda. TVP lambda was ported separately in
#     update_tvcoef.R, not wired in here.

#' One Gibbs sweep with MA(q) errors. Returns updated `state`.
#'
#' @param y       T x N observation matrix (no NA — monthly-only mode).
#' @param state   Current state list. Must include theta (N x q matrix)
#'                in addition to the standard fields.
#' @param config  Config list with extra fields q_MA and theta_prec_prior.
#' @return Updated state.
gibbs_sweep_ma <- function(y, state, config) {
  T_   <- config$T_
  N    <- config$N
  q_MA <- config$q_MA

  # --- 1. Build augmented-state SSM ---------------------------------
  ssm <- build_mct_ssm_ma(
    rho = state$rho, lambda = state$lambda,
    sigma_c = state$sigma_c, sigma_s = state$sigma_s,
    sigma_eps = state$sigma_eps,
    theta = state$theta, T_ = T_, N = N, q_MA = q_MA,
    var_s_init = config$var_s_init
  )

  # --- 2. Joint draw (c, s, eps_lags) via FFBS ----------------------
  Y <- t(y)
  smooth <- simulation_smoother(Y, ssm)
  states <- smooth$states                    # D x T_

  # --- 3. Extract c, s_i, and eps_{i, t..t-q} from augmented state -
  meta <- ssm$.meta
  c_path <- as.numeric(states[meta$c_pos, ])
  s_path <- t(states[meta$s_pos, , drop = FALSE])    # T_ x N

  # eps for each sector: T_ x (q_MA + 1) matrix (current + q lags at each t)
  # We need the "current" eps for SV update: positions meta$eps_start[i] across t
  eps_current <- matrix(0, T_, N)
  for (i in seq_len(N)) {
    eps_current[, i] <- as.numeric(states[meta$eps_start[i], ])
  }

  # --- 4. SV updates for c, s_i, eps_i ------------------------------
  nu_c <- c_path[-1L] - state$rho * c_path[-T_]
  nu_c_full <- c(c_path[1L] * sqrt(max(1 - state$rho^2, 1e-6)), nu_c)
  state$sigma_c <- update_vol(nu_c_full, state$sigma_c, state$gamma_c)

  nu_s <- rbind(s_path[1L, ],
                s_path[-1L, , drop = FALSE] - s_path[-T_, , drop = FALSE])
  for (i in seq_len(N)) {
    state$sigma_s[, i] <- update_vol(nu_s[, i], state$sigma_s[, i],
                                     state$gamma_s[i])
  }
  # eps_current is the "scaled noise" sigma_eps * unit_eps at each t.
  # Directly applicable to update_vol.
  for (i in seq_len(N)) {
    state$sigma_eps[, i] <- update_vol(eps_current[, i], state$sigma_eps[, i],
                                       state$gamma_eps[i])
  }

  # --- 5. Gamma updates --------------------------------------------
  state$gamma_c <- update_gam(
    2 * diff(log(state$sigma_c)),
    nu_prior = config$nu_gam_prior, s2_prior = config$s2_gam_prior
  )
  state$gamma_s <- update_gam(
    2 * diff(log(state$sigma_s)),
    nu_prior = rep(config$nu_gam_prior, N),
    s2_prior = rep(config$s2_gam_prior, N)
  )
  state$gamma_eps <- update_gam(
    2 * diff(log(state$sigma_eps)),
    nu_prior = rep(config$nu_gam_prior, N),
    s2_prior = rep(config$s2_gam_prior, N)
  )

  # --- 6. rho update (conditional Gibbs — MH not yet wired for MA path)
  state$rho <- .update_rho(
    c_path, state$sigma_c,
    prior_mean = config$rho_prior_mean,
    prior_sd   = config$rho_prior_sd
  )

  # --- 7. Lambda update (constant, conjugate) ----------------------
  # Use the "denoised" obs (no MA contribution): obs - MA_contribution
  # = lambda_i * c_t + s_{i,t}. With current eps draws we can subtract
  # the MA piece from y to get a clean lambda regression.
  ma_residual <- matrix(0, T_, N)
  for (i in seq_len(N)) {
    p <- meta$eps_start[i]
    # MA contribution at obs (t, i):
    #   eps_{i,t} + theta_1*eps_{i,t-1} + ... + theta_q*eps_{i,t-q}
    #   = sum over l=0..q of theta_l * states[p+l, t] (with theta_0 = 1)
    for (l in 0:q_MA) {
      coef <- if (l == 0L) 1 else state$theta[i, l]
      ma_residual[, i] <- ma_residual[, i] + coef * as.numeric(states[p + l, ])
    }
  }
  y_denoised <- y - ma_residual
  state$lambda <- .update_lambda(
    y_denoised, c_path, s_path,
    sigma_eps = matrix(1, T_, N),     # noise has been removed; use unit sigma
    ref        = config$ref,
    prior_mean = config$lambda_prior_mean,
    prior_sd   = config$lambda_prior_sd
  )

  # --- 8. Theta update (per sector, conjugate Normal + invertibility)
  # Regression: u_t = sum_l theta_l * lag_l + unit_noise, where
  #   u_t = (y_obs_t - lambda_i * c_t - s_{i,t}) / sigma_eps_{i,t}
  #   lag_l = eps_state_{i, t-l} / sigma_eps_{i, t}   (NY Fed convention)
  #   unit_noise = eps_state_{i, t} / sigma_eps_{i, t} ~ N(0, 1)
  theta_new <- matrix(0, N, q_MA)
  for (i in seq_len(N)) {
    p <- meta$eps_start[i]
    sigma_t <- state$sigma_eps[, i]
    # u_t = full MA-noise residual / sigma_t
    u_t <- (y[, i] - state$lambda[i] * c_path - s_path[, i]) / sigma_t
    # Lagged eps state values, rescaled by current sigma
    x_th <- matrix(0, T_, q_MA)
    for (l in seq_len(q_MA)) {
      lag_eps <- numeric(T_)
      for (t in seq_len(T_)) {
        if (t - l >= 1L) lag_eps[t] <- as.numeric(states[p + l, t])
      }
      x_th[, l] <- lag_eps / sigma_t
    }
    # Drop first q_MA rows where lags are zero (incomplete data window)
    keep <- (q_MA + 1L):T_
    theta_new[i, ] <- update_theta(
      u_t[keep], x_th[keep, , drop = FALSE],
      prec_prior = config$theta_prec_prior
    )
  }
  state$theta <- theta_new

  state$c <- c_path
  state$s <- s_path
  state
}
