# Mixed-frequency Gibbs sweep. Step 9a of gibbs_port_plan.md.
#
# Difference from gibbs_sweep() (Step 7, monthly-only):
#   - Uses build_mct_ssm_mixed() with augmented 3(N+1)-dim state to
#     handle the quarterly-average observation equation in a single
#     linear step.
#   - Extracts c, s from the "current-t" positions of the augmented
#     state vector (positions 1 and 3i+1).
#   - Measurement-noise residuals eps[i, t] only computed at monthly
#     observation times; quarterly and missing positions are NA and
#     skipped by update_vol's NA fallback.
#   - Loadings lambda updated using monthly residuals only (v1
#     simplification; quarterly contributions can be added later).

# Lambda update using ONLY monthly residuals. Quarterly contributions
# could be added by including the quarterly-residual regression too,
# but that requires summing three lags of c — punted to a follow-up.
.update_lambda_mixed <- function(y, obs_type, c_path, s_path, sigma_eps,
                                 ref, prior_mean = 1, prior_sd = 0.3) {
  N      <- ncol(y)
  lambda <- numeric(N)
  lambda[ref] <- 1
  for (i in seq_len(N)) {
    if (i == ref) next
    mask <- (obs_type[, i] == 1L)
    if (!any(mask)) {
      # No monthly obs for this sector — fall back to prior
      lambda[i] <- rnorm(1, prior_mean, prior_sd)
      next
    }
    wy <- (y[mask, i] - s_path[mask, i]) / sigma_eps[mask, i]
    wx <- c_path[mask] / sigma_eps[mask, i]
    prec_post <- 1 / prior_sd^2 + sum(wx^2)
    mean_post <- (prior_mean / prior_sd^2 + sum(wx * wy)) / prec_post
    lambda[i] <- rnorm(1, mean_post, sqrt(1 / prec_post))
  }
  lambda
}

#' One Gibbs sweep for the mixed-frequency MCT Variant A model.
#'
#' @param y         T x N matrix; obs values where present, NA where
#'                  missing. Quarterly values appear at the last month
#'                  of each quarter, monthly values at every observed t.
#' @param obs_type  T x N integer matrix; 0 = missing, 1 = monthly,
#'                  2 = quarterly average (requires t >= 3).
#' @param state     Current state list (same shape as gibbs_sweep).
#' @param config    Config list (same shape as gibbs_sweep, plus
#'                  no new entries).
#' @return Updated state list.
gibbs_sweep_mixed <- function(y, obs_type, state, config) {
  T_ <- config$T_
  N  <- config$N
  use_outliers <- isTRUE(config$use_outliers)

  # --- 1. Build augmented-state SSM ---------------------------------
  # When outliers enabled, the EFFECTIVE measurement-noise SD is
  # sigma_eps * s_outlier (per-obs scale inflation). This is what
  # the KF sees as obs noise.
  effective_sigma_eps <- if (use_outliers) {
    state$sigma_eps * state$s_outlier
  } else {
    state$sigma_eps
  }

  ssm <- build_mct_ssm_mixed(
    rho = state$rho, lambda = state$lambda,
    sigma_c = state$sigma_c, sigma_s = state$sigma_s,
    sigma_eps = effective_sigma_eps,
    obs_type = obs_type, T_ = T_, N = N,
    var_s_init = config$var_s_init
  )

  # --- 2. Joint draw of augmented state via FFBS --------------------
  # Pass Y with NA for missing obs (kalman_filter handles those).
  Y <- t(y)
  smooth <- simulation_smoother(Y, ssm)
  states <- smooth$states                     # D x T_

  # --- 3. Extract c (pos 1) and s_i (pos 3i+1) ----------------------
  c_path <- as.numeric(states[1L, ])
  s_pos  <- 3L * seq_len(N) + 1L              # positions of s_{i,t}
  s_path <- t(states[s_pos, , drop = FALSE])  # T_ x N

  # --- 4. SV updates ------------------------------------------------
  nu_c <- c_path[-1L] - state$rho * c_path[-T_]
  nu_c_full <- c(c_path[1L] * sqrt(max(1 - state$rho^2, 1e-6)), nu_c)
  state$sigma_c <- update_vol(nu_c_full, state$sigma_c, state$gamma_c)

  nu_s <- rbind(s_path[1L, ],
                s_path[-1L, , drop = FALSE] - s_path[-T_, , drop = FALSE])
  for (i in seq_len(N)) {
    state$sigma_s[, i] <- update_vol(nu_s[, i], state$sigma_s[, i],
                                     state$gamma_s[i])
  }

  # Measurement-noise residuals: ONLY at monthly observations.
  # Quarterly + missing positions get NA so update_vol's missing-data
  # fallback kicks in (no update; propagate the RW state).
  c_mat   <- matrix(c_path, T_, N)
  lam_mat <- matrix(state$lambda, T_, N, byrow = TRUE)
  eps_mat <- y - lam_mat * c_mat - s_path
  eps_mat[obs_type != 1L] <- NA

  if (use_outliers) {
    # Step A: update per-obs scale s_outlier given residuals + current
    # sigma_eps. Whitened residual w_t = residual_t / sigma_eps_t should
    # be ~ N(0, s_outlier_t^2).
    for (i in seq_len(N)) {
      probs_i <- c(state$ps[i],
                   rep((1 - state$ps[i]) /
                       (length(config$s_vals) - 1L),
                       length(config$s_vals) - 1L))
      w_res_i <- eps_mat[, i] / state$sigma_eps[, i]
      state$s_outlier[, i] <- update_scl(w_res_i, config$s_vals, probs_i)
    }

    # Step B: update sigma_eps using residuals DIVIDED BY scale
    # (so eps_t / s_t = sigma_t * unit_eps_t matches update_vol's model).
    for (i in seq_len(N)) {
      state$sigma_eps[, i] <- update_vol(
        eps_mat[, i] / state$s_outlier[, i],
        state$sigma_eps[, i], state$gamma_eps[i]
      )
    }

    # Step C: update ps (probability of normal obs) given current
    # s_outlier indicators. Use only monthly-obs positions for the
    # Bernoulli count.
    indicator <- (state$s_outlier == config$s_vals[1L]) & (obs_type == 1L)
    # Restrict to rows with any monthly obs (cleaner counts)
    state$ps <- update_ps(indicator, config$ps_a_prior, config$ps_b_prior)
  } else {
    for (i in seq_len(N)) {
      state$sigma_eps[, i] <- update_vol(eps_mat[, i], state$sigma_eps[, i],
                                         state$gamma_eps[i])
    }
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

  # --- 6. rho update ------------------------------------------------
  # Marginal MH is the unbiased path; conditional Gibbs has a known
  # FFBS-noise downward bias (see update_rho_mh.R). Config flag lets
  # us A/B test if needed.
  if (isTRUE(config$use_marginal_mh_rho)) {
    rho_step <- update_rho_marginal_mh(
      state$rho, y, obs_type, state, config,
      prop_sd    = config$rho_prop_sd,
      prior_mean = config$rho_prior_mean,
      prior_sd   = config$rho_prior_sd
    )
    state$rho             <- rho_step$rho
    state$last_rho_accept <- rho_step$accepted
  } else {
    state$rho <- .update_rho(
      c_path, state$sigma_c,
      prior_mean = config$rho_prior_mean,
      prior_sd   = config$rho_prior_sd
    )
  }

  # --- 7. Lambda update (monthly residuals only for v1) -------------
  state$lambda <- .update_lambda_mixed(
    y, obs_type, c_path, s_path, state$sigma_eps,
    ref        = config$ref,
    prior_mean = config$lambda_prior_mean,
    prior_sd   = config$lambda_prior_sd
  )

  state$c <- c_path
  state$s <- s_path
  state
}
