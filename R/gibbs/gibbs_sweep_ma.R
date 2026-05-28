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
  T_           <- config$T_
  N            <- config$N
  q_MA         <- config$q_MA
  variant      <- if (is.null(config$variant)) "A" else config$variant
  use_outliers <- isTRUE(config$use_outliers)

  # --- 1. Build augmented-state SSM ---------------------------------
  # When use_outliers, the eps-state innovation variance is inflated
  # to (sigma_eps * s_outlier)^2 at each (t, i). With outliers off,
  # s_outlier is implicitly 1 and the SSM is unchanged.
  s_outlier_eff <- if (use_outliers) state$s_outlier else NULL
  ssm <- build_mct_ssm_ma(
    rho = state$rho, lambda = state$lambda,
    sigma_c = state$sigma_c, sigma_s = state$sigma_s,
    sigma_eps = state$sigma_eps,
    theta = state$theta, T_ = T_, N = N, q_MA = q_MA,
    var_s_init = config$var_s_init,
    variant    = variant,
    s_outlier  = s_outlier_eff
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
  # Variant A: AR(1) innovations with stationary-prior c_1 contribution
  # Variant B: RW innovations; c_1 has diffuse prior — pass NA so
  #            update_vol's missing-data fallback leaves sigma_c[1]
  #            pulled toward its SV prior.
  rho_path <- if (variant == "B") 1 else state$rho
  nu_c <- c_path[-1L] - rho_path * c_path[-T_]
  nu_c_full <- if (variant == "B") {
    c(NA_real_, nu_c)
  } else {
    c(c_path[1L] * sqrt(max(1 - rho_path^2, 1e-6)), nu_c)
  }
  state$sigma_c <- update_vol(nu_c_full, state$sigma_c, state$gamma_c)

  nu_s <- rbind(s_path[1L, ],
                s_path[-1L, , drop = FALSE] - s_path[-T_, , drop = FALSE])
  for (i in seq_len(N)) {
    state$sigma_s[, i] <- update_vol(nu_s[, i], state$sigma_s[, i],
                                     state$gamma_s[i])
  }
  # eps_current is the scaled noise — with outliers, this equals
  # sigma_eps * s_outlier * unit_eps; without, sigma_eps * unit_eps.
  if (use_outliers) {
    obs_mask <- !is.na(y)

    # Step A: update per-obs s_outlier from whitened residual.
    # Under the model, w_t = eps_current_t / sigma_eps_t ~ N(0, s_t^2).
    state$s_outlier[!obs_mask] <- 1
    for (i in seq_len(N)) {
      mask_i <- obs_mask[, i]
      if (!any(mask_i)) next
      probs_i <- c(state$ps[i],
                   rep((1 - state$ps[i]) /
                         (length(config$s_vals) - 1L),
                       length(config$s_vals) - 1L))
      w_res_i <- eps_current[mask_i, i] / state$sigma_eps[mask_i, i]
      state$s_outlier[mask_i, i] <-
        update_scl(w_res_i, config$s_vals, probs_i)
    }

    # Step B: sigma_eps SV update on eps_current / s_outlier
    # (removes the outlier inflation, leaving sigma_eps * unit_eps).
    for (i in seq_len(N)) {
      state$sigma_eps[, i] <- update_vol(
        eps_current[, i] / state$s_outlier[, i],
        state$sigma_eps[, i], state$gamma_eps[i]
      )
    }

    # Step C: Beta-Bernoulli update of ps over OBSERVED cells.
    indicator <- matrix(NA, T_, N)
    indicator[obs_mask] <-
      (state$s_outlier[obs_mask] == config$s_vals[1L])
    state$ps <- update_ps(indicator, config$ps_a_prior, config$ps_b_prior)
  } else {
    for (i in seq_len(N)) {
      state$sigma_eps[, i] <- update_vol(eps_current[, i],
                                         state$sigma_eps[, i],
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
  # Variant A: conditional Gibbs (MH not yet wired for MA path).
  # Variant B: no rho update — leave at sentinel value of 1.
  if (variant == "B") {
    state$rho <- 1
  } else {
    state$rho <- .update_rho(
      c_path, state$sigma_c,
      prior_mean = config$rho_prior_mean,
      prior_sd   = config$rho_prior_sd
    )
  }

  # --- 7. Lambda update (constant, conjugate) ----------------------
  # Use the "denoised" obs (no MA contribution): obs - MA_contribution
  # = lambda_i * c_t + s_{i,t}. With current eps draws we can subtract
  # the MA piece from y to get a clean lambda regression. NA cells in
  # y carry through to y_denoised — `.update_lambda` skips NA rows
  # per-sector (added for the sub-group quarterly panel).
  ma_residual <- matrix(0, T_, N)
  for (i in seq_len(N)) {
    p <- meta$eps_start[i]
    for (l in 0:q_MA) {
      coef <- if (l == 0L) 1 else state$theta[i, l]
      ma_residual[, i] <- ma_residual[, i] + coef * as.numeric(states[p + l, ])
    }
  }
  y_denoised <- y - ma_residual
  state$lambda <- .update_lambda(
    y_denoised, c_path, s_path,
    sigma_eps = matrix(1, T_, N),     # noise removed; use unit sigma
    ref        = config$ref,
    prior_mean = config$lambda_prior_mean,
    prior_sd   = config$lambda_prior_sd
  )

  # --- 8. Theta update (per sector, conjugate Normal + invertibility)
  # With outliers ON, divide by sigma_eps * s_outlier (the effective
  # per-period noise scale) so the regression target is unit-variance.
  # Drop NA cells (staggered series starts) from the regression.
  theta_new <- matrix(0, N, q_MA)
  for (i in seq_len(N)) {
    p <- meta$eps_start[i]
    scale_t <- if (use_outliers) {
      state$sigma_eps[, i] * state$s_outlier[, i]
    } else {
      state$sigma_eps[, i]
    }
    u_t <- (y[, i] - state$lambda[i] * c_path - s_path[, i]) / scale_t
    x_th <- matrix(0, T_, q_MA)
    for (l in seq_len(q_MA)) {
      lag_eps <- numeric(T_)
      for (t in seq_len(T_)) {
        if (t - l >= 1L) lag_eps[t] <- as.numeric(states[p + l, t])
      }
      x_th[, l] <- lag_eps / scale_t
    }
    eligible <- (q_MA + 1L):T_
    keep <- eligible[!is.na(u_t[eligible])]
    if (length(keep) < q_MA + 1L) {
      theta_new[i, ] <- rnorm(q_MA, 0,
                              sqrt(1 / config$theta_prec_prior))
    } else {
      theta_new[i, ] <- update_theta(
        u_t[keep], x_th[keep, , drop = FALSE],
        prec_prior = config$theta_prec_prior
      )
    }
  }
  state$theta <- theta_new

  state$c <- c_path
  state$s <- s_path
  state
}
