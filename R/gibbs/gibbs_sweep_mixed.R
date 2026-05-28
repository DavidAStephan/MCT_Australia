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

#' One Gibbs sweep for the mixed-frequency MCT model.
#'
#' Supports Variant A (AR(1) common factor — default) and Variant B
#' (RW common factor) via `config$variant`.
#'
#' @param y         T x N matrix; obs values where present, NA where
#'                  missing. Quarterly values appear at the last month
#'                  of each quarter, monthly values at every observed t.
#' @param obs_type  T x N integer matrix; 0 = missing, 1 = monthly,
#'                  2 = quarterly average (requires t >= 3).
#' @param state     Current state list (same shape as gibbs_sweep).
#' @param config    Config list (same shape as gibbs_sweep). Pass
#'                  `variant = "B"` to switch to the RW common-trend
#'                  model; default is "A".
#' @return Updated state list.
gibbs_sweep_mixed <- function(y, obs_type, state, config) {
  T_ <- config$T_
  N  <- config$N
  use_outliers <- isTRUE(config$use_outliers)
  variant      <- if (is.null(config$variant)) "A" else config$variant

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
    var_s_init = config$var_s_init,
    variant = variant
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
  # Common-factor innovation residuals:
  #   Variant A: c_t - rho * c_{t-1}, with c_1 scaled by sqrt(1 - rho^2)
  #              (matches the stationary AR(1) prior on c_1).
  #   Variant B: c_t - c_{t-1} (RW innovation). c_1 has a diffuse prior
  #              (Sigma_1[1,1] = 4), so its innovation isn't tied to
  #              sigma_c[1] — use sigma_c[1]'s prior draw uninformed by
  #              c_1: pass c_1 itself as the "residual" with the diffuse
  #              prior scaling already absorbed (this just contributes a
  #              constant to the SV update at t=1, which is fine).
  rho_path <- if (variant == "B") 1 else state$rho
  nu_c <- c_path[-1L] - rho_path * c_path[-T_]
  nu_c_full <- if (variant == "B") {
    # c_1 has a diffuse prior (Sigma_1[1,1] = 4) — it carries no info
    # about sigma_c[1]. Pass NA so update_vol's missing-data fallback
    # leaves sigma_c[1] pulled toward its SV prior.
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

  # Measurement-noise residuals: ONLY at monthly observations.
  # Quarterly + missing positions get NA so update_vol's missing-data
  # fallback kicks in (no update; propagate the RW state).
  c_mat   <- matrix(c_path, T_, N)
  lam_mat <- matrix(state$lambda, T_, N, byrow = TRUE)
  eps_mat <- y - lam_mat * c_mat - s_path
  eps_mat[obs_type != 1L] <- NA

  if (use_outliers) {
    # v1 design: only monthly obs participate in the scale mixture.
    # Quarterly obs and missing cells get s_outlier = 1 (no inflation)
    # so the SSM quarterly noise variance stays clean and the ps
    # Bernoulli count only sees the cells that actually evaluated the
    # update. Extending the mixture to quarterly obs requires a
    # 3-lag-aware update_scl variant and is deferred.

    # Step A: update per-obs scale s_outlier given monthly residuals +
    # current sigma_eps. Quarterly + missing cells: pin s_outlier = 1.
    state$s_outlier[obs_type != 1L] <- 1
    for (i in seq_len(N)) {
      mask_m <- (obs_type[, i] == 1L)
      if (!any(mask_m)) next
      probs_i <- c(state$ps[i],
                   rep((1 - state$ps[i]) /
                       (length(config$s_vals) - 1L),
                       length(config$s_vals) - 1L))
      w_res_i <- eps_mat[mask_m, i] / state$sigma_eps[mask_m, i]
      state$s_outlier[mask_m, i] <-
        update_scl(w_res_i, config$s_vals, probs_i)
    }

    # Step B: update sigma_eps using residuals DIVIDED BY scale
    # (so eps_t / s_t = sigma_t * unit_eps_t matches update_vol's model).
    for (i in seq_len(N)) {
      state$sigma_eps[, i] <- update_vol(
        eps_mat[, i] / state$s_outlier[, i],
        state$sigma_eps[, i], state$gamma_eps[i]
      )
    }

    # Step C: update ps via Beta-Bernoulli — only over monthly obs.
    # Set quarterly + missing cells to NA so update_ps's na.rm drops
    # them from the count.
    indicator <- matrix(NA, T_, N)
    mask_m <- (obs_type == 1L)
    indicator[mask_m] <- (state$s_outlier[mask_m] == config$s_vals[1L])
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
  # Variant A: update rho via marginal MH (unbiased) or conditional
  #   Gibbs (faster but FFBS-noise biased — see update_rho_mh.R).
  # Variant B: no rho update — the common-factor transition is fixed
  #   at 1 (RW). Leave state$rho at its sentinel value of 1.
  if (variant == "B") {
    state$rho <- 1
  } else if (isTRUE(config$use_marginal_mh_rho)) {
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
