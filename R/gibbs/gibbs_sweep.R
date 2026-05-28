# One full Gibbs sweep of the MCT Variant A model. Updates all
# parameters in order:
#
#   1. (c, s) joint draw via simulation_smoother (FFBS)
#   2. SV processes sigma_c, sigma_s, sigma_eps via Kim-Shephard
#   3. RW step sizes gamma_c, gamma_s, gamma_eps via update_gam
#   4. AR(1) coef rho via truncated-normal conjugate update
#   5. Constant loadings lambda via Gaussian conjugate update
#
# Each call returns the full updated state, ready for the next sweep
# or for collection into the posterior draws.

#' Update rho (AR(1) coef on common factor) from conjugate truncated
#' Normal posterior given c and sigma_c.
#'
#' Model: c_t = rho * c_{t-1} + nu_t, nu_t ~ N(0, sigma_c[t]^2)
#' Prior: rho ~ N(prior_mean, prior_sd^2) truncated to (-1, 1)
.update_rho <- function(c, sigma_c, prior_mean = 0.5, prior_sd = 0.3) {
  T_ <- length(c)
  # Regress c[2:T] on c[1:(T-1)] with weights 1/sigma_c[2:T]^2
  y <- c[2:T_] / sigma_c[2:T_]
  x <- c[1:(T_ - 1L)] / sigma_c[2:T_]
  # Posterior precision = prior precision + sum(x^2)
  prec_post <- 1 / prior_sd^2 + sum(x^2)
  mean_post <- (prior_mean / prior_sd^2 + sum(x * y)) / prec_post
  sd_post   <- sqrt(1 / prec_post)
  # Sample truncated normal — accept/reject from N(mean_post, sd_post^2)
  # within (-1, 1). For our prior, posterior is usually well inside.
  repeat {
    draw <- rnorm(1, mean_post, sd_post)
    if (draw > -1 && draw < 1) return(draw)
  }
}

#' Update constant loadings lambda. For sectors i != ref, regress
#' (y[i,t] - s[i,t]) on c[t] with weights 1/sigma_eps[i,t]. Reference
#' sector pinned at 1. NA cells in y[, i] are skipped (handles staggered
#' series starts in panels where not every sector is observed throughout).
.update_lambda <- function(y, c, s, sigma_eps, ref,
                           prior_mean = 1, prior_sd = 0.3) {
  N      <- ncol(y)
  lambda <- numeric(N)
  lambda[ref] <- 1
  for (i in seq_len(N)) {
    if (i == ref) next
    obs <- !is.na(y[, i])
    if (!any(obs)) {
      lambda[i] <- rnorm(1, prior_mean, prior_sd)
      next
    }
    wy <- (y[obs, i] - s[obs, i]) / sigma_eps[obs, i]
    wx <- c[obs]                / sigma_eps[obs, i]
    prec_post <- 1 / prior_sd^2 + sum(wx^2)
    mean_post <- (prior_mean / prior_sd^2 + sum(wx * wy)) / prec_post
    lambda[i] <- rnorm(1, mean_post, sqrt(1 / prec_post))
  }
  lambda
}

#' One Gibbs sweep. Returns updated `state` list.
#'
#' Supports Variant A (AR(1) common factor — default for back-compat)
#' and Variant B (RW common factor) via `config$variant`. NA cells in
#' `y` are handled (simulation smoother integrates over missing obs;
#' lambda update skips them; sigma_eps RW propagates without info).
#'
#' @param y         T x N matrix of observed inflations. NA cells are
#'                  treated as missing observations.
#' @param state     List with current draws of c, s, sigma_c, sigma_s,
#'                  sigma_eps, gamma_c, gamma_s, gamma_eps, rho, lambda.
#' @param config    List of fixed config: N, T_, ref, nu_gam_prior,
#'                  s2_gam_prior, rho_prior_mean, rho_prior_sd,
#'                  lambda_prior_mean, lambda_prior_sd, var_s_init,
#'                  variant (default "A").
#' @return Updated `state` list.
gibbs_sweep <- function(y, state, config) {
  T_ <- config$T_
  N  <- config$N
  variant      <- if (is.null(config$variant)) "A" else config$variant
  use_outliers <- isTRUE(config$use_outliers)

  # --- 1. Build SSM from current params -----------------------------
  # With outliers ON, the EFFECTIVE measurement-noise SD seen by the
  # Kalman filter is sigma_eps * s_outlier (per-obs scale inflation).
  effective_sigma_eps <- if (use_outliers) {
    state$sigma_eps * state$s_outlier
  } else {
    state$sigma_eps
  }
  ssm <- build_mct_ssm(
    rho = state$rho, lambda = state$lambda,
    sigma_c = state$sigma_c, sigma_s = state$sigma_s,
    sigma_eps = effective_sigma_eps, T_ = T_, N = N,
    var_s_init = config$var_s_init,
    variant = variant
  )

  # --- 2. Joint draw of (c, s) via simulation_smoother --------------
  # simulate_ssm/kalman_filter expect Y as N x T_ matrix.
  # NA cells in y are missing observations — kalman_filter handles those.
  Y <- t(y)
  smooth <- simulation_smoother(Y, ssm)
  states <- smooth$states                     # (N+1) x T_
  c_path <- as.numeric(states[1L, ])
  s_path <- t(states[-1L, , drop = FALSE])    # T_ x N

  # --- 3. Update SV processes via Kim-Shephard ----------------------
  # Common factor innovations:
  #   Variant A: nu_c[t] = c[t] - rho * c[t-1], with c[1] scaled by
  #              sqrt(1 - rho^2) to use the stationary AR(1) prior.
  #   Variant B: nu_c[t] = c[t] - c[t-1] (RW). c_1 has a diffuse prior
  #              (Sigma_1[1,1] = 4) so it carries no info on sigma_c[1];
  #              pass NA there and let update_vol's missing-data
  #              fallback leave sigma_c[1] pulled toward its SV prior.
  rho_path <- if (variant == "B") 1 else state$rho
  nu_c <- c_path[-1L] - rho_path * c_path[-T_]
  nu_c_full <- if (variant == "B") {
    c(NA_real_, nu_c)
  } else {
    c(c_path[1L] * sqrt(max(1 - rho_path^2, 1e-6)), nu_c)
  }
  state$sigma_c <- update_vol(nu_c_full, state$sigma_c, state$gamma_c)

  # Sector trend RW innovations nu_s[i, t] = s[i, t] - s[i, t-1]
  nu_s <- rbind(s_path[1L, ],
                s_path[-1L, , drop = FALSE] - s_path[-T_, , drop = FALSE])
  for (i in seq_len(N)) {
    state$sigma_s[, i] <- update_vol(nu_s[, i], state$sigma_s[, i],
                                     state$gamma_s[i])
  }

  # Measurement noise eps[i, t] = y[i, t] - lambda[i] c[t] - s[i, t].
  # NA cells in y propagate to NA in eps_mat; update_vol's missing-data
  # fallback leaves sigma_eps unchanged there.
  c_mat   <- matrix(c_path, T_, N)
  lam_mat <- matrix(state$lambda, T_, N, byrow = TRUE)
  eps_mat <- y - lam_mat * c_mat - s_path

  if (use_outliers) {
    obs_mask <- !is.na(y)

    # Step A: update per-obs s_outlier from whitened residuals.
    # Missing cells get s_outlier = 1 (no inflation, irrelevant — the
    # SSM doesn't read those entries anyway since y is NA there).
    state$s_outlier[!obs_mask] <- 1
    for (i in seq_len(N)) {
      mask_i <- obs_mask[, i]
      if (!any(mask_i)) next
      probs_i <- c(state$ps[i],
                   rep((1 - state$ps[i]) /
                         (length(config$s_vals) - 1L),
                       length(config$s_vals) - 1L))
      w_res_i <- eps_mat[mask_i, i] / state$sigma_eps[mask_i, i]
      state$s_outlier[mask_i, i] <-
        update_scl(w_res_i, config$s_vals, probs_i)
    }

    # Step B: update sigma_eps using residuals divided by scale
    # (so eps_t / s_t = sigma_t * unit_eps_t matches update_vol's model).
    for (i in seq_len(N)) {
      state$sigma_eps[, i] <- update_vol(
        eps_mat[, i] / state$s_outlier[, i],
        state$sigma_eps[, i], state$gamma_eps[i]
      )
    }

    # Step C: update ps (probability of normal obs) via Beta-Bernoulli,
    # counting only OBSERVED cells. NA elsewhere drops them from the
    # count via update_ps's na.rm.
    indicator <- matrix(NA, T_, N)
    indicator[obs_mask] <-
      (state$s_outlier[obs_mask] == config$s_vals[1L])
    state$ps <- update_ps(indicator, config$ps_a_prior, config$ps_b_prior)
  } else {
    for (i in seq_len(N)) {
      state$sigma_eps[, i] <- update_vol(eps_mat[, i], state$sigma_eps[, i],
                                         state$gamma_eps[i])
    }
  }

  # --- 4. Update gamma step sizes via update_gam --------------------
  # gamma update target: 2 * diff(log(sigma)) is approximately N(0, gamma^2)
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

  # --- 5. Update rho via truncated-normal conjugate -----------------
  # Variant B: no rho update (transition is fixed at 1 / RW). Leave
  # state$rho at its sentinel value of 1.
  if (variant == "B") {
    state$rho <- 1
  } else {
    state$rho <- .update_rho(
      c_path, state$sigma_c,
      prior_mean = config$rho_prior_mean,
      prior_sd   = config$rho_prior_sd
    )
  }

  # --- 6. Update constant loadings lambda ---------------------------
  state$lambda <- .update_lambda(
    y, c_path, s_path, state$sigma_eps,
    ref        = config$ref,
    prior_mean = config$lambda_prior_mean,
    prior_sd   = config$lambda_prior_sd
  )

  # Stash the latest state paths for output collection
  state$c <- c_path
  state$s <- s_path
  state
}
