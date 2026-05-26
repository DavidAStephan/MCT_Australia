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
#' sector pinned at 1.
.update_lambda <- function(y, c, s, sigma_eps, ref,
                           prior_mean = 1, prior_sd = 0.3) {
  N      <- ncol(y)
  lambda <- numeric(N)
  lambda[ref] <- 1
  for (i in seq_len(N)) {
    if (i == ref) next
    # Regression: (y[, i] - s[, i]) = lambda[i] * c + noise(sd = sigma_eps[, i])
    wy <- (y[, i] - s[, i]) / sigma_eps[, i]
    wx <- c / sigma_eps[, i]
    prec_post <- 1 / prior_sd^2 + sum(wx^2)
    mean_post <- (prior_mean / prior_sd^2 + sum(wx * wy)) / prec_post
    lambda[i] <- rnorm(1, mean_post, sqrt(1 / prec_post))
  }
  lambda
}

#' One Gibbs sweep. Returns updated `state` list.
#'
#' @param y         T x N matrix of observed inflations (monthly only).
#' @param state     List with current draws of c, s, sigma_c, sigma_s,
#'                  sigma_eps, gamma_c, gamma_s, gamma_eps, rho, lambda.
#' @param config    List of fixed config: N, T_, ref, nu_gam_prior,
#'                  s2_gam_prior, rho_prior_mean, rho_prior_sd,
#'                  lambda_prior_mean, lambda_prior_sd, var_s_init.
#' @return Updated `state` list.
gibbs_sweep <- function(y, state, config) {
  T_ <- config$T_
  N  <- config$N

  # --- 1. Build SSM from current params -----------------------------
  ssm <- build_mct_ssm(
    rho = state$rho, lambda = state$lambda,
    sigma_c = state$sigma_c, sigma_s = state$sigma_s,
    sigma_eps = state$sigma_eps, T_ = T_, N = N,
    var_s_init = config$var_s_init
  )

  # --- 2. Joint draw of (c, s) via simulation_smoother --------------
  # simulate_ssm/kalman_filter expect Y as N x T_ matrix
  Y <- t(y)
  smooth <- simulation_smoother(Y, ssm)
  states <- smooth$states                     # (N+1) x T_
  c_path <- as.numeric(states[1L, ])
  s_path <- t(states[-1L, , drop = FALSE])    # T_ x N

  # --- 3. Update SV processes via Kim-Shephard ----------------------
  # Common factor: AR(1) innovations nu_c[t] = c[t] - rho * c[t-1]
  nu_c <- c_path[-1L] - state$rho * c_path[-T_]
  # update_vol expects a length-T input; pad nu_c[1] (no t=0 lag) with
  # the marginal c[1] under stationary distribution
  nu_c_full <- c(c_path[1L] * sqrt(1 - state$rho^2), nu_c)
  state$sigma_c <- update_vol(nu_c_full, state$sigma_c, state$gamma_c)

  # Sector trend RW innovations nu_s[i, t] = s[i, t] - s[i, t-1]
  nu_s <- rbind(s_path[1L, ],
                s_path[-1L, , drop = FALSE] - s_path[-T_, , drop = FALSE])
  for (i in seq_len(N)) {
    state$sigma_s[, i] <- update_vol(nu_s[, i], state$sigma_s[, i],
                                     state$gamma_s[i])
  }

  # Measurement noise eps[i, t] = y[i, t] - lambda[i] c[t] - s[i, t]
  c_mat   <- matrix(c_path, T_, N)
  lam_mat <- matrix(state$lambda, T_, N, byrow = TRUE)
  eps_mat <- y - lam_mat * c_mat - s_path
  for (i in seq_len(N)) {
    state$sigma_eps[, i] <- update_vol(eps_mat[, i], state$sigma_eps[, i],
                                       state$gamma_eps[i])
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
  state$rho <- .update_rho(
    c_path, state$sigma_c,
    prior_mean = config$rho_prior_mean,
    prior_sd   = config$rho_prior_sd
  )

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
