# S3 wrapper around fit_mct_gibbs() output that mimics the cmdstanr
# CmdStanMCMC API (`fit$draws(...)`, `fit$diagnostic_summary()`,
# `fit$save_object(...)`). This lets all existing downstream code in
# R/postprocess.R, R/compare_variants.R, the dashboard etc. consume a
# Gibbs fit without modification.
#
# Step 8 of gibbs_port_plan.md.

#' Construct an mct_gibbs_fit object from a fit_mct_gibbs() result.
#'
#' @param gibbs_output Output of fit_mct_gibbs(): list with `draws`,
#'   `config`, `n_iter`, `elapsed_sec`.
#' @param stan_data List in the shape build_stan_data() returns
#'   (must include `dates`, `groups`, `w`).
#' @return An mct_gibbs_fit S3 object with closure-based methods.
mct_gibbs_fit <- function(gibbs_output, stan_data) {
  draws <- gibbs_output$draws
  config <- gibbs_output$config
  T_ <- config$T_
  N  <- config$N

  # Pre-compute derived (generated-quantities-equivalent) draws so we
  # only do the work once, not on every $draws() call.
  derived <- .compute_derived_gq(draws, stan_data)

  n_keep <- length(draws$rho)

  # Build the closure-based "method" list.
  obj <- list()

  obj$draws <- function(variables = NULL, format = "draws_matrix",
                       ...) {
    # Honour cmdstanr's main supported formats; for anything else we
    # convert from draws_matrix.
    if (is.null(variables)) {
      stop("mct_gibbs_fit$draws requires 'variables' to be specified.")
    }
    out <- .build_draws_matrix(variables, draws, derived, T_, N, n_keep)
    if (format %in% c("draws_matrix", "matrix")) {
      return(out)
    }
    if (format %in% c("draws_df", "df")) {
      return(posterior::as_draws_df(out))
    }
    if (format == "draws_array") {
      return(posterior::as_draws_array(out))
    }
    stop("Unsupported draws format for mct_gibbs_fit: ", format)
  }

  obj$diagnostic_summary <- function(...) {
    # Gibbs has no divergences and no treedepth. ebfmi is a HMC concept;
    # return 1.0 as a no-op so downstream watchdogs treat the chain as
    # healthy. n_chains = 1 because fit_mct_gibbs is single-chain.
    list(
      num_divergent     = 0L,
      num_max_treedepth = 0L,
      ebfmi             = 1.0
    )
  }

  obj$save_object <- function(file) {
    saveRDS(obj, file)
    invisible(file)
  }

  # Convenience accessors
  obj$n_iter      <- gibbs_output$n_iter
  obj$elapsed_sec <- gibbs_output$elapsed_sec
  obj$config      <- config
  obj$stan_data   <- stan_data

  class(obj) <- c("mct_gibbs_fit", "list")
  obj
}

# Compute derived quantities (trend, common_share, log_lik) for each
# draw. These match Stan's generated_quantities block exactly so the
# downstream postprocess code is sampler-agnostic.
.compute_derived_gq <- function(draws, stan_data) {
  T_ <- nrow(draws$c)
  N  <- dim(draws$s)[2L]
  n_keep <- ncol(draws$c)
  w  <- stan_data$w

  # trend[T_, n_keep]:
  #   Variant B (RW common trend — production model):
  #     trend = sum_i w_i * (lambda_i * c_t + s_{i,t})
  #           = (sum_i w_i * lambda_i) * c_t + sum_i w_i * s_{i,t}
  #     The common-factor c_t is part of the trend (Stock-Watson common
  #     drift), so we MUST include it.
  #   Variant A (AR(1) cycle):
  #     c_t is a transitory cycle; trend = sum_i w_i * s_{i,t} only.
  # We detect Variant B by all(rho == 1) — the gibbs_sweep pins rho at
  # 1 throughout for B.
  variant_B <- all(draws$rho == 1)
  trend <- matrix(NA_real_, T_, n_keep)
  for (k in seq_len(n_keep)) {
    sec <- as.numeric(draws$s[, , k] %*% w)
    if (variant_B) {
      wlam <- sum(w * draws$lambda[, k])           # scalar (lambda const)
      trend[, k] <- wlam * draws$c[, k] + sec
    } else {
      trend[, k] <- sec
    }
  }

  # common_share[t, k]:
  #   wlam = sum_i w_i * lambda_i (note: lambda is constant across t)
  #   com_var = wlam^2 * sigma_c[t]^2
  #   sec_var = sum_i w_i^2 * sigma_s[t, i]^2
  #   share = com_var / (com_var + sec_var)
  common_share <- matrix(NA_real_, T_, n_keep)
  for (k in seq_len(n_keep)) {
    wlam2 <- sum(w * draws$lambda[, k])^2
    com_var <- wlam2 * draws$sigma_c[, k]^2
    sec_var <- as.numeric((draws$sigma_s[, , k]^2) %*% (w^2))
    common_share[, k] <- com_var / (com_var + sec_var)
  }

  # log_lik per observation per draw.
  #
  # Monthly obs (i_m, t_m, y_m):
  #   mu = lambda[i] * c[t] + s[t, i]
  #   sigma = sigma_eps[t, i] * s_outlier[t, i]   (last factor = 1 if no outliers)
  #
  # Quarterly obs (i_q, t_q, y_q) — y is the avg of 3 monthly latents:
  #   mu = lambda[i] * (c[t-2] + c[t-1] + c[t])/3 + (s[t-2,i] + s[t-1,i] + s[t,i])/3
  #   var = (sigma_eff[t-2,i]^2 + sigma_eff[t-1,i]^2 + sigma_eff[t,i]^2) / 9
  # matching build_mct_ssm_mixed's observation row for ot == 2.
  use_outliers <- !is.null(draws$s_outlier)

  n_obs_m <- if (!is.null(stan_data$n_obs_m)) stan_data$n_obs_m else 0L
  log_lik_m <- matrix(NA_real_, n_obs_m, n_keep)
  if (n_obs_m > 0L) {
    for (k in seq_len(n_keep)) {
      mu  <- draws$lambda[stan_data$i_m, k] *
               draws$c[cbind(stan_data$t_m, k)] +
               draws$s[cbind(stan_data$t_m, stan_data$i_m, k)]
      sig <- draws$sigma_eps[cbind(stan_data$t_m, stan_data$i_m, k)]
      if (use_outliers) {
        sig <- sig * draws$s_outlier[cbind(stan_data$t_m,
                                           stan_data$i_m, k)]
      }
      log_lik_m[, k] <- dnorm(stan_data$y_m, mu, sig, log = TRUE)
    }
  }

  n_obs_q <- if (!is.null(stan_data$n_obs_q)) stan_data$n_obs_q else 0L
  log_lik_q <- matrix(NA_real_, n_obs_q, n_keep)
  if (n_obs_q > 0L) {
    t_q  <- stan_data$t_q
    i_q  <- stan_data$i_q
    y_q  <- stan_data$y_q
    # All quarterly obs are placed at t >= 3 (last month of quarter), but
    # guard in case of bad inputs.
    stopifnot(all(t_q >= 3L))
    for (k in seq_len(n_keep)) {
      c_t   <- draws$c[cbind(t_q,        k)]
      c_tm1 <- draws$c[cbind(t_q - 1L,   k)]
      c_tm2 <- draws$c[cbind(t_q - 2L,   k)]
      s_t   <- draws$s[cbind(t_q,        i_q, k)]
      s_tm1 <- draws$s[cbind(t_q - 1L,   i_q, k)]
      s_tm2 <- draws$s[cbind(t_q - 2L,   i_q, k)]
      lam_i <- draws$lambda[i_q, k]
      mu <- lam_i * (c_t + c_tm1 + c_tm2) / 3 +
              (s_t + s_tm1 + s_tm2) / 3

      e_t   <- draws$sigma_eps[cbind(t_q,        i_q, k)]
      e_tm1 <- draws$sigma_eps[cbind(t_q - 1L,   i_q, k)]
      e_tm2 <- draws$sigma_eps[cbind(t_q - 2L,   i_q, k)]
      if (use_outliers) {
        e_t   <- e_t   * draws$s_outlier[cbind(t_q,      i_q, k)]
        e_tm1 <- e_tm1 * draws$s_outlier[cbind(t_q - 1L, i_q, k)]
        e_tm2 <- e_tm2 * draws$s_outlier[cbind(t_q - 2L, i_q, k)]
      }
      var_q <- (e_t^2 + e_tm1^2 + e_tm2^2) / 9
      log_lik_q[, k] <- dnorm(y_q, mu, sqrt(var_q), log = TRUE)
    }
  }
  log_lik <- rbind(log_lik_m, log_lik_q)

  list(
    trend         = trend,
    common_share  = common_share,
    log_lik       = log_lik,
    sector_trend  = draws$s,            # alias: same as s
    common_transitory = draws$c         # alias: same as c (Variant A)
  )
}

# Assemble a posterior::draws_matrix for the requested variable name(s).
.build_draws_matrix <- function(variables, draws, derived, T_, N, n_keep) {
  if (length(variables) > 1L) {
    parts <- lapply(variables, function(v)
      .build_draws_matrix(v, draws, derived, T_, N, n_keep))
    return(do.call(posterior::bind_draws,
                   c(parts, list(along = "variable"))))
  }
  v <- variables

  # Scalars
  if (v %in% c("rho")) {
    m <- matrix(draws$rho, ncol = 1L,
                dimnames = list(NULL, v))
    return(posterior::as_draws_matrix(m))
  }

  # Length-N vectors: lambda, gamma_s, gamma_eps
  if (v %in% c("lambda", "gamma_s", "gamma_eps")) {
    src <- switch(v,
                  lambda     = draws$lambda,
                  gamma_s    = draws$gamma_s,
                  gamma_eps  = draws$gamma_eps)
    # src is N x n_keep; t(src) gives n_keep x N with names "v[i]"
    m <- t(src)
    colnames(m) <- paste0(v, "[", seq_len(N), "]")
    return(posterior::as_draws_matrix(m))
  }

  # Length-T vectors: trend, common_share, c, sigma_c, common_transitory
  if (v %in% c("trend", "common_share", "c", "sigma_c",
               "common_transitory")) {
    src <- switch(v,
                  trend             = derived$trend,
                  common_share      = derived$common_share,
                  c                 = draws$c,
                  common_transitory = derived$common_transitory,
                  sigma_c           = draws$sigma_c)
    m <- t(src)
    colnames(m) <- paste0(v, "[", seq_len(T_), "]")
    return(posterior::as_draws_matrix(m))
  }

  # T x N matrix-valued: s, sigma_s, sigma_eps, sector_trend
  if (v %in% c("s", "sigma_s", "sigma_eps", "sector_trend")) {
    src <- switch(v,
                  s            = draws$s,
                  sigma_s      = draws$sigma_s,
                  sigma_eps    = draws$sigma_eps,
                  sector_trend = derived$sector_trend)
    # src is T x N x n_keep; flatten to n_keep x (T*N) col-major over (T,N)
    m <- matrix(NA_real_, n_keep, T_ * N)
    for (k in seq_len(n_keep)) m[k, ] <- as.numeric(src[, , k])
    # Column names: "v[t,i]" in column-major order
    nms <- as.vector(outer(seq_len(T_), seq_len(N),
                           function(t, i) sprintf("%s[%d,%d]", v, t, i)))
    colnames(m) <- nms
    return(posterior::as_draws_matrix(m))
  }

  # log_lik: a (n_obs_m + n_obs_q) x n_keep matrix
  if (v == "log_lik") {
    src <- derived$log_lik   # n_obs x n_keep
    m <- t(src)
    colnames(m) <- paste0(v, "[", seq_len(ncol(m)), "]")
    return(posterior::as_draws_matrix(m))
  }

  stop("mct_gibbs_fit: unknown variable '", v, "'.")
}
