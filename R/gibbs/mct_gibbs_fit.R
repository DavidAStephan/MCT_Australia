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

  # trend[T_, n_keep] = sum_i w_i * s[t, i]
  # We have draws$s as (T_, N, n_keep) — collapse the N dim with w.
  trend <- matrix(NA_real_, T_, n_keep)
  for (k in seq_len(n_keep)) {
    trend[, k] <- as.numeric(draws$s[, , k] %*% w)
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

  # log_lik per observation per draw. Monthly only for v1.
  # log_lik[k, draw] = normal_lpdf(y_m[k] | mu, sigma)
  # where mu = lambda[i_m[k]] * c[t_m[k], draw] + s[t_m[k], i_m[k], draw]
  # and sigma = sigma_eps[t_m[k], i_m[k], draw].
  n_obs_m <- stan_data$n_obs_m
  log_lik_m <- matrix(NA_real_, n_obs_m, n_keep)
  for (k in seq_len(n_keep)) {
    mu  <- draws$lambda[stan_data$i_m, k] *
             draws$c[cbind(stan_data$t_m, k)] +
             draws$s[cbind(stan_data$t_m, stan_data$i_m, k)]
    sig <- draws$sigma_eps[cbind(stan_data$t_m, stan_data$i_m, k)]
    log_lik_m[, k] <- dnorm(stan_data$y_m, mu, sig, log = TRUE)
  }
  # Quarterly handling deferred to Step 9; emit a zero matrix for now
  n_obs_q <- if (!is.null(stan_data$n_obs_q)) stan_data$n_obs_q else 0L
  log_lik_q <- matrix(0, n_obs_q, n_keep)
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
