# Main entry point for the Gibbs port of the MCT Variant A estimator.
# Mirrors fit_mct() in API shape so downstream code can opt into the
# Gibbs sampler with minimal changes.
#
# Per gibbs_port_plan.md Step 7 — first usable end-to-end Gibbs fit.
# Step 9 will add mixed-frequency observations and AR(1)-specific
# refinements; Step 10 will validate against the Stan Variant Ac fit
# on real ABS data.

#' Fit the MCT Variant A model via the Gibbs sampler.
#'
#' @param y         T x N matrix of inflations (monthly only for now;
#'                  mixed-freq added in Step 9).
#' @param ref       Reference sector (1-based index); its loading is
#'                  pinned at 1 for identification.
#' @param n_burn    Number of burn-in iterations.
#' @param n_draw   Number of post-burn draws to collect.
#' @param n_thin    Keep every n_thin-th draw post-burn (default 1).
#' @param init      Optional initial-state list. If NULL, sensible
#'                  defaults are used.
#' @param config    Optional config overrides; see defaults below.
#' @param verbose   If TRUE, print progress every 100 iters.
#' @return List with `draws` (collected samples) and `config`.
fit_mct_gibbs <- function(y,
                          obs_type = NULL,
                          ref = 1L,
                          n_burn = 1000L,
                          n_draw = 2000L,
                          n_thin = 1L,
                          init = NULL,
                          config = NULL,
                          verbose = TRUE) {
  T_ <- nrow(y)
  N  <- ncol(y)

  # Dispatch: if obs_type is given (Step 9 mixed-freq path), use
  # gibbs_sweep_mixed + augmented SSM. Otherwise (Step 7 monthly-only
  # path), use the simpler gibbs_sweep + N+1-dim SSM.
  mixed_freq <- !is.null(obs_type)
  if (mixed_freq) {
    stopifnot(all(dim(obs_type) == c(T_, N)),
              all(obs_type %in% c(0L, 1L, 2L)))
  }

  # Default config — closely mirrors plan.md/Stan Variant Ac priors
  cfg <- list(
    T_ = T_, N = N, ref = ref,
    nu_gam_prior      = 60,
    s2_gam_prior      = 1 / (60 * 12),
    rho_prior_mean    = 0.5,
    rho_prior_sd      = 0.3,
    lambda_prior_mean = 1,
    lambda_prior_sd   = 0.3,
    # var_s_init matches Stan Variant Ac prior s_init ~ N(0, 0.5^2).
    # Larger values (e.g. 100) make s_i and c trade signal in early
    # iterations — c gets absorbed into s, biasing rho downward.
    var_s_init        = 0.25,
    # rho update: marginal MH (recommended) avoids the FFBS-noise
    # bias of the conditional Gibbs path. Set FALSE to recover the
    # biased-but-faster conditional-Gibbs behaviour for A/B testing.
    use_marginal_mh_rho = TRUE,
    rho_prop_sd         = 0.05
  )
  if (!is.null(config)) {
    for (k in names(config)) cfg[[k]] <- config[[k]]
  }

  # Sensible initial state. Sigmas start at 0.5 (close to typical SV
  # scale for demeaned monthly inflation); gammas at 0.05; rho at 0.5;
  # lambda at 1.
  if (is.null(init)) {
    state <- list(
      c         = rep(0, T_),
      s         = matrix(0, T_, N),
      sigma_c   = rep(0.5, T_),
      sigma_s   = matrix(0.5, T_, N),
      sigma_eps = matrix(0.5, T_, N),
      gamma_c   = 0.05,
      gamma_s   = rep(0.05, N),
      gamma_eps = rep(0.05, N),
      rho       = 0.5,
      lambda    = rep(1, N)
    )
  } else {
    state <- init
  }
  state$lambda[ref] <- 1  # enforce identification

  # Pre-allocate posterior storage. Keep every n_thin-th post-burn draw.
  n_keep <- n_draw %/% n_thin
  draws <- list(
    c         = matrix(NA_real_, T_, n_keep),
    s         = array(NA_real_, c(T_, N, n_keep)),
    sigma_c   = matrix(NA_real_, T_, n_keep),
    sigma_s   = array(NA_real_, c(T_, N, n_keep)),
    sigma_eps = array(NA_real_, c(T_, N, n_keep)),
    gamma_c   = numeric(n_keep),
    gamma_s   = matrix(NA_real_, N, n_keep),
    gamma_eps = matrix(NA_real_, N, n_keep),
    rho       = numeric(n_keep),
    lambda    = matrix(NA_real_, N, n_keep)
  )

  # Track rho-MH acceptance rate (post-burn-in only)
  rho_n_accept <- 0L
  rho_n_total  <- 0L

  t0 <- Sys.time()
  total_iter <- n_burn + n_draw
  for (i in seq_len(total_iter)) {
    state <- if (mixed_freq) {
      gibbs_sweep_mixed(y, obs_type, state, cfg)
    } else {
      gibbs_sweep(y, state, cfg)
    }

    # Track rho-MH acceptance (post-burn-in only, mixed-freq path only)
    if (mixed_freq && isTRUE(cfg$use_marginal_mh_rho) && i > n_burn) {
      rho_n_total  <- rho_n_total + 1L
      if (isTRUE(state$last_rho_accept)) rho_n_accept <- rho_n_accept + 1L
    }

    # Collect draws after burn-in, every n_thin-th
    if (i > n_burn) {
      i_post <- i - n_burn
      if (i_post %% n_thin == 0L) {
        k <- i_post %/% n_thin
        draws$c[, k]          <- state$c
        draws$s[, , k]        <- state$s
        draws$sigma_c[, k]    <- state$sigma_c
        draws$sigma_s[, , k]  <- state$sigma_s
        draws$sigma_eps[, , k] <- state$sigma_eps
        draws$gamma_c[k]      <- state$gamma_c
        draws$gamma_s[, k]    <- state$gamma_s
        draws$gamma_eps[, k]  <- state$gamma_eps
        draws$rho[k]          <- state$rho
        draws$lambda[, k]     <- state$lambda
      }
    }

    if (verbose && (i %% 100L == 0L || i == total_iter)) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta <- elapsed / i * (total_iter - i)
      message(sprintf(
        "[fit_mct_gibbs] iter %d/%d (%.1f%%)  elapsed=%.1fs  eta=%.1fs",
        i, total_iter, 100 * i / total_iter, elapsed, eta
      ))
    }
  }

  out <- list(draws = draws, config = cfg, n_iter = total_iter,
              elapsed_sec = as.numeric(difftime(Sys.time(),
                                                t0, units = "secs")))
  if (rho_n_total > 0L) {
    out$rho_accept_rate <- rho_n_accept / rho_n_total
    if (verbose) {
      message(sprintf("[fit_mct_gibbs] rho-MH acceptance: %.1f%% (%d / %d post-burn)",
                      100 * out$rho_accept_rate, rho_n_accept, rho_n_total))
    }
  }
  out
}
