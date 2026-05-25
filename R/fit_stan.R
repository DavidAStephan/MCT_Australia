# cmdstanr wrappers for the two MCT Stan models. Used by both the one-off
# runner script (scripts/run_real_fits.R) and the targets pipeline (Step 8).

#' Default deterministic init for `fit_mct()`. Starts every chain with all
#' z innovations at 0, small RW step sizes, and rho centred at 0.5 — keeps
#' the sampler in a well-conditioned region during warmup. The crucial bit
#' is that all chains start near the same point (not at random U(-2,2) on
#' the unconstrained scale), so they explore one mode and mix.
default_init <- function(stan_data, variant) {
  T_ <- stan_data$T; N <- stan_data$N
  if (variant == "C") {
    return(list(
      mu_hc = 0, mu_hec = 0,
      mu_hs = rep(0, N), mu_he = rep(0, N),
      sigma_hc = 0.02, sigma_hec = 0.02,
      sigma_hs = rep(0.02, N), sigma_he = rep(0.02, N),
      tau_c_init = 0,
      s_init = rep(0, N),
      lambda_tau_fixed = rep(1, N - 1),
      lambda_eps_fixed = rep(1, N - 1),
      z_hc = rep(0, T_),
      z_hec = rep(0, T_),
      z_hs = matrix(0, T_, N),
      z_he = matrix(0, T_, N),
      z_s  = matrix(0, T_, N),
      z_tau_c = rep(0, T_),
      z_c     = rep(0, T_)
    ))
  }
  if (variant %in% c("Ac", "Act")) {
    # Variants Ac and Act (constant-lambda A; Act = Ac + reduce_sum
    # within-chain parallelism). Same parameter shape — no z_lambda, no
    # sigma_lambda, no lambda_init; loadings are a time-invariant N-1
    # vector lambda_const.
    return(list(
      mu_hc = 0,
      mu_hs = rep(0, N),
      mu_he = rep(0, N),
      sigma_hc = 0.03,
      sigma_hs = rep(0.03, N),
      sigma_he = rep(0.03, N),
      rho = 0.5,
      s_init = rep(0, N),
      lambda_const = rep(1, N - 1),
      z_hc = rep(0, T_),
      z_hs = matrix(0, T_, N),
      z_he = matrix(0, T_, N),
      z_s  = matrix(0, T_, N),
      z_c  = rep(0, T_)
    ))
  }
  if (variant %in% c("Akf", "Akfs")) {
    # Variants Akf and Akfs (Phase 3: KF marginalization; Akfs = sparse
    # implementation of the same model). No z_c, no z_s — these are
    # marginalized by the Kalman filter. Same SV innovations as Ac.
    return(list(
      mu_hc = 0,
      mu_hs = rep(0, N),
      mu_he = rep(0, N),
      sigma_hc = 0.03,
      sigma_hs = rep(0.03, N),
      sigma_he = rep(0.03, N),
      rho = 0.5,
      s_init = rep(0, N),
      lambda_const = rep(1, N - 1),
      z_hc = rep(0, T_),
      z_hs = matrix(0, T_, N),
      z_he = matrix(0, T_, N)
    ))
  }
  init <- list(
    mu_hc = 0,
    mu_hs = rep(0, N),
    mu_he = rep(0, N),
    sigma_hc = 0.03,
    sigma_hs = rep(0.03, N),
    sigma_he = rep(0.03, N),
    sigma_lambda = rep(0.005, N - 1),
    s_init = rep(0, N),
    lambda_init = rep(1, N - 1),
    z_hc = rep(0, T_),
    z_hs = matrix(0, T_, N),
    z_he = matrix(0, T_, N),
    z_s  = matrix(0, T_, N),
    z_lambda = matrix(0, T_, N - 1),
    z_c  = rep(0, T_)
  )
  if (variant == "A") init$rho <- 0.5
  init
}

#' Fit a Variant A or B MCT model to a `stan_data` list from
#' `build_stan_data()`. The Stan-required fields are forwarded; the
#' diagnostic fields (`y`, `obs_type`, `dates`, `groups`) are stripped.
#'
#' @param stan_data List from `build_stan_data()`.
#' @param variant "A" (AR(1) common) or "B" (RW common).
#' @param model_dir Directory holding the .stan files.
#' @param chains, parallel_chains, iter_warmup, iter_sampling,
#'   adapt_delta, max_treedepth, seed Standard cmdstanr/NUTS arguments.
#'   `adapt_delta` defaulted up from the brief's 0.95 to 0.97 after the first
#'   real-data run showed mode-trapping; combined with the deterministic
#'   `default_init` this gave clean mixing on T=435 N=11.
#' @param init Either a list (same for every chain), a function taking
#'   `chain_id` and returning a list, or NULL to use Stan's random U(-2,2).
#'   Default `default_init(stan_data, variant)`.
#' @param save_dir If non-NULL, save the fit via `$save_object()` to
#'   `{save_dir}/fit_{variant}.rds`. Default `outputs/draws`.
#' @return A CmdStanMCMC object.
fit_mct <- function(stan_data,
                    variant = c("A", "B", "C", "Ac", "Act", "Akf", "Akfs"),
                    model_dir = "stan",
                    chains = 4,
                    parallel_chains = chains,
                    iter_warmup = 1500,
                    iter_sampling = 1500,
                    adapt_delta = 0.97,
                    max_treedepth = 12,
                    seed = 20260522,
                    init = NULL,
                    # Mass-matrix adaptation windowing (plan.md #5). NULL =
                    # use Stan defaults (75/50/25). Larger init+term buffers
                    # give the dense mass matrix more time to converge,
                    # which often pays for itself on stiff geometries like
                    # the SV+TVP state-space here.
                    init_buffer = NULL,
                    term_buffer = NULL,
                    window = NULL,
                    # Within-chain parallelism (plan.md #3). Only used by
                    # threaded variants ("Act"). NULL = 1 = no threading.
                    threads_per_chain = NULL,
                    save_dir = "outputs/draws") {
  variant <- match.arg(variant)
  model_path <- file.path(model_dir, paste0("mct_aus_", variant, ".stan"))
  stopifnot("Stan model file not found" = file.exists(model_path))
  threaded_variants <- c("Act")

  stan_input <- stan_data[c(
    "T", "N", "ref", "w",
    "n_obs_m", "t_m", "i_m", "y_m",
    "n_obs_q", "t_q", "i_q", "y_q"
  )]
  # cmdstanr accepts a function (called once per chain) or a list-of-lists.
  # Default: one shared deterministic init per chain.
  if (is.null(init)) {
    init_list <- default_init(stan_data, variant)
    init <- rep(list(init_list), chains)
  } else if (is.list(init) && !is.list(init[[1]])) {
    init <- rep(list(init), chains)
  }

  cpp_opts <- if (variant %in% threaded_variants) {
    list(stan_threads = TRUE)
  } else {
    NULL
  }
  m <- cmdstanr::cmdstan_model(model_path, cpp_options = cpp_opts)

  message(sprintf(
    "[fit_mct/%s] T=%d N=%d  n_obs_m=%d n_obs_q=%d  chains=%d  warmup=%d  sampling=%d  adapt_delta=%.2f",
    variant, stan_data$T, stan_data$N,
    stan_data$n_obs_m, stan_data$n_obs_q,
    chains, iter_warmup, iter_sampling, adapt_delta
  ))

  t0 <- Sys.time()
  sample_args <- list(
    data = stan_input,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    init = init,
    refresh = max(50, (iter_warmup + iter_sampling) %/% 20),
    show_messages = TRUE
  )
  if (!is.null(init_buffer))       sample_args$init_buffer       <- init_buffer
  if (!is.null(term_buffer))       sample_args$term_buffer       <- term_buffer
  if (!is.null(window))            sample_args$window            <- window
  if (!is.null(threads_per_chain)) sample_args$threads_per_chain <- threads_per_chain
  fit <- do.call(m$sample, sample_args)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("[fit_mct/%s] sampled in %.1f min", variant, elapsed))

  # Multi-modality watchdog. Variant C in particular has occasional warmup
  # failures where one chain gets stuck in a saddle-point and never explores
  # the rest of the posterior. ebfmi < 0.3 is the conventional cutoff for
  # "this chain is not mixing"; bad ebfmi shows up before bad R-hat does.
  diag <- fit$diagnostic_summary()
  bad <- which(diag$ebfmi < 0.3)
  if (length(bad) > 0) {
    message(sprintf(
      "[fit_mct/%s] WARNING: %d/%d chains have ebfmi < 0.3 (chain ids: %s; ebfmi values: %s). Inspect $diagnostic_summary() and consider discarding these chains before summarising.",
      variant, length(bad), length(diag$ebfmi),
      paste(bad, collapse = ","),
      paste(sprintf("%.3f", diag$ebfmi[bad]), collapse = ",")
    ))
  }

  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    out_path <- file.path(save_dir, paste0("fit_", variant, ".rds"))
    fit$save_object(file = out_path)
    message("[fit_mct/", variant, "] saved fit object to ", out_path)
  }

  fit
}

#' Export a slim parquet of the headline generated quantities for the
#' dashboard and downstream analysis. Drops the bulky per-month per-sector
#' state arrays (h_c, h_s, h_e, lambda, s_trend, z_*) and keeps only:
#'   - `trend[T]`
#'   - `common_share[T]`
#'   - `common_transitory[T]`
#'   - `sector_trend[T, N]` (re-exported s_trend)
#' Plus `log_lik` for the LOO comparison in Step 7.
export_fit_parquet <- function(fit, variant,
                               out_dir = "outputs/draws") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  d <- fit$draws(
    variables = c("trend", "common_share", "common_transitory",
                  "sector_trend", "log_lik"),
    format = "draws_df"
  )
  # draws_df carries posterior-specific attributes; strip them so arrow
  # serialises the columns as a plain table.
  d <- as.data.frame(d)
  out_path <- file.path(out_dir, paste0("fit_", variant, "_slim.parquet"))
  arrow::write_parquet(d, out_path)
  message(sprintf(
    "[export/%s] wrote %d draws x %d cols to %s",
    variant, nrow(d), ncol(d), out_path
  ))
  invisible(out_path)
}

#' Compose fetch + prep + Stan-data assembly into a single call. Convenient
#' for both the one-off runner and the targets pipeline.
build_real_stan_data <- function() {
  q <- fetch_cpi_quarterly()
  m <- fetch_cpi_monthly()
  w <- fetch_weights()
  infl <- demean_series(build_inflation_series(dplyr::bind_rows(q, m)))
  build_stan_data(infl, weights = w)
}
