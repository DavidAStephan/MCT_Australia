# LOO-CV based variant comparison + side-by-side viz wrappers.

#' Compute approximate leave-one-out cross-validation for a fit. Uses the
#' `log_lik` array from the model's generated quantities.
#'
#' @param fit A CmdStanMCMC.
#' @param r_eff_cores Cores for effective sample size estimation (default 1).
#' @return A `loo` object.
loo_for_fit <- function(fit, r_eff_cores = 1L) {
  ll <- fit$draws(variables = "log_lik", format = "matrix")
  # `loo::loo.array()` wants (iter, chain, n_obs) or `loo.matrix()` (draws, n_obs).
  # Use the matrix form — `fit$draws(format="matrix")` collapses chains.
  ll <- as.matrix(ll)
  # Drop columns where log_lik is exactly 0 across all draws (i.e., padded
  # cells that the model never evaluated). Shouldn't happen with the flat
  # observation arrays but guard anyway.
  keep <- colSums(abs(ll)) > 0
  ll <- ll[, keep, drop = FALSE]

  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(seq_len(fit$num_chains()),
                                                    each = fit$metadata()$iter_sampling),
                             cores = r_eff_cores)
  loo::loo(ll, r_eff = r_eff, cores = r_eff_cores)
}

#' Side-by-side LOO comparison across any number of fitted variants.
#'
#' @param loo_list Named list of `loo` objects from `loo_for_fit()`. Names
#'   become row labels in the returned `compare.loo` matrix. Must contain at
#'   least two entries.
#' @return A `compare.loo` matrix; the top row (`elpd_diff = 0`) is the
#'   best-fitting model, and `elpd_diff` for subsequent rows is the gap to
#'   the best model (negative — they fit worse).
compare_fits_loo <- function(loo_list, ...) {
  # Back-compat: accept compare_fits_loo(loo_A, loo_B) as well as a list.
  if (!is.list(loo_list) || inherits(loo_list, "loo")) {
    loo_list <- c(list(loo_list), list(...))
    if (is.null(names(loo_list)) || any(names(loo_list) == "")) {
      names(loo_list) <- LETTERS[seq_along(loo_list)]
    }
  }
  stopifnot("Need at least two loo objects" = length(loo_list) >= 2)
  loo::loo_compare(loo_list)
}

#' Human-readable one-paragraph write-up of the LOO comparison across any
#' number of variants. Returns a character scalar (markdown).
#'
#' Caveats encoded:
#'   - flags high Pareto-k counts per variant
#'   - reports the elpd gap between the best and the runner-up plus its SE
#'   - states whether that gap is within 2 SE.
#'
#' @param loo_list Named list of `loo` objects. Variant descriptions for the
#'   trailing sentence are taken from `variant_descriptions`; unrecognised
#'   names fall back to a generic blurb.
variant_comparison_text <- function(loo_list, ...,
                                    variant_descriptions = .default_variant_descriptions()) {
  # Back-compat: accept variant_comparison_text(loo_A, loo_B) as well as a list.
  if (!is.list(loo_list) || inherits(loo_list, "loo")) {
    loo_list <- c(list(loo_list), list(...))
    if (is.null(names(loo_list)) || any(names(loo_list) == "")) {
      names(loo_list) <- LETTERS[seq_along(loo_list)]
    }
  }
  cmp <- compare_fits_loo(loo_list)
  preferred <- rownames(cmp)[1]
  runner_up <- rownames(cmp)[2]
  elpd_diff <- cmp[2, "elpd_diff"]
  se_diff <- cmp[2, "se_diff"]
  ratio <- abs(elpd_diff) / se_diff

  judgement <- if (ratio < 2) {
    sprintf("The gap between Variant %s and runner-up Variant %s (%.1f elpd, SE %.1f) is within 2 SE — the LOO ranking is not decisive.",
            preferred, runner_up, abs(elpd_diff), se_diff)
  } else {
    sprintf("Variant %s is preferred over runner-up Variant %s by **%.1f elpd points (SE %.1f, |Δ|/SE = %.1f).**",
            preferred, runner_up, abs(elpd_diff), se_diff, ratio)
  }

  high_k_counts <- vapply(loo_list, function(x) sum(x$diagnostics$pareto_k > 0.7),
                          integer(1))
  high_k_note <- if (sum(high_k_counts) > 0) {
    parts <- mapply(function(n, k) sprintf("%d (Variant %s)", k, n),
                    names(high_k_counts), high_k_counts,
                    USE.NAMES = FALSE)
    sprintf(" PSIS-LOO Pareto-k diagnostics flagged %s observations with k > 0.7; treat the elpd comparison as indicative rather than decisive.",
            paste(parts, collapse = ", "))
  } else {
    ""
  }

  descs <- vapply(names(loo_list), function(n) {
    if (!is.null(variant_descriptions[[n]])) variant_descriptions[[n]]
    else sprintf("Variant %s", n)
  }, character(1))

  paste0(
    "**LOO comparison.** ", judgement, high_k_note,
    " The decomposition each variant offers is also distinct: ",
    paste(descs, collapse = "; "), "."
  )
}

.default_variant_descriptions <- function() {
  list(
    A = "Variant A treats the common factor as a stationary AR(1) cycle, so its trend is the weighted sum of sector-specific trends",
    B = "Variant B treats the common factor as a random-walk trend (the canonical Stock-Watson decomposition)",
    C = "Variant C combines both — a random-walk common trend AND an i.i.d. common cycle, each with stochastic volatility, with time-invariant loadings (Gaussian noise; MA(3) and Student-t innovations from the NY Fed spec were deferred — see README)"
  )
}

#' Headline-trend comparison plot. Accepts any number of fits as `...` and a
#' matching `labels` vector. Wraps the trend_path summarisation + viz.
#'
#' @examples
#' compare_trend_plot(fit_A, fit_B, dates = stan_data$dates)
#' compare_trend_plot(fit_A, fit_B, fit_C, dates = stan_data$dates,
#'                    labels = c("A: AR(1)", "B: RW trend", "C: combined"))
compare_trend_plot <- function(..., dates,
                               labels = NULL,
                               probs = c(0.16, 0.5, 0.84),
                               title = "Trend inflation by variant") {
  fits <- list(...)
  if (is.null(labels)) {
    labels <- if (length(fits) == 2) c("A: AR(1) common", "B: RW common trend")
              else if (length(fits) == 3) c("A: AR(1) common", "B: RW common trend",
                                            "C: combined trend + cycle")
              else LETTERS[seq_along(fits)]
  }
  summaries <- lapply(fits, trend_path, dates = dates, probs = probs)
  plot_compare_paths(summaries, labels = labels,
                     title = title,
                     subtitle = "Median + 68% credible band")
}

#' Common-share comparison plot. Same dispatch as compare_trend_plot.
#'
#' `since` mirrors `plot_common_share()` — common_share is poorly
#' identified on quarterly-only data, so we typically restrict to the
#' monthly-data era (Apr 2024 onward).
compare_common_share_plot <- function(..., dates,
                                      labels = NULL,
                                      probs = c(0.16, 0.5, 0.84),
                                      title = "Common share of inflation variance",
                                      since = NULL) {
  fits <- list(...)
  if (is.null(labels)) {
    labels <- if (length(fits) == 2) c("A: AR(1) common", "B: RW common trend")
              else if (length(fits) == 3) c("A: AR(1) common", "B: RW common trend",
                                            "C: combined trend + cycle")
              else LETTERS[seq_along(fits)]
  }
  summaries <- lapply(fits, common_share_path, dates = dates, probs = probs)
  if (!is.null(since)) {
    summaries <- lapply(summaries, function(s)
      s[s$date >= as.Date(since), , drop = FALSE])
  }
  plot_compare_paths(summaries, labels = labels,
                     title = title, ylab = NULL,
                     subtitle = "Median + 68% credible band") +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                labels = scales::percent_format(accuracy = 1))
}
