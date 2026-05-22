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

#' Side-by-side LOO comparison of variants A and B.
#'
#' @param loo_A,loo_B `loo` objects from `loo_for_fit()`.
#' @return A `compare.loo` matrix; positive elpd_diff in row 1 means that
#'   row's model is preferred.
compare_fits_loo <- function(loo_A, loo_B) {
  loo::loo_compare(list(A = loo_A, B = loo_B))
}

#' Human-readable one-paragraph write-up of the variant comparison, ready
#' to drop into the dashboard. Returns a character scalar (markdown).
#'
#' Caveats encoded:
#'   - flags high Pareto-k counts as a reason to distrust the comparison
#'   - reports elpd_diff and its SE
#'   - states which variant is preferred (or that the diff is within 2 SE).
variant_comparison_text <- function(loo_A, loo_B) {
  cmp <- compare_fits_loo(loo_A, loo_B)
  preferred <- rownames(cmp)[1]
  elpd_diff <- cmp[2, "elpd_diff"]
  se_diff <- cmp[2, "se_diff"]
  n_high_k_A <- sum(loo_A$diagnostics$pareto_k > 0.7)
  n_high_k_B <- sum(loo_B$diagnostics$pareto_k > 0.7)
  ratio <- elpd_diff / se_diff
  judgement <- if (abs(ratio) < 2) {
    sprintf("The difference (%.1f, SE %.1f) is within 2 SE — the LOO ranking is not decisive.",
            elpd_diff, se_diff)
  } else {
    sprintf("Variant %s is preferred by **%.1f elpd points (SE %.1f, |Δ|/SE = %.1f).**",
            preferred, abs(elpd_diff), se_diff, abs(ratio))
  }
  high_k_note <- if (n_high_k_A + n_high_k_B > 0) {
    sprintf(" PSIS-LOO Pareto-k diagnostics flagged %d (Variant A) and %d (Variant B) observations with k > 0.7; the elpd comparison should be treated as indicative rather than decisive.",
            n_high_k_A, n_high_k_B)
  } else {
    ""
  }
  paste0(
    "**LOO comparison.** ", judgement, high_k_note,
    " The decomposition each variant offers is also distinct: Variant A treats the common factor as transitory, so its trend is the weighted sum of sector-specific trends; Variant B treats the common factor as the drifting common trend itself, the canonical Stock-Watson decomposition."
  )
}

#' Headline-trend comparison plot wrapper.
compare_trend_plot <- function(fit_A, fit_B, dates,
                               probs = c(0.16, 0.5, 0.84),
                               title = "Trend inflation: Variant A vs Variant B") {
  sA <- trend_path(fit_A, dates, probs)
  sB <- trend_path(fit_B, dates, probs)
  plot_compare_paths(sA, sB, labels = c("A: AR(1) common", "B: RW common trend"),
                     title = title,
                     subtitle = "Median + 68% credible band")
}

#' Common-share comparison plot wrapper.
compare_common_share_plot <- function(fit_A, fit_B, dates,
                                      probs = c(0.16, 0.5, 0.84),
                                      title = "Common share of inflation variance") {
  sA <- common_share_path(fit_A, dates, probs)
  sB <- common_share_path(fit_B, dates, probs)
  plot_compare_paths(sA, sB, labels = c("A: AR(1) common", "B: RW common trend"),
                     title = title, ylab = NULL,
                     subtitle = "Median + 68% credible band") +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                labels = scales::percent_format(accuracy = 1))
}
