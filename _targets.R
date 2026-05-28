# targets pipeline for MCT_Australia. Run with:
#   targets::tar_make()
#
# Production model (as of 2026-05-28):
#   - 33 ABS CPI sub-groups (closer NY-Fed-like granularity)
#   - Quarterly latent clock (matches the data — sub-groups are quarterly only)
#   - Variant B: Stock-Watson RW common trend + sector RW trends + SV on all
#   - Gibbs sampler (~4-5 min wall time, vs Stan's hours at this N)
#
# Earlier vintages (11-group monthly mixed-freq A/B/C with Stan, and the
# 1993-onward 11-group fits) are archived in outputs/draws/legacy_11group/
# and outputs/draws/legacy_1990/.

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "dplyr", "tibble", "tidyr", "lubridate",
    "readabs", "readxl",
    "posterior", "loo",
    "ggplot2", "scales",
    "arrow", "MASS"
  ),
  format = "rds"
)

tar_source("R")
tar_source("R/gibbs")

list(
  # -------------------- Data ingest --------------------
  # SA sub-group quarterly inflation, aggregated from Appendix 1a's
  # 90 expenditure-class SA indices using the 2025 Weighting Pattern
  # EC weights. SA removes the regular annual seasonality (most
  # notably Education's Q1 fee re-indexing) that the model would
  # otherwise misattribute to outliers/sector trend.
  tar_target(cpi_q_raw,        fetch_cpi_quarterly_subgroup_sa()),
  tar_target(weights,          fetch_weights_subgroup()),
  tar_target(headline_yoy,     fetch_cpi_headline_yoy()),
  tar_target(trimmed_mean_yoy, fetch_cpi_trimmed_mean_yoy()),

  # -------------------- Inflation series & cleaning --------------------
  tar_target(cpi_q_clean, build_inflation_series(cpi_q_raw)),
  tar_target(infl_demeaned, demean_series(cpi_q_clean)),
  tar_target(stan_data,
             build_stan_data_quarterly(infl_demeaned, weights = weights)),

  # Persist stan_data to outputs/draws/ as a side-effect so the cron's
  # "rotate prior vintage" step can pick it up alongside fit_B.rds.
  tar_target(
    stan_data_file, {
      path <- "outputs/draws/stan_data.rds"
      saveRDS(stan_data, path)
      normalizePath(path)
    },
    format = "file"
  ),

  # -------------------- Variant B fit (Gibbs, quarterly clock) --------
  # File-target — value is the absolute path string, not the loaded object.
  # On cache miss the command fits-and-saves; on cache hit just returns
  # the path. ~4-5 min wall time at N=33, T=133.
  tar_target(
    fit_B, {
      path <- "outputs/draws/fit_B.rds"
      if (!file.exists(path)) {
        T_ <- stan_data$T;  N <- stan_data$N
        y  <- stan_data$y
        set.seed(20260528)
        fit <- fit_mct_gibbs(
          y, q_MA = 1L, ref = stan_data$ref,
          n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
          config = list(variant          = "B",
                        s2_gam_prior     = 1 / (60 * 4),  # quarterly cadence
                        theta_prec_prior = 10,            # tight MA prior
                        use_outliers     = TRUE)
        )
        saveRDS(fit, path)
      }
      normalizePath(path)
    },
    format = "file",
    cue = tar_cue(command = FALSE)
  ),

  # -------------------- Post-processing --------------------
  # fit_B target value is the path string; load + wrap inside each target.
  tar_target(fit_B_wrap, {
    raw <- readRDS(fit_B)
    mct_gibbs_fit(raw, stan_data)
  }),
  tar_target(trend_B, trend_path(fit_B_wrap, stan_data$dates)),
  tar_target(common_share_B,
             common_share_path(fit_B_wrap, stan_data$dates)),
  tar_target(sector_contribs_B,
             sector_contributions(fit_B_wrap, weights$weight,
                                  stan_data$dates, stan_data$groups)),

  # -------------------- LOO + fit diagnostics --------------------
  # Single-variant model: no inter-variant LOO comparison; report LOO
  # absolute fit + Pareto-k diagnostics for the production fit.
  tar_target(loo_B, {
    ll <- as.matrix(fit_B_wrap$draws(variables = "log_lik"))
    keep <- colSums(abs(ll)) > 0
    ll <- ll[, keep, drop = FALSE]
    r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1L, nrow(ll)))
    loo::loo(ll, r_eff = r_eff, cores = 1L)
  }),

  # -------------------- Dashboard --------------------
  # Manual target (instead of tar_quarto) so we can explicitly declare
  # the upstream dependencies. tar_quarto()'s static dependency scan
  # missed our multi-line `tar_load(c(...))` block on a clean runner,
  # dispatching the render before fit_B/trend_B/etc. existed.
  tar_target(
    dashboard, {
      # Touch every upstream target so it becomes a dependency.
      .deps <- list(
        stan_data, stan_data_file, weights, infl_demeaned,
        headline_yoy, trimmed_mean_yoy,
        fit_B, fit_B_wrap,
        trend_B, common_share_B, sector_contribs_B, loo_B
      )
      quarto::quarto_render("dashboard", as_job = FALSE, quiet = FALSE)
      normalizePath("dashboard/_site/index.html")
    },
    format = "file"
  )
)
