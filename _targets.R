# targets pipeline for MCT_Australia. Run with:
#   targets::tar_make()
# Re-fits are content-addressed: the Stan fits (~30-90 min each) only re-run
# when stan_data hashes differ.

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "dplyr", "tibble", "tidyr", "lubridate",
    "readabs", "readxl",
    "cmdstanr", "posterior", "loo",
    "ggplot2", "scales",
    "arrow"
  ),
  format = "rds"  # default; specific targets override (see fit_*)
)

# Source every helper in R/ into the targets workers.
tar_source("R")

list(
  # -------------------- Data ingest --------------------
  tar_target(cpi_q_raw, fetch_cpi_quarterly()),
  tar_target(cpi_m_raw, fetch_cpi_monthly()),
  tar_target(weights,   fetch_weights()),

  # -------------------- Inflation series & cleaning --------------------
  tar_target(cpi_q_clean, build_inflation_series(cpi_q_raw)),
  tar_target(cpi_m_clean, build_inflation_series(cpi_m_raw)),
  tar_target(infl_demeaned, {
    combined <- dplyr::bind_rows(cpi_q_clean, cpi_m_clean)
    demean_series(combined)
  }),
  tar_target(stan_data, build_stan_data(infl_demeaned, weights = weights)),

  # -------------------- Stan fits --------------------
  # fit_X targets are `format = "file"` — their *value* is the absolute path
  # string, not the loaded CmdStanMCMC object. This avoids duplicating
  # ~7 GB of fits into _targets/objects/. Downstream targets (and the
  # dashboard) call readRDS(fit_X) when they need the object; it's held in
  # worker memory only for the duration of that target.
  #
  # The command fits-and-saves on cache miss. On the local workstation
  # (where the rds files already exist from standalone runs) targets just
  # validates the file hash and skips. On a clean runner the command runs
  # fit_mct() and saves the rds. When stan_data changes (new ABS vintage)
  # targets invalidates and refits.
  #
  # `cue = tar_cue(command = FALSE)` prevents command-text changes from
  # triggering refits — useful during pipeline refactors so cached fits
  # aren't invalidated by cosmetic edits to this block.
  tar_target(
    fit_A, {
      path <- "outputs/draws/fit_A.rds"
      if (!file.exists(path)) fit_mct(stan_data, variant = "A",
                                      save_dir = "outputs/draws")
      normalizePath(path)
    },
    format = "file",
    cue = tar_cue(command = FALSE)
  ),
  tar_target(
    fit_B, {
      path <- "outputs/draws/fit_B.rds"
      if (!file.exists(path)) fit_mct(stan_data, variant = "B",
                                      save_dir = "outputs/draws")
      normalizePath(path)
    },
    format = "file",
    cue = tar_cue(command = FALSE)
  ),
  tar_target(
    fit_C, {
      path <- "outputs/draws/fit_C.rds"
      if (!file.exists(path)) fit_mct(stan_data, variant = "C",
                                      save_dir = "outputs/draws")
      normalizePath(path)
    },
    format = "file",
    cue = tar_cue(command = FALSE)
  ),

  # -------------------- Post-processing --------------------
  tar_target(trend_A, trend_path(readRDS(fit_A), stan_data$dates)),
  tar_target(trend_B, trend_path(readRDS(fit_B), stan_data$dates)),
  tar_target(trend_C, trend_path(readRDS(fit_C), stan_data$dates)),
  tar_target(common_share_A,
             common_share_path(readRDS(fit_A), stan_data$dates)),
  tar_target(common_share_B,
             common_share_path(readRDS(fit_B), stan_data$dates)),
  tar_target(common_share_C,
             common_share_path(readRDS(fit_C), stan_data$dates)),
  tar_target(sector_contribs_A,
             sector_contributions(readRDS(fit_A), weights$weight,
                                  stan_data$dates, stan_data$groups)),
  tar_target(sector_contribs_B,
             sector_contributions(readRDS(fit_B), weights$weight,
                                  stan_data$dates, stan_data$groups)),
  tar_target(sector_contribs_C,
             sector_contributions(readRDS(fit_C), weights$weight,
                                  stan_data$dates, stan_data$groups)),

  # -------------------- LOO variant comparison --------------------
  tar_target(loo_A, loo_for_fit(readRDS(fit_A))),
  tar_target(loo_B, loo_for_fit(readRDS(fit_B))),
  tar_target(loo_C, loo_for_fit(readRDS(fit_C))),
  tar_target(variant_comparison, {
    loo_list <- list(A = loo_A, B = loo_B, C = loo_C)
    list(
      loo_compare = compare_fits_loo(loo_list),
      text = variant_comparison_text(loo_list)
    )
  }),

  # -------------------- Dashboard --------------------
  # Quarto reads the targets it needs via `tar_read()` inside the .qmd file.
  tar_quarto(dashboard, path = "dashboard", working_directory = "dashboard")
)
