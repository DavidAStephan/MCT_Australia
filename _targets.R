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

  # -------------------- Stan fits (slow; cached by hash of stan_data) ------
  # `save_dir = NULL` because targets handles caching itself.
  tar_target(fit_A,
             fit_mct(stan_data, variant = "A", save_dir = NULL),
             format = "rds"),
  tar_target(fit_B,
             fit_mct(stan_data, variant = "B", save_dir = NULL),
             format = "rds"),

  # -------------------- Post-processing --------------------
  tar_target(trend_A, trend_path(fit_A, stan_data$dates)),
  tar_target(trend_B, trend_path(fit_B, stan_data$dates)),
  tar_target(common_share_A, common_share_path(fit_A, stan_data$dates)),
  tar_target(common_share_B, common_share_path(fit_B, stan_data$dates)),
  tar_target(sector_contribs_A,
             sector_contributions(fit_A, weights$weight,
                                  stan_data$dates, stan_data$groups)),
  tar_target(sector_contribs_B,
             sector_contributions(fit_B, weights$weight,
                                  stan_data$dates, stan_data$groups)),

  # -------------------- LOO variant comparison --------------------
  tar_target(loo_A, loo_for_fit(fit_A)),
  tar_target(loo_B, loo_for_fit(fit_B)),
  tar_target(variant_comparison, list(
    loo_compare = compare_fits_loo(loo_A, loo_B),
    text = variant_comparison_text(loo_A, loo_B)
  )),

  # -------------------- Dashboard --------------------
  # Quarto reads the targets it needs via `tar_read()` inside the .qmd file.
  tar_quarto(dashboard, path = "dashboard/index.qmd")
)
