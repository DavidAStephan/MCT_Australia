# Fit both MCT variants on the full real ABS panel. Designed to run as a
# background batch:
#
#   nohup Rscript scripts/run_real_fits.R > outputs/fit_run.log 2>&1 &
#
# Writes:
#   outputs/draws/fit_A.rds              (full CmdStanMCMC)
#   outputs/draws/fit_A_slim.parquet     (key generated quantities only)
#   outputs/draws/fit_B.rds
#   outputs/draws/fit_B_slim.parquet
#   outputs/draws/stan_data.rds          (the assembled input, for later
#                                         post-processing without re-fetching)

suppressMessages({
  library(dplyr); library(tibble); library(lubridate)
  library(readabs); library(readxl)
  library(cmdstanr); library(posterior); library(arrow)
})

source("R/data_fetch.R")
source("R/data_prep.R")
source("R/fit_stan.R")

dir.create("outputs/draws", showWarnings = FALSE, recursive = TRUE)

message("=== Building stan_data from ABS ===  ", Sys.time())
stan_data <- build_real_stan_data()
saveRDS(stan_data, "outputs/draws/stan_data.rds")
message(sprintf(
  "  T=%d  N=%d  n_obs_m=%d  n_obs_q=%d  date range %s to %s",
  stan_data$T, stan_data$N, stan_data$n_obs_m, stan_data$n_obs_q,
  as.character(min(stan_data$dates)), as.character(max(stan_data$dates))
))

run_one <- function(variant) {
  message("\n=== Fitting Variant ", variant, " ===  ", Sys.time())
  fit <- tryCatch(
    fit_mct(stan_data, variant = variant),
    error = function(e) {
      message("[", variant, "] FIT FAILED: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(fit)) {
    tryCatch(
      export_fit_parquet(fit, variant = variant),
      error = function(e) message("[", variant, "] parquet export failed: ",
                                  conditionMessage(e))
    )
    # Sampling diagnostics summary
    diag <- fit$diagnostic_summary()
    n_div <- sum(diag$num_divergent)
    n_tree <- sum(diag$num_max_treedepth)
    summary_rhat <- fit$summary(c("trend", "common_share"),
                                rhat = posterior::rhat)
    max_rhat <- max(summary_rhat$rhat, na.rm = TRUE)
    message(sprintf(
      "[%s] divergent=%d  max_treedepth_hits=%d  max_rhat(trend,common_share)=%.4f",
      variant, n_div, n_tree, max_rhat
    ))
  }
  fit
}

fit_A <- run_one("A")
fit_B <- run_one("B")

message("\n=== Done at ", Sys.time(), " ===")
