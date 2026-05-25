# Phase 3 validation: fit Akf (KF-marginalized) on the same sim data we
# used for A/Ac/Act. Goal: confirm Akf math is correct AND measure
# real wall-time gain vs Ac.
#
# Run with: Rscript scripts/validate_Akf_sim.R
#
# Two correctness signals:
#   1. Akf must converge (no divergences, R-hat near 1, ebfmi > 0.3).
#   2. Akf's posterior medians + intervals for the HYPERPARAMETERS
#      (rho, sigma_hc, mu_hc, lambda_const, etc.) should match Ac's
#      within Monte Carlo error. Both target the same marginal density.

suppressMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

{
  .lines <- readLines("tests/testthat/test-stan-simulated.R")
  .start <- min(grep("^simulate_mct_A <-", .lines))
  .end <- min(grep("^test_that", .lines)) - 1
  eval(parse(text = paste(.lines[.start:.end], collapse = "\n")))
}

T_ <- 120; N <- 5; seed <- 42
cat(sprintf("[Akf-val] simulating T=%d N=%d seed=%d\n", T_, N, seed))
sim <- simulate_mct_A(T_ = T_, N = N, seed = seed)
sd_full <- sim$stan_data
sd_full$dates  <- seq(as.Date("1990-01-01"), by = "month", length.out = T_)
sd_full$groups <- paste0("g", seq_len(N))

# Verify the Akf pre-condition: obs arrays sorted ascending by t.
stopifnot("monthly obs not sorted by t" =
            !is.unsorted(sd_full$t_m, strictly = FALSE))
stopifnot("quarterly obs not sorted by t" =
            !is.unsorted(sd_full$t_q, strictly = FALSE))

run <- function(variant_) {
  cat(sprintf("\n[Akf-val] fitting %s ...\n", variant_))
  t0 <- Sys.time()
  fit <- fit_mct(
    sd_full, variant = variant_,
    chains = 4, parallel_chains = 4,
    iter_warmup = 750, iter_sampling = 750,
    adapt_delta = 0.97, max_treedepth = 12,
    seed = 20260525, save_dir = NULL
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  diag <- fit$diagnostic_summary()
  list(
    variant = variant_, elapsed_min = elapsed,
    n_divergent = sum(diag$num_divergent),
    n_max_td    = sum(diag$num_max_treedepth),
    ebfmi_min   = min(diag$ebfmi),
    fit = fit
  )
}

res_Ac  <- run("Ac")
res_Akf <- run("Akf")

cat("\n[Akf-val] =========================================================\n")
cat(sprintf("                            %10s  %10s  %10s\n",
            "Ac", "Akf", "speedup"))
cat(sprintf("  wall time (min):          %10.2f  %10.2f  %9.2fx\n",
            res_Ac$elapsed_min, res_Akf$elapsed_min,
            res_Ac$elapsed_min / res_Akf$elapsed_min))
cat(sprintf("  divergent transitions:    %10d  %10d\n",
            res_Ac$n_divergent, res_Akf$n_divergent))
cat(sprintf("  max-treedepth hits:       %10d  %10d\n",
            res_Ac$n_max_td, res_Akf$n_max_td))
cat(sprintf("  min ebfmi across chains:  %10.3f  %10.3f\n",
            res_Ac$ebfmi_min, res_Akf$ebfmi_min))

# Hyperparameter posterior comparison: both Akf and Ac target the same
# marginal density over (hyperparams, h, lambda). If math is right,
# medians + intervals should be within MC noise of each other.
cat("\n[Akf-val] Hyperparameter posterior medians [16%, 84%] — should match:\n")
params_to_compare <- c("rho", "mu_hc", "sigma_hc",
                       "lambda_const[1]", "lambda_const[2]")
for (p in params_to_compare) {
  a  <- posterior::summarise_draws(res_Ac$fit$draws(variables = p),
                                    median, ~quantile(.x, c(0.16, 0.84)))
  k  <- posterior::summarise_draws(res_Akf$fit$draws(variables = p),
                                    median, ~quantile(.x, c(0.16, 0.84)))
  cat(sprintf("  %-22s  Ac:%6.3f [%6.3f, %6.3f]   Akf:%6.3f [%6.3f, %6.3f]\n",
              p,
              a$median, a$`16%`, a$`84%`,
              k$median, k$`16%`, k$`84%`))
}

# R-hat check on the same params (highest indicator of mixing problems)
cat("\n[Akf-val] R-hat (should be < 1.01):\n")
for (p in params_to_compare) {
  rh_a <- posterior::rhat(res_Ac$fit$draws(variables = p))
  rh_k <- posterior::rhat(res_Akf$fit$draws(variables = p))
  cat(sprintf("  %-22s  Ac: %5.3f   Akf: %5.3f\n", p, rh_a, rh_k))
}

cat("\n[Akf-val] done.\n")
