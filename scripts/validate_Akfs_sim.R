# Phase 3a optimization validation: Ac vs Akf vs Akfs on the same sim
# data. Akfs is the SPARSE implementation of the Akf marginalization —
# same math, different inner loops. Target: 5-20x speedup over Akf,
# matching posteriors.
#
# Run with: Rscript scripts/validate_Akfs_sim.R

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
cat(sprintf("[Akfs-val] simulating T=%d N=%d seed=%d\n", T_, N, seed))
sim <- simulate_mct_A(T_ = T_, N = N, seed = seed)
sd_full <- sim$stan_data
sd_full$dates  <- seq(as.Date("1990-01-01"), by = "month", length.out = T_)
sd_full$groups <- paste0("g", seq_len(N))

stopifnot(!is.unsorted(sd_full$t_m), !is.unsorted(sd_full$t_q))

run <- function(variant_) {
  cat(sprintf("\n[Akfs-val] fitting %s ...\n", variant_))
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

# Only run Akfs — Ac and Akf timings are cached from the previous
# validate_Akf_sim.R run (Ac=5.28 min, Akf=174.83 min). Same data, same
# seed, same config — so the comparison is still apples-to-apples.
# Hyperparameter posteriors for Ac/Akf from that log:
#   rho        Ac: 0.704 [0.480, 0.851]   Akf: 0.702 [0.462, 0.850]
#   mu_hc      Ac:-0.981 [-1.550,-0.419]  Akf:-0.984 [-1.532,-0.409]
#   sigma_hc   Ac: 0.054 [0.017, 0.115]   Akf: 0.055 [0.018, 0.117]
#   lambda[1]  Ac: 1.033 [0.769, 1.320]   Akf: 1.034 [0.763, 1.323]
#   lambda[2]  Ac: 0.909 [0.643, 1.197]   Akf: 0.904 [0.642, 1.180]
res_Akfs <- run("Akfs")

Ac_min  <- 5.28
Akf_min <- 174.83

cat("\n[Akfs-val] ========================================================\n")
cat(sprintf("                       %10s  %10s  %10s\n",
            "Ac (log)", "Akf (log)", "Akfs"))
cat(sprintf("  wall time (min):     %10.2f  %10.2f  %10.2f\n",
            Ac_min, Akf_min, res_Akfs$elapsed_min))
cat(sprintf("  speedup vs Akf:      %10.2fx %10s  %10.2fx\n",
            Akf_min / Ac_min, "1.00x", Akf_min / res_Akfs$elapsed_min))
cat(sprintf("  speedup vs Ac:       %10s  %10.2fx  %10.2fx\n",
            "1.00x", Ac_min / Akf_min, Ac_min / res_Akfs$elapsed_min))
cat(sprintf("  divergent transitions: %23d  %10d\n",
            0, res_Akfs$n_divergent))    # Akf had 0
cat(sprintf("  max-treedepth hits:    %23d  %10d\n",
            0, res_Akfs$n_max_td))
cat(sprintf("  min ebfmi:             %23.3f  %10.3f\n",
            0.920, res_Akfs$ebfmi_min))   # Akf was 0.920

cat("\n[Akfs-val] Akfs hyperparameter medians [16%, 84%]:\n")
cat("  (cached Akf values in header above for comparison)\n")
params_to_compare <- c("rho", "mu_hc", "sigma_hc",
                       "lambda_const[1]", "lambda_const[2]")
for (p in params_to_compare) {
  s_Akfs <- posterior::summarise_draws(res_Akfs$fit$draws(variables = p),
                                       median, ~quantile(.x, c(0.16, 0.84)))
  cat(sprintf("  %-18s  Akfs: %6.3f [%6.3f, %6.3f]\n",
              p, s_Akfs$median, s_Akfs$`16%`, s_Akfs$`84%`))
}

cat("\n[Akfs-val] done.\n")
