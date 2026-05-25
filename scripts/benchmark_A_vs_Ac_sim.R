# Head-to-head: Variant A vs Variant Ac on the SAME small sim dataset.
# Fast (~5-8 min total) and apples-to-apples — same data, same seed for the
# sampler, same chains × iters. The only difference is the model structure
# (A: time-varying lambda RW; Ac: constant lambda).
#
# Run with: Rscript scripts/benchmark_A_vs_Ac_sim.R
#
# Why this beats the real-data benchmark:
#   - 5 min not 90 min, so we can iterate fast on speed hypotheses
#   - Same data → identical posterior geometry baseline → wall-time delta
#     attributable purely to model structure
#   - T=120 is the same size as the existing test-stan-simulated.R suite

suppressMessages({
  library(cmdstanr)
  library(posterior)
})
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

# Extract just the simulator function defs from the test file (skip the
# library/skip preamble and the test_that blocks that would auto-fit).
{
  .lines <- readLines("tests/testthat/test-stan-simulated.R")
  .start <- min(grep("^simulate_mct_A <-", .lines))
  .end <- min(grep("^test_that", .lines)) - 1
  eval(parse(text = paste(.lines[.start:.end], collapse = "\n")))
}

T_ <- 120; N <- 5; seed <- 42
cat(sprintf("[benchmark] simulating T=%d N=%d seed=%d (Variant A truth)\n",
            T_, N, seed))
sim <- simulate_mct_A(T_ = T_, N = N, seed = seed)
sd_full <- sim$stan_data
# build_stan_data() returns a list with extra fields (dates, groups);
# the sim returns only the Stan-required ones, so add a placeholder.
sd_full$dates  <- seq(as.Date("1990-01-01"), by = "month", length.out = T_)
sd_full$groups <- paste0("g", seq_len(N))

# Short config — enough for timing + basic diagnostics.
chains <- 4
iter_warmup <- 750
iter_sampling <- 750

run <- function(variant_, threads_per_chain = NULL) {
  label <- if (is.null(threads_per_chain)) variant_ else
    sprintf("%s (x%d threads/chain)", variant_, threads_per_chain)
  cat(sprintf("\n[benchmark] fitting variant %s ...\n", label))
  t0 <- Sys.time()
  fit <- fit_mct(
    sd_full,
    variant = variant_,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = 0.97,
    max_treedepth = 12,
    seed = 20260525,
    threads_per_chain = threads_per_chain,
    save_dir = NULL
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  diag <- fit$diagnostic_summary()
  list(
    variant = label,
    elapsed_min = elapsed,
    n_divergent = sum(diag$num_divergent),
    n_max_td    = sum(diag$num_max_treedepth),
    ebfmi_min   = min(diag$ebfmi),
    fit = fit
  )
}

res_A   <- run("A")
res_Ac  <- run("Ac")
res_Act <- run("Act", threads_per_chain = 2)

cat("\n[benchmark] =========================================================\n")
cat(sprintf("                            %10s  %10s  %22s\n",
            "A", "Ac", "Act (2 threads)"))
cat(sprintf("  wall time (min):          %10.2f  %10.2f  %22.2f\n",
            res_A$elapsed_min, res_Ac$elapsed_min, res_Act$elapsed_min))
cat(sprintf("  speedup vs A:             %10s  %10.2fx  %22.2fx\n",
            "1.00x",
            res_A$elapsed_min / res_Ac$elapsed_min,
            res_A$elapsed_min / res_Act$elapsed_min))
cat(sprintf("  divergent transitions:    %10d  %10d  %22d\n",
            res_A$n_divergent, res_Ac$n_divergent, res_Act$n_divergent))
cat(sprintf("  max-treedepth hits:       %10d  %10d  %22d\n",
            res_A$n_max_td, res_Ac$n_max_td, res_Act$n_max_td))
cat(sprintf("  min ebfmi across chains:  %10.3f  %10.3f  %22.3f\n",
            res_A$ebfmi_min, res_Ac$ebfmi_min, res_Act$ebfmi_min))

cat("\n[benchmark] Recovery check (latest-t trend, truth = ",
    sprintf("%.3f", tail(sim$truth$trend, 1)), "):\n", sep = "")
for (r in list(res_A, res_Ac, res_Act)) {
  tr <- trend_path(r$fit, sd_full$dates)
  cat(sprintf("  %s: median=%6.3f  68%% band=[%6.3f, %6.3f]\n",
              r$variant, tail(tr$median, 1),
              tail(tr$lower, 1), tail(tr$upper, 1)))
}

cat("\n[benchmark] done.\n")
