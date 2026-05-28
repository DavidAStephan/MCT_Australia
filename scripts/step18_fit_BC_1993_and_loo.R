# Step 18 — Fit Variants B and C on the 1993-onward sample, then run
# the 3-way LOO comparison. Companion to step17 which already produced
# fit_A_1993.rds.
#
# Saves under non-default names so the production fits (fit_B.rds /
# fit_C.rds, used by the monthly cron) remain untouched.

suppressMessages({
  library(cmdstanr)
  library(posterior)
  library(loo)
})
for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

stan_data_1993 <- readRDS("outputs/draws/stan_data_1993.rds")
cat(sprintf(
  "[step18] stan_data_1993: T=%d  N=%d  n_obs_m=%d  n_obs_q=%d\n",
  stan_data_1993$T, stan_data_1993$N,
  stan_data_1993$n_obs_m, stan_data_1993$n_obs_q
))

# --------------------------------------------------------------
# Variant B (RW common trend) — usually faster than A (no rho)
# --------------------------------------------------------------
fit_B_path <- "outputs/draws/fit_B_1993.rds"
if (file.exists(fit_B_path)) {
  cat("\n[step18] B cached at ", fit_B_path, " — loading.\n", sep = "")
  fit_B_1993 <- readRDS(fit_B_path)
} else {
  cat("\n[step18] Fitting Variant B on 1993 sample...\n")
  t0 <- Sys.time()
  fit_B_1993 <- fit_mct(
    stan_data_1993,
    variant         = "B",
    chains          = 4, parallel_chains = 4,
    iter_warmup     = 1500, iter_sampling = 1500,
    adapt_delta     = 0.95, max_treedepth = 12,
    seed            = 19930102,
    save_dir        = NULL
  )
  fit_B_1993$save_object(fit_B_path)
  cat(sprintf("[step18] B done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  cat("[step18] B diagnostics:\n"); print(fit_B_1993$diagnostic_summary())
}

# --------------------------------------------------------------
# Variant C (combined RW common trend + white-noise common cycle)
# --------------------------------------------------------------
fit_C_path <- "outputs/draws/fit_C_1993.rds"
if (file.exists(fit_C_path)) {
  cat("\n[step18] C cached at ", fit_C_path, " — loading.\n", sep = "")
  fit_C_1993 <- readRDS(fit_C_path)
} else {
  cat("\n[step18] Fitting Variant C on 1993 sample...\n")
  t0 <- Sys.time()
  fit_C_1993 <- fit_mct(
    stan_data_1993,
    variant         = "C",
    chains          = 4, parallel_chains = 4,
    iter_warmup     = 1500, iter_sampling = 1500,
    adapt_delta     = 0.95, max_treedepth = 12,
    seed            = 19930103,
    save_dir        = NULL
  )
  fit_C_1993$save_object(fit_C_path)
  cat(sprintf("[step18] C done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  cat("[step18] C diagnostics:\n"); print(fit_C_1993$diagnostic_summary())
}

# --------------------------------------------------------------
# 3-way LOO on the 1993 sample
# --------------------------------------------------------------
fit_A_1993 <- readRDS("outputs/draws/fit_A_1993.rds")

cat("\n[step18] Computing LOO for A_1993, B_1993, C_1993 ...\n")
loo_A_1993 <- loo_for_fit(fit_A_1993)
loo_B_1993 <- loo_for_fit(fit_B_1993)
loo_C_1993 <- loo_for_fit(fit_C_1993)

cat("\n=== loo_compare on 1993 sample ===\n")
loo_list_1993 <- list(A = loo_A_1993, B = loo_B_1993, C = loo_C_1993)
cmp_1993 <- compare_fits_loo(loo_list_1993)
print(cmp_1993)

cat("\n", variant_comparison_text(loo_list_1993), "\n", sep = "")

# --------------------------------------------------------------
# Side-by-side: 1990 sample vs 1993 sample
# --------------------------------------------------------------
cat("\n=== 1990-sample LOO (from cached targets store) ===\n")
loo_A_1990 <- targets::tar_read(loo_A)
loo_B_1990 <- targets::tar_read(loo_B)
loo_C_1990 <- targets::tar_read(loo_C)
cmp_1990 <- compare_fits_loo(list(A = loo_A_1990, B = loo_B_1990,
                                  C = loo_C_1990))
print(cmp_1990)

# Save the 1993 LOO objects for re-inspection without recomputing
saveRDS(list(A = loo_A_1993, B = loo_B_1993, C = loo_C_1993,
             compare = cmp_1993),
        "outputs/draws/loo_1993.rds")
cat("\n[step18] Saved: outputs/draws/loo_1993.rds\n")

# Pareto-k summary per variant
cat("\n=== Pareto-k diagnostics (k > 0.7 flagged) ===\n")
for (name in c("A", "B", "C")) {
  k_1990 <- list(A = loo_A_1990, B = loo_B_1990, C = loo_C_1990)[[name]]$diagnostics$pareto_k
  k_1993 <- loo_list_1993[[name]]$diagnostics$pareto_k
  cat(sprintf("  %s:  1990 sample = %d / %d  |  1993 sample = %d / %d\n",
              name, sum(k_1990 > 0.7), length(k_1990),
              sum(k_1993 > 0.7), length(k_1993)))
}

cat("\n[step18] done.\n")
