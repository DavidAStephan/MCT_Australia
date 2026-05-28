# Step 14b — LOO comparison: baseline Gibbs vs +outliers.
#
# Wraps both fits with mct_gibbs_fit (which exposes $draws(log_lik) with
# both monthly + quarterly entries) and uses loo::loo + loo_compare.
#
# Single-chain Gibbs, so we compute r_eff with chain_id = 1.
#
# Usage: Rscript scripts/step14b_loo_baseline_vs_outliers.R

suppressMessages({
  library(MASS)
  library(loo)
  library(posterior)
})

for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("R",       "\\.R$", full.names = TRUE)) source(f)

sd_real <- readRDS("outputs/draws/stan_data.rds")
fit_base_raw <- readRDS("outputs/draws/fit_gibbs_baseline.rds")
fit_out_raw  <- readRDS("outputs/draws/fit_gibbs_outliers.rds")

cat("=== Step 14b: LOO baseline vs +outliers ===\n")
cat("n_obs_m =", sd_real$n_obs_m, "  n_obs_q =", sd_real$n_obs_q,
    "  total =", sd_real$n_obs_m + sd_real$n_obs_q, "\n")

wrap_base <- mct_gibbs_fit(fit_base_raw, sd_real)
wrap_out  <- mct_gibbs_fit(fit_out_raw,  sd_real)

ll_base <- as.matrix(wrap_base$draws(variables = "log_lik"))
ll_out  <- as.matrix(wrap_out$draws(variables  = "log_lik"))

cat("log_lik dims: baseline ", paste(dim(ll_base), collapse = " x "),
    " | outliers ", paste(dim(ll_out), collapse = " x "), "\n", sep = "")

# Drop padded zero-only columns (defensive — shouldn't happen now that
# quarterly log_lik is computed).
keep <- (colSums(abs(ll_base)) > 0) & (colSums(abs(ll_out)) > 0)
cat("Observations with non-zero log_lik in both fits: ",
    sum(keep), " / ", length(keep), "\n", sep = "")
ll_base <- ll_base[, keep, drop = FALSE]
ll_out  <- ll_out [, keep, drop = FALSE]

# Single-chain Gibbs: chain_id = 1 throughout
chain_id_base <- rep(1L, nrow(ll_base))
chain_id_out  <- rep(1L, nrow(ll_out))

cat("\nComputing r_eff + LOO for baseline...\n")
r_eff_base <- loo::relative_eff(exp(ll_base), chain_id = chain_id_base,
                                cores = 1L)
loo_base <- loo::loo(ll_base, r_eff = r_eff_base, cores = 1L)
print(loo_base)

cat("\nComputing r_eff + LOO for outliers...\n")
r_eff_out <- loo::relative_eff(exp(ll_out), chain_id = chain_id_out,
                               cores = 1L)
loo_out <- loo::loo(ll_out, r_eff = r_eff_out, cores = 1L)
print(loo_out)

cat("\n=== loo_compare(baseline, outliers) ===\n")
cmp <- loo::loo_compare(list(baseline = loo_base, outliers = loo_out))
print(cmp)

# Headline numbers
preferred <- rownames(cmp)[1L]
elpd_diff <- cmp[2L, "elpd_diff"]
se_diff   <- cmp[2L, "se_diff"]
cat(sprintf(
  "\nHeadline: %s preferred. elpd_diff = %.1f (SE %.1f, |Δ|/SE = %.2f)\n",
  preferred, elpd_diff, se_diff,
  if (se_diff > 0) abs(elpd_diff) / se_diff else NA_real_
))

# Pareto-k summary
cat("\nPareto-k diagnostics:\n")
cat(sprintf("  baseline: k > 0.7 in %d / %d obs\n",
            sum(loo_base$diagnostics$pareto_k > 0.7),
            length(loo_base$diagnostics$pareto_k)))
cat(sprintf("  outliers: k > 0.7 in %d / %d obs\n",
            sum(loo_out$diagnostics$pareto_k > 0.7),
            length(loo_out$diagnostics$pareto_k)))

# Cache the loo objects for re-inspection
saveRDS(list(baseline = loo_base, outliers = loo_out, compare = cmp),
        "outputs/draws/loo_baseline_vs_outliers.rds")
cat("\nSaved -> outputs/draws/loo_baseline_vs_outliers.rds\n")
