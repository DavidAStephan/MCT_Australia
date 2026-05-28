# Step 20 — Build sub-group stan_data (quarterly-only, 1993-onward).
# Verify obs_type distribution and sample shape before kicking off any fits.

for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)
suppressMessages({library(dplyr); library(lubridate)})

cpi_sub  <- fetch_cpi_quarterly_subgroup()
w_sub    <- fetch_weights_subgroup()
infl_sub <- build_inflation_series(cpi_sub) |>
  demean_series(window = c(1993, 2019))

stan_data_sg <- build_stan_data(
  infl_sub,
  weights = w_sub,
  ref = "Rents",
  start = as.Date("1993-01-01")
)
saveRDS(stan_data_sg, "outputs/draws/stan_data_subgroup.rds")

cat(sprintf(
  "[step20] Sub-group stan_data: T=%d  N=%d  n_obs_m=%d  n_obs_q=%d\n",
  stan_data_sg$T, stan_data_sg$N,
  stan_data_sg$n_obs_m, stan_data_sg$n_obs_q
))
cat("Date range:", as.character(min(stan_data_sg$dates)),
    "->", as.character(max(stan_data_sg$dates)), "\n")
cat("Reference sub-group (ref index ", stan_data_sg$ref, "): ",
    stan_data_sg$groups[stan_data_sg$ref], "\n", sep = "")

# obs counts per sub-group
n_obs_per <- table(stan_data_sg$i_q)
cat("\nQuarterly obs per sub-group (any >0 missing data?):\n")
for (i in seq_len(stan_data_sg$N)) {
  cnt <- if (as.character(i) %in% names(n_obs_per)) n_obs_per[as.character(i)] else 0L
  cat(sprintf("  %-50s  n_q_obs=%3d\n",
              stan_data_sg$groups[i], cnt))
}
cat("\nTotal n_obs_q =", stan_data_sg$n_obs_q,
    " (cf. max possible at T=", stan_data_sg$T,
    ", N=", stan_data_sg$N, ": about ",
    stan_data_sg$N * floor(stan_data_sg$T / 3), ")\n", sep = "")
