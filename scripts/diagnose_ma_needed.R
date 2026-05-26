# Do AU CPI obs residuals exhibit MA(q) structure?
# If YES → invest in MA(3) port; the model gap is real for AU.
# If NO → MA(3) port is unnecessary; skip directly to outliers.
#
# Diagnostic: load the cached Stan Variant A posterior, compute
# residuals u_{i,t} = (y_{i,t} - lambda_{i,t} c_t - s_{i,t}) / sigma_eps_{i,t},
# and inspect their autocorrelations (Ljung-Box test + sample ACF).
#
# Strategy: compare ACFs to a 95% noise band 1.96/sqrt(T_obs).
# If multiple lags significantly outside the band per sector, MA(3)
# is justified. If most ACFs sit inside the band, skip the port.

suppressMessages({
  library(MASS); library(posterior); library(dplyr)
})
for (f in list.files("R", "\\.R$", full.names = TRUE, recursive = TRUE)) {
  source(f)
}

cat("[ma-diag] loading real stan_data + Stan fit_A ...\n")
sd_real <- readRDS("outputs/draws/stan_data.rds")
fit_A   <- readRDS("outputs/draws/fit_A.rds")
T_ <- sd_real$T; N <- sd_real$N

# Use POSTERIOR MEDIANS for the latent paths (point estimate of residuals)
lambda_med <- apply(
  fit_A$draws(variables = "lambda", format = "draws_matrix"), 2, median
)
# Reshape: cols are "lambda[1,1]", "lambda[2,1]", ..., "lambda[T,1]", "lambda[1,2]", ...
lambda_mat <- matrix(lambda_med, T_, N, byrow = FALSE)

c_med <- apply(
  fit_A$draws(variables = "c", format = "draws_matrix"), 2, median
)
s_med <- apply(
  fit_A$draws(variables = "sector_trend", format = "draws_matrix"), 2, median
)
s_mat <- matrix(s_med, T_, N, byrow = FALSE)

# h_e is the log-variance for measurement noise; sigma_eps = exp(h_e/2)
he_med <- apply(
  fit_A$draws(variables = "h_e", format = "draws_matrix"), 2, median
)
sigma_eps_mat <- matrix(exp(he_med / 2), T_, N, byrow = FALSE)

# Compute whitened residuals at MONTHLY obs times only
# (residuals at quarterly obs times conflate multi-month avg, less clean)
u_mat <- matrix(NA_real_, T_, N)
for (k in seq_along(sd_real$y_m)) {
  t <- sd_real$t_m[k]; i <- sd_real$i_m[k]
  u_mat[t, i] <- (sd_real$y_m[k] - lambda_mat[t, i] * c_med[t] - s_mat[t, i]) /
                  sigma_eps_mat[t, i]
}

cat(sprintf("[ma-diag] computed %d whitened residuals across %d sectors\n",
            sum(!is.na(u_mat)), N))

cat("\n[ma-diag] === Per-sector ACF analysis (only sectors with >=12 obs) ===\n")
cat("Lag-1 to lag-5 ACF; |ACF| > 0.20 is suggestive at typical sample sizes.\n\n")
cat(sprintf("%-22s %5s   %s\n", "Sector",
            "Nobs", "ACF[1..5]"))
for (i in seq_len(N)) {
  u_i <- u_mat[, i]
  u_i <- u_i[!is.na(u_i)]
  if (length(u_i) < 12) next
  acf_vals <- as.numeric(acf(u_i, lag.max = 5, plot = FALSE)$acf)[-1]
  cat(sprintf("%-22s %5d   %s\n",
              sd_real$groups[i], length(u_i),
              paste(sprintf("% .2f", acf_vals), collapse = " ")))
}

# Aggregate Ljung-Box across sectors
cat("\n[ma-diag] === Aggregate Ljung-Box test for autocorrelation through lag 3 ===\n")
cat("p-value < 0.05 ⇒ significant autocorrelation ⇒ MA(3) port justified\n")
for (i in seq_len(N)) {
  u_i <- u_mat[, i]
  u_i <- u_i[!is.na(u_i)]
  if (length(u_i) < 12) next
  lb <- Box.test(u_i, lag = 3, type = "Ljung-Box")
  marker <- ifelse(lb$p.value < 0.05, "*", " ")
  cat(sprintf("%s  %-22s  Q = %5.2f   p = %.3f\n",
              marker, sd_real$groups[i], lb$statistic, lb$p.value))
}
