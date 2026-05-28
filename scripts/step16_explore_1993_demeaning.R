# Step 16 — Explore the 1993-onward sample and alternative demeaning
# windows.
#
# Goal: see how the per-sector demeaned series change under different
# windows when the latent monthly clock starts at 1993-01 instead of
# 1990-01. Cheap (no Stan fitting) — informs the choice of demeaning
# window before we re-fit Variant A on the new sample.

suppressMessages({
  library(targets)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)

# Pull cached inflation series from the targets store.
cpi_q_clean <- tar_read(cpi_q_clean)
cpi_m_clean <- tar_read(cpi_m_clean)
weights     <- tar_read(weights)

combined <- dplyr::bind_rows(cpi_q_clean, cpi_m_clean)

windows <- list(
  "2000-2019 (current default)" = c(2000, 2019),
  "1993-2019 (full pre-pandemic)" = c(1993, 2019),
  "1995-2019" = c(1995, 2019),
  "2000-2024 (excludes only 2020-21 pandemic shock)" = c(2000, 2024)
)

cat("=== Step 16: per-sector means under candidate demeaning windows ===\n")
cat("(Means computed on quarterly series — these are the values that\n",
    " get subtracted from every obs of the corresponding group)\n\n",
    sep = "")

for (wname in names(windows)) {
  w <- windows[[wname]]
  m <- combined |>
    dplyr::filter(.data$source == "quarterly",
                  lubridate::year(.data$date) >= w[1],
                  lubridate::year(.data$date) <= w[2]) |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(mean_infl = mean(.data$inflation, na.rm = TRUE),
                     n = dplyr::n(),
                     .groups = "drop")
  cat(sprintf("--- Window: %s (n_obs per series = %d) ---\n", wname, m$n[1L]))
  for (i in seq_len(nrow(m))) {
    cat(sprintf("  %-50s %6.3f\n",
                as.character(m$group[i]), m$mean_infl[i]))
  }
  cat(sprintf("  %-50s %6.3f (weighted)\n", "[weighted avg]",
              sum(m$mean_infl * (weights$weight / 100))))
  cat("\n")
}

cat("\n=== Build a 1993-onward stan_data with the current demeaning ===\n")
infl_demeaned_current <- demean_series(combined, window = c(2000, 2019))
sd_1993 <- build_stan_data(
  infl_demeaned_current, weights = weights,
  start = as.Date("1993-01-01")
)
cat(sprintf("1993-onward: T=%d  N=%d  n_obs_m=%d  n_obs_q=%d  dates: %s -> %s\n",
            sd_1993$T, sd_1993$N, sd_1993$n_obs_m, sd_1993$n_obs_q,
            min(sd_1993$dates), max(sd_1993$dates)))

# Also build the original 1990-onward for comparison
sd_1990 <- build_stan_data(
  infl_demeaned_current, weights = weights,
  start = as.Date("1990-01-01")
)
cat(sprintf("1990-onward: T=%d  N=%d  n_obs_m=%d  n_obs_q=%d  dates: %s -> %s\n",
            sd_1990$T, sd_1990$N, sd_1990$n_obs_m, sd_1990$n_obs_q,
            min(sd_1990$dates), max(sd_1990$dates)))

cat(sprintf("\nObservations dropped by 1990->1993 cut: %d quarterly + %d monthly = %d total\n",
            sd_1990$n_obs_q - sd_1993$n_obs_q,
            sd_1990$n_obs_m - sd_1993$n_obs_m,
            (sd_1990$n_obs_q + sd_1990$n_obs_m) -
              (sd_1993$n_obs_q + sd_1993$n_obs_m)))

# Sanity check: distribution of inflation values in the dropped window
dropped <- combined |>
  dplyr::filter(.data$date >= as.Date("1990-01-01"),
                .data$date <  as.Date("1993-01-01"),
                .data$source == "quarterly") |>
  dplyr::summarise(
    n_obs   = dplyr::n(),
    mean    = mean(.data$inflation, na.rm = TRUE),
    median  = median(.data$inflation, na.rm = TRUE),
    sd      = sd(.data$inflation, na.rm = TRUE),
    p05     = quantile(.data$inflation, 0.05, na.rm = TRUE),
    p95     = quantile(.data$inflation, 0.95, na.rm = TRUE),
    max_abs = max(abs(.data$inflation), na.rm = TRUE)
  )
cat("\n1990-1992 quarterly inflation stats (the dropped period):\n")
print(dropped)

# Same for 1993-1996 (the "buffer" period after the cut)
buffer <- combined |>
  dplyr::filter(.data$date >= as.Date("1993-01-01"),
                .data$date <  as.Date("1996-01-01"),
                .data$source == "quarterly") |>
  dplyr::summarise(
    n_obs   = dplyr::n(),
    mean    = mean(.data$inflation, na.rm = TRUE),
    median  = median(.data$inflation, na.rm = TRUE),
    sd      = sd(.data$inflation, na.rm = TRUE),
    p05     = quantile(.data$inflation, 0.05, na.rm = TRUE),
    p95     = quantile(.data$inflation, 0.95, na.rm = TRUE)
  )
cat("\n1993-1995 quarterly inflation stats (kept, the immediate post-targeting):\n")
print(buffer)

# Cache the candidate stan_data for the next step
dir.create("outputs/draws", showWarnings = FALSE, recursive = TRUE)
saveRDS(sd_1993, "outputs/draws/stan_data_1993.rds")
cat("\nSaved: outputs/draws/stan_data_1993.rds (start=1993-01, demean=2000-2019)\n")
