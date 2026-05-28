# Step 19 — Quick verification of sub-group data fetchers + coverage.

for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)
suppressMessages(library(dplyr))

w_sub <- fetch_weights_subgroup()
cat("Sub-group weights: n=", nrow(w_sub),
    "  sum=", sum(w_sub$weight), "\n")

cpi_sub <- fetch_cpi_quarterly_subgroup()
cat("Sub-group CPI obs: nrow=", nrow(cpi_sub), "\n")
cat("Distinct sub-groups present:", length(unique(cpi_sub$group)), "\n")

cov <- cpi_sub |>
  dplyr::group_by(.data$group) |>
  dplyr::summarise(
    start = min(.data$date), end = max(.data$date),
    n     = dplyr::n(),
    n_NA  = sum(is.na(.data$index)),
    .groups = "drop"
  )
cat("\nPer-subgroup coverage (start / end / n_obs / n_NA):\n")
for (i in seq_len(nrow(cov))) {
  cat(sprintf("  %-50s  %s -> %s  n=%3d  NA=%d\n",
              cov$group[i], cov$start[i], cov$end[i], cov$n[i], cov$n_NA[i]))
}
cat("\nLatest sub-group start date:", as.character(max(cov$start)),
    " (this is what limits the sample)\n")
cat("Sub-groups that start AFTER 1993-01:\n")
late <- cov[cov$start > as.Date("1993-01-01"), ]
for (i in seq_len(nrow(late))) {
  cat(sprintf("  %-50s  starts %s\n", late$group[i], late$start[i]))
}
