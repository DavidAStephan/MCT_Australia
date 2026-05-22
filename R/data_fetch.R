# Data fetchers for the MCT_Australia pipeline.
#
# All Australian CPI data now lives in ABS catalogue 6401.0:
#   - Table 2  : monthly group index numbers, weighted average of 8 capital
#                cities (the new complete Monthly CPI, back-series to Apr 2024)
#   - Table 18 : quarterly group, sub-group and expenditure class index
#                numbers, weighted average of 8 capital cities (history back
#                to 1980s; what we use for the pre-Apr-2024 sample)
#
# 2025 expenditure weights are pulled separately from the ABS "Weighting
# Pattern" workbook (already cached in data/raw/cpi_weights_2025.xlsx; see
# fetch_weights()).
#
# We rely on {readabs} for the time-series downloads — it handles the ABS
# spreadsheet plumbing and caches under Sys.getenv("R_READABS_PATH").

#' The 11 standard ABS CPI groups, ordered with Housing first so it can be
#' used as the identification reference sector ($\lambda_{Housing,t} = 1$).
abs_cpi_groups <- function() {
  c(
    "Housing",
    "Food and non-alcoholic beverages",
    "Alcohol and tobacco",
    "Clothing and footwear",
    "Furnishings, household equipment and services",
    "Health",
    "Transport",
    "Communication",
    "Recreation and culture",
    "Education",
    "Insurance and financial services"
  )
}

#' Pull quarterly group-level CPI index numbers (weighted average of 8 capital
#' cities) from ABS 6401.0 Table 18.
#'
#' @return A tibble with columns `date` (quarter-end date, class Date), `group`
#'   (one of `abs_cpi_groups()`), `index` (numeric), `source` ("quarterly").
#'   One row per (group, quarter).
fetch_cpi_quarterly <- function(cache_dir = "data/raw") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(R_READABS_PATH = normalizePath(cache_dir, mustWork = FALSE))

  raw <- readabs::read_abs(cat_no = "6401.0", tables = 18, retain_files = TRUE)

  # Table 18 contains group, sub-group and expenditure class series. We want
  # only the 11 top-level groups. readabs's `series` column carries the full
  # ABS series label, typically of the form
  #   "Index Numbers ;  <Group name> ;  Australia ;"
  # We match on the group name appearing as the second ";"-delimited element.
  groups <- abs_cpi_groups()

  out <- raw |>
    dplyr::filter(
      .data$series_type == "Original",
      grepl("Index Numbers", .data$series, fixed = TRUE)
    ) |>
    dplyr::mutate(
      group_field = trimws(vapply(
        strsplit(.data$series, ";", fixed = TRUE),
        function(parts) if (length(parts) >= 2) parts[[2]] else NA_character_,
        character(1)
      ))
    ) |>
    dplyr::filter(.data$group_field %in% groups) |>
    # Communication and Education each appear twice in Table 18 — legacy
    # linked series carried alongside the current series. The pairs are
    # numerically identical apart from a handful of rounding-scale differences
    # in 1998–99 for Communication. Collapse to one row per (date, group).
    dplyr::group_by(.data$date, .data$group_field) |>
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(value = dplyr::if_else(is.nan(.data$value), NA_real_, .data$value)) |>
    dplyr::transmute(
      date = as.Date(.data$date),
      group = factor(.data$group_field, levels = groups),
      index = as.numeric(.data$value),
      source = "quarterly"
    ) |>
    dplyr::arrange(.data$group, .data$date) |>
    tibble::as_tibble()

  stopifnot(
    "Expected 11 distinct groups in quarterly data" =
      length(unique(out$group)) == length(groups)
  )
  out
}

#' Pull monthly group-level CPI index numbers (weighted average of 8 capital
#' cities) from ABS 6401.0 Table 2. Back-series begins April 2024.
#'
#' @return A tibble with columns `date`, `group`, `index`, `source` ("monthly").
fetch_cpi_monthly <- function(cache_dir = "data/raw") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(R_READABS_PATH = normalizePath(cache_dir, mustWork = FALSE))

  raw <- readabs::read_abs(cat_no = "6401.0", tables = 2, retain_files = TRUE)

  groups <- abs_cpi_groups()

  out <- raw |>
    dplyr::filter(
      .data$series_type == "Original",
      grepl("Index Numbers", .data$series, fixed = TRUE)
    ) |>
    dplyr::mutate(
      group_field = trimws(vapply(
        strsplit(.data$series, ";", fixed = TRUE),
        function(parts) if (length(parts) >= 2) parts[[2]] else NA_character_,
        character(1)
      ))
    ) |>
    dplyr::filter(.data$group_field %in% groups) |>
    dplyr::transmute(
      date = as.Date(.data$date),
      group = factor(.data$group_field, levels = groups),
      index = as.numeric(.data$value),
      source = "monthly"
    ) |>
    dplyr::arrange(.data$group, .data$date) |>
    tibble::as_tibble()

  stopifnot(
    "Expected 11 distinct groups in monthly data" =
      length(unique(out$group)) == length(groups)
  )
  out
}

#' Read the 2025 expenditure weights (weighted average of 8 capital cities)
#' from the cached ABS Weighting Pattern workbook, Table 1.
#'
#' @return A tibble with columns `group` (factor in the same order as
#'   `abs_cpi_groups()`) and `weight` (numeric, summing to 100).
fetch_weights <- function(path = "data/raw/cpi_weights_2025.xlsx") {
  if (!file.exists(path)) {
    stop("Weights workbook not found at ", path, ". Download from:\n",
         "  https://www.abs.gov.au/system/files/83cf3a8e846b6b3f91ae7ea99ef968ab/",
         "Consumer%20Price%20Index%20-%202025%20Weighting%20Pattern.xlsx")
  }

  # Table 1 has a staggered layout — label and weight both shift right as
  # you descend the hierarchy:
  #   col 1 = group label,        col 4 = group weight
  #   col 2 = sub-group label,    col 5 = sub-group weight
  #   col 3 = expenditure class,  col 6 = EC weight
  # Header rows occupy r1:r7.
  raw <- suppressMessages(
    readxl::read_excel(path, sheet = "Table 1", skip = 6, col_names = FALSE)
  )
  stopifnot("Expected at least 4 columns in Table 1" = ncol(raw) >= 4)
  names(raw)[1] <- "label"
  names(raw)[4] <- "group_weight"

  groups <- abs_cpi_groups()
  groups_upper <- toupper(groups)

  out <- raw |>
    dplyr::filter(.data$label %in% groups_upper) |>
    dplyr::transmute(
      group = factor(
        groups[match(.data$label, groups_upper)],
        levels = groups
      ),
      weight = as.numeric(.data$group_weight)
    ) |>
    dplyr::arrange(.data$group)

  stopifnot(
    "Expected 11 group weights" = nrow(out) == length(groups),
    "Group weights should sum to ~100" = abs(sum(out$weight) - 100) < 0.05
  )

  out
}
