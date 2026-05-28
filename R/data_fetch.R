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

#' Map each ABS CPI sub-group to a 3-way macro bucket
#' (Goods / Services ex. housing / Housing) — matches the NY Fed MCT
#' decomposition's broad sector aggregation.
#'
#' Splits chosen by dominant character of each sub-group:
#'   - Housing bucket: 4 Housing sub-groups (rents, new dwellings,
#'     other housing, utilities).
#'   - Services ex. housing bucket: services-dominated sub-groups
#'     (Communication, Education, Insurance, Financial services,
#'     Medical services, Meals out, Holiday travel, Domestic services,
#'     Urban transport fares, plus services within Recreation).
#'   - Goods bucket: everything else (food items, alcohol, tobacco,
#'     garments, furniture, motor vehicles + fuel, etc.).
#' This mirrors ABS's analytical goods/services split at sub-group
#' granularity; it's approximate because some sub-groups have a mix.
abs_subgroup_to_macrobucket <- function() {
  c(
    # Housing
    "Rents"                                              = "Housing",
    "New dwelling purchase by owner-occupiers"           = "Housing",
    "Other housing"                                      = "Housing",
    "Utilities"                                          = "Housing",
    # Services ex. housing
    "Meals out and take away foods"                      = "Services ex. housing",
    "Domestic and household services"                    = "Services ex. housing",
    "Medical, dental and hospital services"              = "Services ex. housing",
    "Urban transport fares"                              = "Services ex. housing",
    "Communication"                                      = "Services ex. housing",
    "Holiday travel and accommodation"                   = "Services ex. housing",
    "Other recreation, sport and culture"                = "Services ex. housing",
    "Education"                                          = "Services ex. housing",
    "Insurance"                                          = "Services ex. housing",
    "Financial services"                                 = "Services ex. housing",
    # Goods
    "Bread and cereal products"                          = "Goods",
    "Meat and seafoods"                                  = "Goods",
    "Dairy and related products"                         = "Goods",
    "Fruit and vegetables"                               = "Goods",
    "Food products n.e.c."                               = "Goods",
    "Non-alcoholic beverages"                            = "Goods",
    "Alcoholic beverages"                                = "Goods",
    "Tobacco"                                            = "Goods",
    "Garments"                                           = "Goods",
    "Footwear"                                           = "Goods",
    "Accessories and clothing services"                  = "Goods",
    "Furniture and furnishings"                          = "Goods",
    "Household textiles"                                 = "Goods",
    "Household appliances, utensils and tools"           = "Goods",
    "Non-durable household products"                     = "Goods",
    "Medical products, appliances and equipment"         = "Goods",
    "Private motoring"                                   = "Goods",
    "Audio, visual and computing equipment and services" = "Goods",
    "Newspapers, books and stationery"                   = "Goods"
  )
}

#' Map each ABS CPI sub-group to its parent group. Returns a named
#' character vector: names = sub-group labels (matching
#' `abs_cpi_subgroups()`), values = parent group labels (matching
#' `abs_cpi_groups()`). Used to aggregate the 33 sub-group sector
#' contributions back to the 11 parent groups for the dashboard's
#' stacked-area chart.
abs_subgroup_parent <- function() {
  c(
    # Housing (4)
    "Rents"                                              = "Housing",
    "New dwelling purchase by owner-occupiers"           = "Housing",
    "Other housing"                                      = "Housing",
    "Utilities"                                          = "Housing",
    # Food and non-alcoholic beverages (7)
    "Bread and cereal products"                          = "Food and non-alcoholic beverages",
    "Meat and seafoods"                                  = "Food and non-alcoholic beverages",
    "Dairy and related products"                         = "Food and non-alcoholic beverages",
    "Fruit and vegetables"                               = "Food and non-alcoholic beverages",
    "Food products n.e.c."                               = "Food and non-alcoholic beverages",
    "Non-alcoholic beverages"                            = "Food and non-alcoholic beverages",
    "Meals out and take away foods"                      = "Food and non-alcoholic beverages",
    # Alcohol and tobacco (2)
    "Alcoholic beverages"                                = "Alcohol and tobacco",
    "Tobacco"                                            = "Alcohol and tobacco",
    # Clothing and footwear (3)
    "Garments"                                           = "Clothing and footwear",
    "Footwear"                                           = "Clothing and footwear",
    "Accessories and clothing services"                  = "Clothing and footwear",
    # Furnishings (5)
    "Furniture and furnishings"                          = "Furnishings, household equipment and services",
    "Household textiles"                                 = "Furnishings, household equipment and services",
    "Household appliances, utensils and tools"           = "Furnishings, household equipment and services",
    "Non-durable household products"                     = "Furnishings, household equipment and services",
    "Domestic and household services"                    = "Furnishings, household equipment and services",
    # Health (2)
    "Medical products, appliances and equipment"         = "Health",
    "Medical, dental and hospital services"              = "Health",
    # Transport (2)
    "Private motoring"                                   = "Transport",
    "Urban transport fares"                              = "Transport",
    # Communication (1)
    "Communication"                                      = "Communication",
    # Recreation and culture (4)
    "Audio, visual and computing equipment and services" = "Recreation and culture",
    "Newspapers, books and stationery"                   = "Recreation and culture",
    "Holiday travel and accommodation"                   = "Recreation and culture",
    "Other recreation, sport and culture"                = "Recreation and culture",
    # Education (1)
    "Education"                                          = "Education",
    # Insurance and financial services (2)
    "Insurance"                                          = "Insurance and financial services",
    "Financial services"                                 = "Insurance and financial services"
  )
}

#' The 33 ABS CPI sub-groups, ordered with "Rents" first so it can be used
#' as the identification reference. "Rents" is observed continuously back
#' to 1972 (the largest "New dwelling purchase by owner-occupiers" only
#' starts 1998-06 in Table 18, which would break the `lambda_ref,t = 1`
#' identification on the 1993-onward sample). Rents is the second-largest
#' Housing sub-group (w=6.60 in the 2025 weighting pattern), stable, and
#' housing-adjacent for continuity with the previous Housing-group
#' reference. Other sub-groups follow in the same group order as
#' `abs_cpi_groups()`.
abs_cpi_subgroups <- function() {
  c(
    # Housing (reference first)
    "Rents",
    "New dwelling purchase by owner-occupiers",
    "Other housing",
    "Utilities",
    # Food and non-alcoholic beverages
    "Bread and cereal products",
    "Meat and seafoods",
    "Dairy and related products",
    "Fruit and vegetables",
    "Food products n.e.c.",
    "Non-alcoholic beverages",
    "Meals out and take away foods",
    # Alcohol and tobacco
    "Alcoholic beverages",
    "Tobacco",
    # Clothing and footwear
    "Garments",
    "Footwear",
    "Accessories and clothing services",
    # Furnishings, household equipment and services
    "Furniture and furnishings",
    "Household textiles",
    "Household appliances, utensils and tools",
    "Non-durable household products",
    "Domestic and household services",
    # Health
    "Medical products, appliances and equipment",
    "Medical, dental and hospital services",
    # Transport
    "Private motoring",
    "Urban transport fares",
    # Communication
    "Communication",
    # Recreation and culture
    "Audio, visual and computing equipment and services",
    "Newspapers, books and stationery",
    "Holiday travel and accommodation",
    "Other recreation, sport and culture",
    # Education
    "Education",
    # Insurance and financial services
    "Insurance",
    "Financial services"
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

#' Pull the All-groups CPI YoY series from ABS 6401.0 Table 18 (the
#' first row of the quarterly publication). Returns the year-on-year
#' percentage change in the headline CPI — the most-cited Australian
#' inflation measure. Used as an overlay on the dashboard's trend chart
#' alongside the model's MCT.
#'
#' @return Tibble with `date` (quarter-end date) and `yoy` (numeric, %).
fetch_cpi_headline_yoy <- function(cache_dir = "data/raw") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(R_READABS_PATH = normalizePath(cache_dir, mustWork = FALSE))

  raw <- readabs::read_abs(cat_no = "6401.0", tables = 18,
                           retain_files = TRUE)

  all_groups <- raw |>
    dplyr::filter(
      .data$series_type == "Original",
      grepl("Index Numbers", .data$series, fixed = TRUE),
      grepl("All groups CPI ;  Australia", .data$series, fixed = TRUE)
    ) |>
    dplyr::group_by(.data$date) |>
    dplyr::summarise(index = mean(.data$value, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(.data$date)

  all_groups |>
    dplyr::mutate(
      date = as.Date(.data$date),
      yoy  = 100 * (.data$index / dplyr::lag(.data$index, 4L) - 1)
    ) |>
    dplyr::filter(!is.na(.data$yoy)) |>
    dplyr::select("date", "yoy") |>
    tibble::as_tibble()
}

#' Pull the Trimmed Mean CPI YoY series from ABS 6401.0 Appendix 1a
#' (Analytical Series). This is the RBA's preferred core-inflation
#' measure (year-on-year % change in the trimmed-mean CPI index).
#'
#' The appendix isn't part of the auto-fetched tables 1-18 — it's a
#' separate spreadsheet posted with each quarterly release. We point
#' at a pinned URL; on the next release the URL changes (`apr-2026`
#' segment) so the cron will need to be pointed at the new vintage.
#'
#' @param path Local cache path. If missing, downloads from `url`.
#' @param url  Source URL (default points to the Apr-2026 release).
#' @return Tibble with `date` (quarter-end), `yoy` (numeric, %).
fetch_cpi_trimmed_mean_yoy <- function(
    path = "data/raw/6401.0/64010appendix1a.xlsx",
    url  = paste0("https://www.abs.gov.au/statistics/economy/",
                  "price-indexes-and-inflation/consumer-price-index-australia/",
                  "apr-2026/64010Appendix1a.xlsx")) {
  if (!file.exists(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    utils::download.file(url, destfile = path, mode = "wb", quiet = TRUE)
  }
  raw <- suppressMessages(
    readxl::read_excel(path, sheet = "Data1", col_names = TRUE)
  )
  series_name <- paste0(
    "Percentage Change from Corresponding Quarter of Previous Year",
    " ;  Trimmed Mean ;  Australia ;"
  )
  stopifnot("Trimmed Mean YoY column not found" =
              series_name %in% names(raw))

  # Header rows in ABS spreadsheets: first ~9 rows are metadata. The
  # first column has the date in the data rows; we detect the data
  # start by the first row whose 1st column parses to a Date.
  first_col <- raw[[1L]]
  date_parsed <- suppressWarnings(as.Date(as.numeric(first_col),
                                          origin = "1899-12-30"))
  start_row <- which(!is.na(date_parsed) & date_parsed > as.Date("1900-01-01"))[1L]
  if (is.na(start_row)) {
    stop("Could not locate data start row in appendix 1a")
  }
  data <- raw[start_row:nrow(raw), ]
  dates <- as.Date(as.numeric(data[[1L]]), origin = "1899-12-30")
  yoy   <- suppressWarnings(as.numeric(data[[series_name]]))
  tibble::tibble(date = dates, yoy = yoy) |>
    dplyr::filter(!is.na(.data$date), !is.na(.data$yoy))
}

#' Read the full 3-level (group / sub-group / expenditure-class) hierarchy
#' from the ABS Weighting Pattern workbook. Returns one row per EC with
#' parent group, parent sub-group, EC label and EC weight (sums to ~100).
#'
#' Used as the EC→sub-group mapping for aggregating seasonally-adjusted
#' EC series in Appendix 1a back up to our 33 sub-groups.
fetch_weights_ec_hierarchy <- function(
    path = "data/raw/cpi_weights_2025.xlsx") {
  if (!file.exists(path)) {
    stop("Weights workbook not found at ", path)
  }
  raw <- suppressMessages(
    readxl::read_excel(path, sheet = "Table 1", skip = 6,
                       col_names = FALSE)
  )
  names(raw)[1:6] <- c("group", "subgroup", "ec",
                       "w_group", "w_subgroup", "w_ec")
  groups_upper <- toupper(abs_cpi_groups())
  cur_group <- NA_character_
  cur_sg    <- NA_character_
  out <- list()
  for (i in seq_len(nrow(raw))) {
    if (!is.na(raw$group[i]) && raw$group[i] %in% groups_upper) {
      cur_group <- abs_cpi_groups()[match(raw$group[i], groups_upper)]
    }
    if (!is.na(raw$subgroup[i]) && !is.na(raw$w_subgroup[i])) {
      cur_sg <- raw$subgroup[i]
    }
    if (!is.na(raw$ec[i]) && !is.na(raw$w_ec[i])) {
      out[[length(out) + 1L]] <- data.frame(
        group    = cur_group,
        subgroup = cur_sg,
        ec       = raw$ec[i],
        weight   = as.numeric(raw$w_ec[i])
      )
    }
  }
  tibble::as_tibble(do.call(rbind, out))
}

#' Pull the seasonally-adjusted quarterly CPI sub-group inflation series.
#'
#' Sources expenditure-class SA indices from ABS 6401.0 Appendix 1a
#' (90 ECs, SA back to 1972-09) and aggregates them up to our 33
#' sub-groups using the EC weights from the 2025 Weighting Pattern.
#'
#' This replaces `fetch_cpi_quarterly_subgroup()` for the production
#' pipeline — the SA series strips out seasonal patterns (e.g.
#' Education's annual Q1 fee re-indexing) that the model would
#' otherwise mistake for outliers.
#'
#' @param path Local cache path for the appendix workbook.
#' @param url  Source URL (current ABS release).
#' @return Tibble with `date, group, index, source = "quarterly"`,
#'   where `group` carries the 33 sub-group labels (matching
#'   `abs_cpi_subgroups()`).
fetch_cpi_quarterly_subgroup_sa <- function(
    path = "data/raw/6401.0/64010appendix1a.xlsx",
    url  = paste0("https://www.abs.gov.au/statistics/economy/",
                  "price-indexes-and-inflation/consumer-price-index-australia/",
                  "apr-2026/64010Appendix1a.xlsx")) {
  if (!file.exists(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    utils::download.file(url, destfile = path, mode = "wb", quiet = TRUE)
  }
  raw <- suppressMessages(
    readxl::read_excel(path, sheet = "Data1", col_names = TRUE)
  )

  # Data starts at the first row whose 1st column parses to a Date.
  first_col <- raw[[1L]]
  dp <- suppressWarnings(as.Date(as.numeric(first_col),
                                 origin = "1899-12-30"))
  start_row <- which(!is.na(dp) & dp > as.Date("1900-01-01"))[1L]
  if (is.na(start_row)) {
    stop("Could not locate data start row in Appendix 1a")
  }
  d <- raw[start_row:nrow(raw), ]
  dates <- as.Date(as.numeric(d[[1L]]), origin = "1899-12-30")

  # All "Index Numbers" columns are SA EC indices (plus All groups, TM, WM).
  idx_cols <- grep("Index Numbers", names(raw), value = TRUE)
  clean_ec <- function(s) {
    s <- sub("^Index Numbers ;\\s+", "", s)
    s <- sub("\\s+;\\s+Australia\\s+;\\s*$", "", s)
    trimws(s)
  }
  ec_labels <- vapply(idx_cols, clean_ec, character(1))
  # Drop the 3 aggregates (All groups SA, Trimmed Mean, Weighted Median)
  drop <- c("All groups CPI, seasonally adjusted",
            "Trimmed Mean", "Weighted Median")
  keep <- !ec_labels %in% drop
  idx_cols <- idx_cols[keep]
  ec_labels <- ec_labels[keep]

  # Load the EC→sub-group hierarchy + weights
  hier <- fetch_weights_ec_hierarchy()

  # Reconcile naming differences between Appendix 1a EC labels and the
  # weighting-workbook EC labels. The hierarchy table has labels like
  # "Garments for men" / "Garments for women" / "Garments for infants
  # and children" — Appendix 1a uses the same labels. Any unmatched
  # EC label is dropped with a warning.
  unmatched <- setdiff(ec_labels, hier$ec)
  if (length(unmatched) > 0) {
    warning("Dropping ", length(unmatched),
            " Appendix 1a ECs not found in weights hierarchy: ",
            paste(head(unmatched, 5), collapse = ", "),
            if (length(unmatched) > 5) ", ..." else "")
  }

  # Build long-format SA inflation, then aggregate to sub-groups.
  sg_to_keep <- abs_cpi_subgroups()
  ec_to_sg <- setNames(hier$subgroup, hier$ec)
  ec_to_w  <- setNames(hier$weight,   hier$ec)

  ec_long <- dplyr::bind_rows(lapply(seq_along(idx_cols), function(j) {
    ec    <- ec_labels[j]
    sg    <- ec_to_sg[ec]
    w_ec  <- ec_to_w[ec]
    if (is.na(sg) || !(sg %in% sg_to_keep)) return(NULL)
    idx <- as.numeric(d[[idx_cols[j]]])
    data.frame(
      date     = dates,
      ec       = ec,
      subgroup = sg,
      w_ec     = w_ec,
      index    = idx
    )
  }))

  # Compute per-EC annualised q/q inflation, then weighted-aggregate
  # by EC weight within each sub-group → sub-group inflation rate.
  ec_long <- ec_long |>
    dplyr::arrange(.data$ec, .data$date) |>
    dplyr::group_by(.data$ec) |>
    dplyr::mutate(
      qoq_ann = 4 * (log(.data$index) - dplyr::lag(log(.data$index))) * 100
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$qoq_ann), is.finite(.data$qoq_ann))

  sg_long <- ec_long |>
    dplyr::group_by(.data$date, .data$subgroup) |>
    dplyr::summarise(
      qoq_ann_agg = sum(.data$qoq_ann * .data$w_ec, na.rm = TRUE) /
                    sum(.data$w_ec[!is.na(.data$qoq_ann)]),
      .groups = "drop"
    )

  # The model expects an INDEX (not inflation directly) and the
  # downstream build_inflation_series re-derives inflation from it.
  # Reconstruct a synthetic SA index per sub-group from the aggregated
  # inflation (start at 100 in the first available period, cumulate).
  sg_long |>
    dplyr::arrange(.data$subgroup, .data$date) |>
    dplyr::group_by(.data$subgroup) |>
    dplyr::mutate(
      qoq_log = .data$qoq_ann_agg / 4 / 100,         # back to log q/q
      log_idx = cumsum(.data$qoq_log),
      index   = 100 * exp(.data$log_idx)
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      date   = as.Date(.data$date),
      group  = factor(.data$subgroup, levels = sg_to_keep),
      index  = .data$index,
      source = "quarterly"
    )
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

#' Pull quarterly sub-group CPI index numbers from ABS 6401.0 Table 18.
#'
#' Same source as `fetch_cpi_quarterly()` but extracts the 33 sub-group
#' series (one level below groups) instead of the 11 top-level groups.
#' Sub-group breakdowns are quarterly only — the ABS Monthly CPI Table 2
#' publishes group-level data only.
#'
#' @return Tibble with columns `date, group, index, source = "quarterly"`,
#'   where the `group` column is renamed for compatibility with downstream
#'   code that expects a "group" factor — entries are the 33 sub-group
#'   labels from `abs_cpi_subgroups()`. The actual parent group is
#'   available as an attribute.
fetch_cpi_quarterly_subgroup <- function(cache_dir = "data/raw") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(R_READABS_PATH = normalizePath(cache_dir, mustWork = FALSE))

  raw <- readabs::read_abs(cat_no = "6401.0", tables = 18,
                           retain_files = TRUE)

  subgroups <- abs_cpi_subgroups()

  out <- raw |>
    dplyr::filter(
      .data$series_type == "Original",
      grepl("Index Numbers", .data$series, fixed = TRUE)
    ) |>
    dplyr::mutate(
      label = trimws(vapply(
        strsplit(.data$series, ";", fixed = TRUE),
        function(parts) if (length(parts) >= 2) parts[[2]] else NA_character_,
        character(1)
      ))
    ) |>
    dplyr::filter(.data$label %in% subgroups) |>
    # Communication and Education each appear at both group AND sub-group
    # level with identical numbers (the sub-group IS the only sub-group
    # under that parent). The label-matching here picks them up once; the
    # subsequent group-by collapse any rare duplicates from legacy linked
    # series in Table 18.
    dplyr::group_by(.data$date, .data$label) |>
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(value = dplyr::if_else(is.nan(.data$value),
                                         NA_real_, .data$value)) |>
    dplyr::transmute(
      date = as.Date(.data$date),
      group = factor(.data$label, levels = subgroups),
      index = as.numeric(.data$value),
      source = "quarterly"
    ) |>
    dplyr::arrange(.data$group, .data$date) |>
    tibble::as_tibble()

  found <- unique(as.character(out$group))
  missing <- setdiff(subgroups, found)
  if (length(missing) > 0L) {
    warning("Sub-groups missing from Table 18: ",
            paste(missing, collapse = ", "))
  }

  out
}

#' Read the 2025 sub-group weights from the ABS Weighting Pattern workbook.
#'
#' Walks the Table 1 hierarchy to extract the 33 sub-group weights, in the
#' order of `abs_cpi_subgroups()`.
#'
#' @return Tibble with columns `group` (factor of sub-group labels in
#'   abs_cpi_subgroups() order — column named `group` for compatibility
#'   with downstream code) and `weight` (numeric, summing to ~100).
fetch_weights_subgroup <- function(path = "data/raw/cpi_weights_2025.xlsx") {
  if (!file.exists(path)) {
    stop("Weights workbook not found at ", path, ". Download from:\n",
         "  https://www.abs.gov.au/system/files/83cf3a8e846b6b3f91ae7ea99ef968ab/",
         "Consumer%20Price%20Index%20-%202025%20Weighting%20Pattern.xlsx")
  }

  raw <- suppressMessages(
    readxl::read_excel(path, sheet = "Table 1", skip = 6, col_names = FALSE)
  )
  stopifnot("Expected at least 5 columns in Table 1" = ncol(raw) >= 5)
  names(raw)[2] <- "subgroup_label"
  names(raw)[5] <- "subgroup_weight"

  subgroups <- abs_cpi_subgroups()

  out <- raw |>
    dplyr::filter(!is.na(.data$subgroup_label),
                  !is.na(.data$subgroup_weight),
                  .data$subgroup_label %in% subgroups) |>
    dplyr::transmute(
      group = factor(.data$subgroup_label, levels = subgroups),
      weight = as.numeric(.data$subgroup_weight)
    ) |>
    dplyr::arrange(.data$group)

  stopifnot(
    "Expected 33 sub-group weights" = nrow(out) == length(subgroups),
    "Sub-group weights should sum to ~100" =
      abs(sum(out$weight) - 100) < 0.05
  )

  out
}
