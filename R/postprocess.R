# Post-processing helpers for MCT fits. All functions accept either a
# CmdStanMCMC fit (and pull draws via `$draws()`) or, where convenient, the
# already-extracted draws matrix/array — to keep them unit-testable without
# a live fit.

#' Extract a (D x T) matrix of draws for a vector-valued generated quantity.
extract_path_matrix <- function(fit, param) {
  d <- fit$draws(variables = param, format = "draws_matrix")
  m <- as.matrix(d)
  # Strip posterior attrs, retain column names like "trend[1]"
  attributes(m) <- list(dim = dim(m), dimnames = dimnames(m))
  m
}

#' Extract a (D x T x N) array for a matrix-valued generated quantity. Uses
#' Stan / R's column-major flattening: column k of the draws_matrix is
#' `sector_trend[((k-1) mod T) + 1, ((k-1) %/% T) + 1]`.
extract_matrix_array <- function(fit, param, T_, N) {
  m <- extract_path_matrix(fit, param)
  stopifnot(ncol(m) == T_ * N)
  array(as.numeric(m), dim = c(nrow(m), T_, N))
}

#' Summarise a (D x T) draws matrix on a date axis as a tidy tibble with
#' median + (lower, upper) credible bands.
summarise_path <- function(draws_matrix, dates,
                           probs = c(0.16, 0.5, 0.84)) {
  stopifnot(length(probs) == 3,
            probs[1] < probs[2], probs[2] < probs[3])
  stopifnot(ncol(draws_matrix) == length(dates))
  q <- apply(draws_matrix, 2, quantile, probs = probs, names = FALSE)
  tibble::tibble(
    date = dates,
    lower = q[1, ],
    median = q[2, ],
    upper = q[3, ]
  )
}

#' Trend path with credible bands (median + 68% by default).
trend_path <- function(fit, dates, probs = c(0.16, 0.5, 0.84)) {
  m <- extract_path_matrix(fit, "trend")
  summarise_path(m, dates, probs)
}

#' Common-share path with credible bands.
common_share_path <- function(fit, dates, probs = c(0.16, 0.5, 0.84)) {
  m <- extract_path_matrix(fit, "common_share")
  summarise_path(m, dates, probs)
}

#' Per-sector contributions w_i * s_{i,t} for stacked-area charts.
#'
#' @param fit CmdStanMCMC fit.
#' @param weights Numeric vector length N, normalised to sum 1.
#' @param dates Length-T Date vector.
#' @param groups Length-N character vector (sector display names).
#' @return Long tibble (date, group, lower, median, upper).
sector_contributions <- function(fit, weights, dates, groups,
                                 probs = c(0.16, 0.5, 0.84)) {
  T_ <- length(dates); N <- length(groups)
  stopifnot(length(weights) == N)
  w <- weights / sum(weights)
  s <- extract_matrix_array(fit, "sector_trend", T_, N)

  # For Variant B (RW common trend) we need to include each sector's
  # share of the common-trend movement so per-sector contributions sum
  # to the headline trend. Attribute the common-trend contribution to
  # sector i as w_i * lambda_i * c_t. We can detect Variant B from the
  # underlying fit's stored config (mct_gibbs_fit attaches it).
  variant_B <- !is.null(fit$config) &&
               isTRUE(fit$config$variant == "B")
  if (variant_B) {
    lambda <- extract_path_matrix(fit, "lambda")  # n_draw x N
    c_path <- extract_path_matrix(fit, "c")       # n_draw x T_
    # Per-draw, per-(t, i): w_i * lambda_i * c_t. `s` (and the result)
    # is shape (n_draw, T_, N); build com_arr to match.
    n_draw <- nrow(c_path)
    com_arr <- array(0, dim = c(n_draw, T_, N))
    for (i in seq_len(N)) {
      # lambda[, i] is length n_draw; replicate across T_:
      com_arr[, , i] <- w[i] * (lambda[, i] * c_path)
    }
    contribs <- sweep(s, 3, w, FUN = "*") + com_arr
  } else {
    # Variant A: common factor is a cycle, NOT part of trend. Sector
    # contribution is just w_i * s_{i,t}.
    contribs <- sweep(s, 3, w, FUN = "*")
  }

  dplyr::bind_rows(lapply(seq_len(N), function(i) {
    tib <- summarise_path(contribs[, , i], dates, probs)
    tib$group <- groups[i]
    tib
  })) |>
    dplyr::mutate(group = factor(.data$group, levels = groups)) |>
    dplyr::select("date", "group", "lower", "median", "upper")
}

#' Compute the weighted-average of per-series demean offsets — the
#' constant to add back to a demeaned `trend_path()` output to recover
#' raw % inflation (so the trend line is on the same axis as headline
#' or trimmed-mean CPI YoY series).
#'
#' @param infl_demeaned The same tibble passed to `build_stan_data...()`.
#'   Must have columns `group`, `series_mean` (added by `demean_series()`).
#' @param weights Tibble with `group` and `weight` (the same one passed
#'   to the stan_data builder; normalised internally).
#' @return Scalar numeric — the pp constant to add to demeaned trend.
demean_offset <- function(infl_demeaned, weights) {
  stopifnot("series_mean" %in% names(infl_demeaned))
  means <- infl_demeaned |>
    dplyr::distinct(.data$group, .data$series_mean) |>
    dplyr::filter(!is.na(.data$series_mean))
  w_lookup <- setNames(weights$weight, as.character(weights$group))
  w <- w_lookup[as.character(means$group)] / sum(weights$weight)
  sum(w * means$series_mean)
}

#' Decompose the (demeaned) Variant B trend into per-bucket common +
#' sector-specific components, matching the NY Fed MCT decomposition
#' chart's "Sector-Specific and Common Components" view.
#'
#' For each macro bucket (Goods / Services ex. housing / Housing):
#'   common[t]   = sum_{i in bucket} w_i * lambda_i * c_t
#'   specific[t] = sum_{i in bucket} w_i * s_{i,t}
#'
#' Per-draw, then median + 68% band.
#'
#' @param fit mct_gibbs_fit wrapper (Variant B).
#' @param weights Numeric vector of length N (sub-group weights).
#' @param dates Length-T Date vector.
#' @param sub_groups Length-N character vector — sub-group labels.
#' @param bucket_map Named character (sub-group -> bucket).
#' @return Long tibble: date, bucket, component (common|specific),
#'   lower, median, upper.
decompose_trend_by_bucket <- function(fit, weights, dates, sub_groups,
                                      bucket_map,
                                      probs = c(0.16, 0.5, 0.84)) {
  T_ <- length(dates); N <- length(sub_groups)
  stopifnot(length(weights) == N)
  w <- weights / sum(weights)

  s_arr  <- extract_matrix_array(fit, "sector_trend", T_, N)  # draws x T x N
  lambda <- extract_path_matrix(fit, "lambda")                # draws x N
  c_path <- extract_path_matrix(fit, "c")                     # draws x T

  buckets <- unique(unname(bucket_map[sub_groups]))
  out <- list()
  for (b in buckets) {
    in_bucket <- which(bucket_map[sub_groups] == b)
    # Common contribution: sum_{i in bucket} w_i * lambda_i * c_t
    common_per_draw <- matrix(0, nrow(c_path), T_)
    for (i in in_bucket) {
      common_per_draw <- common_per_draw +
        w[i] * (lambda[, i] * c_path)
    }
    # Specific contribution: sum_{i in bucket} w_i * s_{i,t}
    specific_per_draw <- matrix(0, nrow(c_path), T_)
    for (i in in_bucket) {
      specific_per_draw <- specific_per_draw + w[i] * s_arr[, , i]
    }
    out[[length(out) + 1L]] <- summarise_path(common_per_draw,
                                              dates, probs) |>
      dplyr::mutate(bucket = b, component = "Common")
    out[[length(out) + 1L]] <- summarise_path(specific_per_draw,
                                              dates, probs) |>
      dplyr::mutate(bucket = b, component = "Sector-specific")
  }
  bucket_order <- c("Goods", "Services ex. housing", "Housing")
  dplyr::bind_rows(out) |>
    dplyr::mutate(bucket = factor(.data$bucket, levels = bucket_order),
                  component = factor(.data$component,
                                     levels = c("Common", "Sector-specific"))) |>
    dplyr::select("date", "bucket", "component",
                  "lower", "median", "upper")
}

#' Roll up a sub-group `sector_contributions()` tibble to its parent
#' groups (the 11 standard ABS CPI groups). Sums the medians per
#' parent-group; the upper/lower bands of the summed series are NOT the
#' simple sum of the per-sub-group bands (medians don't sum exactly), so
#' the returned `lower`/`upper` for parent-groups are the simple sum and
#' should be treated as illustrative — they aren't intended for tight
#' uncertainty quantification. For the stacked-area chart, only `median`
#' is used.
#'
#' @param contribs Long tibble from `sector_contributions()`. The `group`
#'   factor levels are the sub-group names.
#' @param parent_map Named character vector (sub-group -> parent), e.g.
#'   `abs_subgroup_parent()`.
#' @param parent_order Optional character vector of parent group names in
#'   the order to factor them. Defaults to `unique(parent_map)`.
aggregate_contribs_by_parent <- function(contribs, parent_map,
                                         parent_order = NULL) {
  parent_levels <- if (is.null(parent_order)) {
    unique(unname(parent_map))
  } else {
    parent_order
  }
  contribs |>
    dplyr::mutate(parent = factor(parent_map[as.character(.data$group)],
                                  levels = parent_levels)) |>
    dplyr::group_by(.data$date, .data$parent) |>
    dplyr::summarise(
      lower  = sum(.data$lower),
      median = sum(.data$median),
      upper  = sum(.data$upper),
      .groups = "drop"
    ) |>
    dplyr::rename(group = "parent") |>
    dplyr::select("date", "group", "lower", "median", "upper")
}

#' Compare trend estimates from this month's fit vs last month's fit at the
#' dates they both cover. Powers the dashboard's "revisions" panel.
#'
#' @param current_fit,prior_fit Two CmdStanMCMC fits of the same variant.
#' @param current_dates,prior_dates The Date vectors used by each fit
#'   (i.e., `stan_data$dates`).
#' @return Tibble with one row per overlapping date and columns
#'   (date, curr_median, curr_lower, curr_upper, prior_median, prior_lower,
#'    prior_upper, diff) where `diff = curr_median - prior_median`.
posterior_revisions <- function(current_fit, prior_fit,
                                current_dates, prior_dates,
                                probs = c(0.16, 0.5, 0.84)) {
  common_dates <- intersect(as.character(current_dates),
                            as.character(prior_dates))
  if (length(common_dates) == 0) {
    return(tibble::tibble(
      date = as.Date(character(0)),
      curr_median = numeric(0), curr_lower = numeric(0), curr_upper = numeric(0),
      prior_median = numeric(0), prior_lower = numeric(0), prior_upper = numeric(0),
      diff = numeric(0)
    ))
  }
  curr_idx <- match(common_dates, as.character(current_dates))
  prev_idx <- match(common_dates, as.character(prior_dates))

  curr_m <- extract_path_matrix(current_fit, "trend")[, curr_idx, drop = FALSE]
  prev_m <- extract_path_matrix(prior_fit, "trend")[, prev_idx, drop = FALSE]
  dates_d <- as.Date(common_dates)

  cs <- summarise_path(curr_m, dates_d, probs) |>
    dplyr::rename(curr_lower = "lower", curr_median = "median",
                  curr_upper = "upper")
  ps <- summarise_path(prev_m, dates_d, probs) |>
    dplyr::rename(prior_lower = "lower", prior_median = "median",
                  prior_upper = "upper")

  dplyr::left_join(cs, ps, by = "date") |>
    dplyr::mutate(diff = .data$curr_median - .data$prior_median)
}

#' Sector-by-sector attribution of the trend revision between two vintages.
#'
#' Decomposes the difference in median trend at overlapping dates into a sum
#' of per-sector contributions: which sector's contribution to the headline
#' trend moved most between this month's vintage and last month's? This is
#' the MCMC-native analogue of the NY Fed MCT's Kalman-update attribution in
#' `decompose_update.m` — same editorial question ("what drove this month's
#' revision?"), different mechanics. Theirs runs a fresh Kalman filter with
#' and without the latest data point at posterior-median parameters; ours
#' compares the per-sector posterior-median contributions across two fits.
#'
#' Both fits should be the same variant (or at least produce a
#' `sector_trend[T, N]` generated quantity).
#'
#' @param current_fit,prior_fit CmdStanMCMC fits.
#' @param current_dates,prior_dates Date vectors used by each fit.
#' @param weights Numeric vector of length N (will be normalised to sum 1).
#' @param groups Length-N character vector of sector display names.
#' @return Long tibble with one row per overlapping (date, group):
#'   `date, group, curr_contrib, prior_contrib, delta_contrib`. Summing
#'   `delta_contrib` across groups at a given date recovers the change in
#'   the headline trend median between the two vintages (up to within-fit
#'   Monte-Carlo error, since each median is from a finite sample of draws).
decompose_revisions <- function(current_fit, prior_fit,
                                current_dates, prior_dates,
                                weights, groups) {
  curr <- sector_contributions(current_fit, weights, current_dates, groups,
                               probs = c(0.16, 0.5, 0.84))
  prev <- sector_contributions(prior_fit, weights, prior_dates, groups,
                               probs = c(0.16, 0.5, 0.84))

  curr <- dplyr::select(curr, "date", "group",
                        curr_contrib = "median")
  prev <- dplyr::select(prev, "date", "group",
                        prior_contrib = "median")

  dplyr::inner_join(curr, prev, by = c("date", "group")) |>
    dplyr::mutate(delta_contrib = .data$curr_contrib - .data$prior_contrib)
}

#' Read a slim parquet of fit output (produced by export_fit_parquet) and
#' rebuild the path summaries the dashboard needs without re-loading the
#' full fit. Handy for downstream consumers.
read_slim_parquet <- function(path, dates, groups,
                              probs = c(0.16, 0.5, 0.84)) {
  d <- arrow::read_parquet(path)
  d <- as.data.frame(d)
  m <- as.matrix(d)
  T_ <- length(dates); N <- length(groups)

  cols <- function(prefix, k) {
    grep(paste0("^", prefix, "\\["), colnames(m), value = TRUE)[seq_len(k)]
  }

  trend_m <- m[, grep("^trend\\[", colnames(m))]
  cs_m    <- m[, grep("^common_share\\[", colnames(m))]
  ct_m    <- m[, grep("^common_transitory\\[", colnames(m))]
  st_cols <- grep("^sector_trend\\[", colnames(m), value = TRUE)
  st_m    <- m[, st_cols, drop = FALSE]
  st_arr  <- array(as.numeric(st_m), dim = c(nrow(st_m), T_, N))

  list(
    trend = summarise_path(trend_m, dates, probs),
    common_share = summarise_path(cs_m, dates, probs),
    common_transitory = summarise_path(ct_m, dates, probs),
    sector_trend_array = st_arr
  )
}
