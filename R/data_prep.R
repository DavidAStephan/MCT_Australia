# Data preparation: from raw CPI index numbers to the (T x N) matrices
# that the Stan model consumes.
#
# Production sample starts at 1993-01 (RBA inflation-targeting era — see
# Step 16/17 sensitivity study). The brief originally specified 1990-01;
# the regime change of the 1990-92 disinflation was inflating the rho
# posterior of Variant A. Each cell (t, n) is one of three states encoded
# by `obs_type`:
#
#   obs_type == 0  -> unobserved (Stan integrates over the latent inflation)
#   obs_type == 1  -> monthly observation (Apr 2024 onward)
#   obs_type == 2  -> quarterly observation placed at the last month of the
#                     quarter; equals the simple average of latent monthly
#                     inflations in months t-2, t-1, t

# Demean window matches the production sample start (1993-01) — see
# Step 16/17 sensitivity study (May 2026): switching to 1993-onward
# cuts the early-90s disinflation regime change, and using a
# 1993-2019 demean window keeps the demeaned series in-sample with
# minimal pandemic contamination.
DEMEAN_WINDOW <- c(1993, 2019)

#' Convert a tidy (date, group, index, source) frame to annualised
#' period-on-period log-differenced inflation.
#'
#' Annualisation factor: 4 for quarterly, 12 for monthly. Inflation is
#' computed within (group, source); the first period of each series is NA.
build_inflation_series <- function(df) {
  stopifnot(all(c("date", "group", "index", "source") %in% names(df)))

  df |>
    dplyr::arrange(.data$group, .data$source, .data$date) |>
    dplyr::group_by(.data$group, .data$source) |>
    dplyr::mutate(
      ann_factor = dplyr::if_else(.data$source == "monthly", 12, 4),
      inflation = .data$ann_factor *
        (log(.data$index) - dplyr::lag(log(.data$index))) * 100
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-"ann_factor") |>
    dplyr::filter(!is.na(.data$inflation))
}

#' Subtract the 2000–2019 mean of *quarterly* inflation from every observation
#' of the same group, regardless of frequency. Quarterly and monthly
#' annualised inflations are on the same scale, so the same mean applies.
#'
#' Returns the same frame with `inflation` replaced by the demeaned series and
#' a new column `series_mean` recording the subtracted value.
demean_series <- function(infl_df,
                          window = DEMEAN_WINDOW) {
  means <- infl_df |>
    dplyr::filter(
      .data$source == "quarterly",
      lubridate::year(.data$date) >= window[1],
      lubridate::year(.data$date) <= window[2]
    ) |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(series_mean = mean(.data$inflation, na.rm = TRUE),
                     .groups = "drop")

  infl_df |>
    dplyr::left_join(means, by = "group") |>
    dplyr::mutate(inflation = .data$inflation - .data$series_mean)
}

#' Build the (T x N) matrices and flat observation arrays the Stan model
#' consumes.
#'
#' Combines quarterly and monthly inflation onto a single monthly latent
#' clock running from `start = "1993-01-01"` (RBA inflation-targeting era)
#' to the most recent observed month. For t >= 2024-04-01 the monthly
#' observation is used (obs_type=1);
#' for earlier dates that coincide with a quarter-end month, the quarterly
#' observation is placed there (obs_type=2); other cells are unobserved
#' (obs_type=0).
#'
#' Quarterly observations require a 3-month look-back window (t >= 3) so the
#' model can average π[t-2], π[t-1], π[t]. Any quarterly observation falling
#' on t=1 or t=2 of the latent clock is dropped from the Stan input.
#'
#' @param infl_demeaned Output of `demean_series()`.
#' @param weights A tibble with columns `group` (factor matching the levels of
#'   `infl_demeaned$group`) and `weight` (numeric, sums to ~100 in
#'   percentage-points or ~1 if already normalised). Internally normalised to
#'   sum 1 before being passed to Stan.
#' @param ref Name of the reference sector for identification (default
#'   "Housing"). The Stan model fixes `lambda[t, ref] = 1` for all t.
#' @param start Earliest month to include in the latent monthly clock.
#' @param monthly_cutover First month where monthly data takes precedence over
#'   quarterly (default 2024-04-01, the ABS Monthly CPI back-series start).
#' @return A list ready to pass to cmdstanr's `$sample(data = ...)`. Fields:
#'   - `T`, `N`, `ref` (integer), `w` (length-N, sums to 1)
#'   - `n_obs_m`, `t_m`, `i_m`, `y_m` (monthly observations)
#'   - `n_obs_q`, `t_q`, `i_q`, `y_q` (quarterly observations)
#'   And, for diagnostics:
#'   - `y` (T x N matrix), `obs_type` (T x N integer matrix),
#'     `dates` (length-T Date vector), `groups` (length-N character vector).
build_stan_data <- function(infl_demeaned,
                            weights,
                            ref = "Housing",
                            start = as.Date("1993-01-01"),
                            monthly_cutover = as.Date("2024-04-01")) {
  groups <- levels(infl_demeaned$group)
  N <- length(groups)
  stopifnot("Expected N = 11 groups or 33 sub-groups" = N %in% c(11L, 33L))

  # Align weights to the group ordering used by infl_demeaned, normalise to 1.
  stopifnot(all(c("group", "weight") %in% names(weights)))
  w_lookup <- setNames(weights$weight, as.character(weights$group))
  stopifnot("All groups must have weights" =
              all(groups %in% names(w_lookup)))
  w <- as.numeric(w_lookup[groups])
  w <- w / sum(w)

  # Reference sector index (1-based, matching the position in `groups`).
  ref_idx <- match(ref, groups)
  stopifnot("ref must be one of the groups" = !is.na(ref_idx))

  # Most recent observed month across all series (monthly takes precedence).
  end <- max(infl_demeaned$date)
  # Align `end` to the first of its month so the seq() below is clean.
  end <- as.Date(format(end, "%Y-%m-01"))

  dates <- seq(start, end, by = "month")
  T_ <- length(dates)

  y <- matrix(0, nrow = T_, ncol = N)
  obs_type <- matrix(0L, nrow = T_, ncol = N)
  colnames(y) <- groups
  colnames(obs_type) <- groups

  # 1) Place monthly observations (obs_type = 1) on or after the cutover.
  m_df <- infl_demeaned |>
    dplyr::filter(.data$source == "monthly", .data$date >= monthly_cutover)
  for (i in seq_len(nrow(m_df))) {
    t_idx <- match(as.Date(format(m_df$date[i], "%Y-%m-01")), dates)
    n_idx <- match(as.character(m_df$group[i]), groups)
    if (!is.na(t_idx) && !is.na(n_idx)) {
      y[t_idx, n_idx] <- m_df$inflation[i]
      obs_type[t_idx, n_idx] <- 1L
    }
  }

  # 2) Place quarterly observations (obs_type = 2) BEFORE the cutover, at the
  #    last month of each quarter. ABS dates a quarter by the first day of its
  #    last month (e.g. Q1 2020 = 2020-03-01) which is already what we want.
  q_df <- infl_demeaned |>
    dplyr::filter(.data$source == "quarterly", .data$date < monthly_cutover)
  for (i in seq_len(nrow(q_df))) {
    t_idx <- match(as.Date(format(q_df$date[i], "%Y-%m-01")), dates)
    n_idx <- match(as.character(q_df$group[i]), groups)
    if (!is.na(t_idx) && !is.na(n_idx) && obs_type[t_idx, n_idx] == 0L) {
      y[t_idx, n_idx] <- q_df$inflation[i]
      obs_type[t_idx, n_idx] <- 2L
    }
  }

  # Flat observation arrays for Stan. Drop quarterly cells with t < 3 (cannot
  # construct the 3-month look-back).
  m_idx <- which(obs_type == 1L, arr.ind = TRUE)
  q_idx <- which(obs_type == 2L & row(obs_type) >= 3L, arr.ind = TRUE)

  list(
    T = T_,
    N = N,
    ref = ref_idx,
    w = w,
    n_obs_m = nrow(m_idx),
    t_m = as.integer(m_idx[, 1]),
    i_m = as.integer(m_idx[, 2]),
    y_m = as.numeric(y[m_idx]),
    n_obs_q = nrow(q_idx),
    t_q = as.integer(q_idx[, 1]),
    i_q = as.integer(q_idx[, 2]),
    y_q = as.numeric(y[q_idx]),
    # Diagnostics:
    y = y,
    obs_type = obs_type,
    dates = dates,
    groups = groups
  )
}

#' Build a QUARTERLY-CLOCKED stan_data list for the sub-group MCT model.
#'
#' Unlike `build_stan_data()` (which uses a monthly latent clock to
#' accommodate mixed-frequency observations at the group level), this
#' builder operates entirely at quarterly resolution. The latent
#' common-factor and sector-trend processes are quarterly, observations
#' are quarterly, and there is no `obs_type` array — every (t, i) cell
#' is either observed (finite y) or missing (NA).
#'
#' This is the right structure for sub-group analysis where the ABS only
#' publishes quarterly data, and avoids the underidentification that
#' arises when you try to estimate 399 monthly latent states from ~125
#' quarterly observations per series.
#'
#' @param infl_demeaned Output of `demean_series()` applied to a
#'   quarterly-only inflation frame (e.g. from
#'   `fetch_cpi_quarterly_subgroup()`).
#' @param weights Tibble with columns `group` and `weight` (sums to ~100).
#' @param ref Name of the reference series for identification
#'   (default "Rents" — the largest sub-group within Housing observed
#'   continuously back to the early 1970s).
#' @param start Earliest quarter to include (default 1993-Q1).
#' @return A list with fields:
#'   - `T`, `N`, `ref` (integer), `w` (length-N, sums to 1)
#'   - `n_obs_q`, `t_q`, `i_q`, `y_q` (flat quarterly obs arrays)
#'   - `y` (T x N matrix with NA for missing), `dates` (length-T Date
#'     vector — first day of each quarter), `groups` (length-N
#'     character vector). For API compatibility with the monthly
#'     `build_stan_data()` we also include zero-length `n_obs_m`, `t_m`,
#'     `i_m`, `y_m` and an `obs_type` matrix with 2 (observed) or
#'     0 (missing).
build_stan_data_quarterly <- function(infl_demeaned,
                                      weights,
                                      ref = "Rents",
                                      start = as.Date("1993-01-01")) {
  groups <- levels(infl_demeaned$group)
  N <- length(groups)
  stopifnot("Expected N = 33 sub-groups" = N == 33L)

  stopifnot(all(c("group", "weight") %in% names(weights)))
  w_lookup <- setNames(weights$weight, as.character(weights$group))
  stopifnot("All groups must have weights" =
              all(groups %in% names(w_lookup)))
  w <- as.numeric(w_lookup[groups]) / sum(as.numeric(w_lookup[groups]))

  ref_idx <- match(ref, groups)
  stopifnot("ref must be one of the groups" = !is.na(ref_idx))

  end <- max(infl_demeaned$date)
  # Align both ends to first-of-quarter month (ABS quarter labels use
  # the LAST month of the quarter — e.g. Q1 = March, Q2 = June). We
  # build a sequence on those last-month-of-quarter dates so the obs
  # join lines up cleanly.
  align_quarter_end <- function(d) {
    m <- as.integer(format(d, "%m"))
    qend_month <- ((m - 1L) %/% 3L) * 3L + 3L
    as.Date(sprintf("%d-%02d-01", as.integer(format(d, "%Y")), qend_month))
  }
  start_q <- align_quarter_end(start)
  end_q   <- align_quarter_end(end)
  dates   <- seq(start_q, end_q, by = "3 months")
  T_      <- length(dates)

  y <- matrix(NA_real_, nrow = T_, ncol = N)
  colnames(y) <- groups

  q_df <- infl_demeaned |>
    dplyr::filter(.data$source == "quarterly")
  for (i in seq_len(nrow(q_df))) {
    t_idx <- match(as.Date(format(q_df$date[i], "%Y-%m-01")), dates)
    n_idx <- match(as.character(q_df$group[i]), groups)
    if (!is.na(t_idx) && !is.na(n_idx)) {
      y[t_idx, n_idx] <- q_df$inflation[i]
    }
  }

  q_idx <- which(!is.na(y), arr.ind = TRUE)

  # obs_type: 2 = quarterly observed, 0 = missing — kept for downstream
  # code that expects this field (postprocess.R reads it for diagnostics).
  obs_type <- matrix(0L, nrow = T_, ncol = N)
  obs_type[!is.na(y)] <- 2L

  list(
    T   = T_,
    N   = N,
    ref = ref_idx,
    w   = w,
    # For a quarterly-clocked panel each "quarterly" observation IS the
    # single-period latent, not a 3-month average. Surface obs in the
    # `*_m` slots so the downstream `.compute_derived_gq` log_lik logic
    # (which expects direct observations there) works without further
    # case-splitting. `n_obs_q = 0` ensures the avg-of-3 quarterly path
    # is skipped.
    n_obs_m = nrow(q_idx),
    t_m = as.integer(q_idx[, 1]),
    i_m = as.integer(q_idx[, 2]),
    y_m = as.numeric(y[q_idx]),
    n_obs_q = 0L,
    t_q = integer(0),
    i_q = integer(0),
    y_q = numeric(0),
    y        = y,
    obs_type = obs_type,
    dates    = dates,
    groups   = groups,
    clock    = "quarterly"
  )
}
