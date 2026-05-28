# ggplot helpers for the dashboard. Visual style is deliberately spare
# (single axis colour, light grid, no chart-junk) — inspired by Liberty
# Street Economics.

#' Theme used across the dashboard.
mct_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(linewidth = 0.2, colour = "grey90"),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.2, colour = "grey90"),
      plot.title         = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
      plot.subtitle      = ggplot2::element_text(colour = "grey40", size = ggplot2::rel(0.95)),
      plot.caption       = ggplot2::element_text(colour = "grey50", hjust = 0, size = ggplot2::rel(0.85)),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_blank()
    )
}

#' Trend line with credible band, optionally overlaid with reference
#' series (e.g. headline CPI YoY). `summary` is `summarise_path()`'s
#' output (date, lower, median, upper). `overlays` is an optional named
#' list of tibbles each with columns (date, value) — drawn as
#' additional lines with a legend.
plot_trend_band <- function(summary, title = "Trend inflation",
                            subtitle = NULL, caption = NULL,
                            colour = "#1f5582",
                            ylab = "Annualised % (demeaned)",
                            overlays = NULL, since = NULL,
                            trend_label = "Multivariate Core Trend",
                            band_label = "MCT 68% credible band",
                            zero_line = TRUE) {
  if (!is.null(since)) {
    summary <- summary[summary$date >= as.Date(since), , drop = FALSE]
    if (!is.null(overlays)) {
      overlays <- lapply(overlays, function(d)
        d[d$date >= as.Date(since), , drop = FALSE])
    }
  }
  p <- ggplot2::ggplot()
  if (zero_line) {
    p <- p + ggplot2::geom_hline(yintercept = 0,
                                 colour = "grey70", linewidth = 0.3)
  }
  p <- p +
    ggplot2::geom_ribbon(
      data = summary,
      ggplot2::aes(.data$date,
                   ymin = .data$lower, ymax = .data$upper,
                   fill = band_label),
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = summary,
      ggplot2::aes(.data$date, .data$median,
                   colour = trend_label, linetype = trend_label),
      linewidth = 0.7
    )

  # Overlay palette: each line gets a distinctly different hue so it's
  # easy to tell apart at a glance (vs. NY Fed's two-greys-by-linetype
  # convention, which is hard to read on screens).
  overlay_styles <- list(
    "Headline CPI (YoY)"     = list(colour = "#b45b1f",  # burnt orange
                                    linetype = "dashed"),
    "Trimmed mean CPI (YoY)" = list(colour = "#2e7d32",  # forest green
                                    linetype = "solid")
  )
  if (!is.null(overlays)) {
    # Build a SINGLE long tibble with a `series` column so ggplot can
    # map colour + linetype + legend correctly. Building geom_line in
    # a for-loop hits a lazy-evaluation bug — the `nm` variable inside
    # aes() captures the LAST iteration's value, collapsing all series
    # to the same colour.
    overlay_long <- dplyr::bind_rows(
      lapply(names(overlays), function(nm) {
        d <- overlays[[nm]]
        stopifnot(all(c("date", "value") %in% names(d)))
        d$series <- nm
        d[, c("date", "value", "series")]
      })
    )
    overlay_long$series <- factor(overlay_long$series,
                                  levels = names(overlays))
    p <- p + ggplot2::geom_line(
      data = overlay_long,
      ggplot2::aes(.data$date, .data$value,
                   colour = .data$series, linetype = .data$series),
      linewidth = 0.7
    )
  }

  # Build colour + linetype + fill scales with explicit ordering so the
  # legend always shows MCT first, then overlays in input order.
  ordered_lines <- c(trend_label, names(overlays))
  line_colours <- c(setNames(colour, trend_label),
                    vapply(names(overlays),
                           function(nm) overlay_styles[[nm]]$colour %||% "#444",
                           character(1)))
  line_ltys    <- c(setNames("solid", trend_label),
                    vapply(names(overlays),
                           function(nm) overlay_styles[[nm]]$linetype %||% "dotted",
                           character(1)))

  p +
    ggplot2::scale_colour_manual(
      name = NULL, breaks = ordered_lines, values = line_colours,
      guide = ggplot2::guide_legend(order = 1)
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL, breaks = ordered_lines, values = line_ltys,
      guide = ggplot2::guide_legend(order = 1)
    ) +
    ggplot2::scale_fill_manual(
      name = NULL, values = setNames(colour, band_label),
      guide = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = ylab) +
    mct_theme() +
    ggplot2::theme(legend.position = "top",
                   legend.key.width = ggplot2::unit(1.4, "lines"))
}

# Local null-default operator for clean `value %||% default` syntax.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' NY-Fed-style sector-decomposition chart: dashed MCT trend (median)
#' overlaid with per-bucket Common + Sector-specific lines (3 buckets ×
#' 2 components = 6 lines + 1 trend = 7 total).
#'
#' @param decomp Long tibble from `decompose_trend_by_bucket()` (with
#'   columns date, bucket, component, median, ...).
#' @param trend  Long tibble (date, median) of the headline trend
#'   (demeaned or raw — pick to match `decomp` for visual consistency).
#' @param since Optional Date — restrict to dates >= since.
plot_sector_decomposition <- function(decomp, trend,
                                      title = "Sector decomposition of trend inflation",
                                      subtitle = NULL, caption = NULL,
                                      since = NULL) {
  if (!is.null(since)) {
    decomp <- decomp[decomp$date >= as.Date(since), , drop = FALSE]
    trend  <- trend [trend$date  >= as.Date(since), , drop = FALSE]
  }

  # Per-bucket palette (matches NY Fed's blue / red / brown).
  bucket_cols <- c(
    "Goods"                = "#1f77b4",   # blue
    "Services ex. housing" = "#a63a3a",   # muted red
    "Housing"              = "#8a6b1f"    # warm brown
  )

  # Build a series-level palette by combining bucket × component:
  # Common = saturated bucket colour, Sector-specific = paler shade.
  to_pale <- function(hex, alpha = 0.45) {
    rgb <- col2rgb(hex)[, 1L] / 255
    white <- c(1, 1, 1)
    pale <- rgb + (white - rgb) * (1 - alpha)
    rgb(pale[1L], pale[2L], pale[3L])
  }
  series_names <- character(0)
  series_vals  <- character(0)
  for (b in names(bucket_cols)) {
    series_names <- c(series_names,
                      paste(b, "Common", sep = " · "),
                      paste(b, "Sector-specific", sep = " · "))
    series_vals  <- c(series_vals,
                      bucket_cols[[b]],
                      to_pale(bucket_cols[[b]]))
  }
  series_names <- c(series_names, "Multivariate Core Trend")
  series_vals  <- c(series_vals,  "#1d1d1f")
  series_cols  <- setNames(series_vals, series_names)

  # Decomp series labels (6 lines: 3 buckets × 2 components).
  # Keep `series` as character (NOT factor) so the trend layer's
  # series can join the same colour scale without level mismatches.
  decomp_only <- decomp |>
    dplyr::mutate(series = paste(.data$bucket, .data$component,
                                 sep = " · "))
  decomp_series_levels <- names(series_cols)[
    names(series_cols) != "Multivariate Core Trend"]

  # Trend layer (one line, dashed, heavier)
  trend_line <- tibble::tibble(date   = trend$date,
                               median = trend$median,
                               series = "Multivariate Core Trend")

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70",
                        linewidth = 0.3) +
    # Decomposition lines (6)
    ggplot2::geom_line(
      data = decomp_only,
      ggplot2::aes(.data$date, .data$median, colour = .data$series),
      linewidth = 0.6
    ) +
    # MCT trend overlay (dashed, heavier)
    ggplot2::geom_line(
      data = trend_line,
      ggplot2::aes(.data$date, .data$median, colour = .data$series),
      linetype = "dashed", linewidth = 0.85
    ) +
    ggplot2::scale_colour_manual(
      name = NULL, values = series_cols,
      breaks = c("Multivariate Core Trend", decomp_series_levels)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = "Percentage points (annualised)") +
    mct_theme() +
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
                   legend.key.width = ggplot2::unit(1.4, "lines")) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 4))
}

#' Common-share path — same shape as trend, scaled to [0, 1].
#' @param since If non-NULL, restrict the plotted series to dates >= `since`.
#'   Useful for the common-share metric on this dataset — pre-monthly-data
#'   (i.e. before Apr 2024) the quarterly-only obs structure can't reliably
#'   separate common-factor from sector-specific volatility, so the SV
#'   processes drift and the chart frequently pins at the 0% / 100% bounds.
plot_common_share <- function(summary, title = "Common share of inflation variance",
                              subtitle = NULL, caption = NULL,
                              colour = "#7a4e7d", since = NULL,
                              ylim = NULL) {
  if (!is.null(since)) {
    summary <- summary[summary$date >= as.Date(since), , drop = FALSE]
  }
  # Default y-limits: a tight window around the actual data range so
  # small movements are visible (the [0%, 100%] full range hides the
  # action when common-share sits in a narrow band like 80-95%).
  if (is.null(ylim)) {
    pad <- 0.02
    lo  <- max(0, floor((min(summary$lower) - pad) * 20) / 20)
    hi  <- min(1, ceiling((max(summary$upper) + pad) * 20) / 20)
    ylim <- c(lo, hi)
  }
  ggplot2::ggplot(summary, ggplot2::aes(.data$date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = colour, alpha = 0.20) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median), colour = colour,
                       linewidth = 0.7) +
    ggplot2::scale_y_continuous(limits = ylim,
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = NULL) +
    mct_theme()
}

#' Compare two-or-more model paths (e.g., trend_A vs trend_B vs trend_C) on
#' the same axes. Accepts either:
#'   - the legacy pairwise signature `(summary_A, summary_B, ...)`, OR
#'   - a named list of summaries as the first argument
#' so the existing two-variant call sites keep working unchanged.
plot_compare_paths <- function(summaries, ..., labels = NULL,
                               title = "Variant comparison",
                               subtitle = NULL, caption = NULL,
                               colours = c("#1f5582", "#b45b1f", "#2e7d32",
                                           "#7a4e7d", "#c0392b"),
                               ylab = "Annualised % (demeaned)") {
  # Back-compat: positional args `(summary_A, summary_B[, summary_C, ...])`
  # land here as `summaries` (= summary_A) and `...` (the rest).
  extras <- list(...)
  if (!is.list(summaries) || all(c("date", "median", "lower", "upper") %in%
                                 names(summaries))) {
    # `summaries` is itself a single summary tibble; collect the others.
    summary_list <- c(list(summaries), extras)
    if (is.null(labels)) labels <- LETTERS[seq_along(summary_list)]
  } else {
    summary_list <- summaries
    if (is.null(labels)) labels <- names(summary_list)
    if (is.null(labels)) labels <- LETTERS[seq_along(summary_list)]
  }
  stopifnot("Need at least two summaries to compare" =
              length(summary_list) >= 2)
  stopifnot("Need a label per summary" =
              length(labels) == length(summary_list))

  df <- dplyr::bind_rows(lapply(seq_along(summary_list), function(i) {
    dplyr::mutate(summary_list[[i]], variant = labels[i])
  }))
  df <- dplyr::mutate(df, variant = factor(.data$variant, levels = labels))

  ggplot2::ggplot(df, ggplot2::aes(.data$date, colour = .data$variant,
                                   fill = .data$variant)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         alpha = 0.13, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median), linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = colours[seq_along(summary_list)]) +
    ggplot2::scale_fill_manual(values = colours[seq_along(summary_list)]) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = ylab) +
    mct_theme()
}

#' Faceted small-multiples chart of group contributions to trend
#' inflation — one panel per group. `contribs` is the long tibble from
#' `sector_contributions()` (or its parent-aggregated version). Stacked
#' areas don't read cleanly when contributions cross zero or when there
#' are >5-6 series; small multiples scale to 11+ groups without losing
#' detail.
plot_sector_contribs <- function(contribs,
                                 since = NULL,
                                 title = "Group contributions to trend inflation",
                                 subtitle = NULL,
                                 caption = NULL,
                                 colour = "#1f5582") {
  if (!is.null(since)) {
    contribs <- dplyr::filter(contribs, .data$date >= since)
  }
  ggplot2::ggplot(contribs,
                  ggplot2::aes(.data$date, .data$median)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_line(colour = colour, linewidth = 0.6) +
    ggplot2::geom_area(fill = colour, alpha = 0.18) +
    ggplot2::facet_wrap(~ .data$group, ncol = 3, scales = "fixed",
                        labeller = ggplot2::label_wrap_gen(width = 28)) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = "Percentage points (annualised)") +
    mct_theme() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.85),
                                         face = "bold", hjust = 0),
      panel.spacing = ggplot2::unit(0.6, "lines")
    )
}

#' Horizontal bar chart of group contributions at a single date.
#' Sorted descending by median contribution. Use as a "what drove the
#' headline this period?" summary alongside the time-series facets.
plot_contribs_bar <- function(contribs, at_date,
                              title = "Group contributions, latest period",
                              subtitle = NULL,
                              caption = NULL,
                              bar_colour = "#1f5582") {
  d <- dplyr::filter(contribs, .data$date == at_date) |>
    dplyr::arrange(.data$median) |>
    dplyr::mutate(group = factor(.data$group, levels = .data$group))
  total <- sum(d$median)
  sub_total <- sprintf("Sums to %+.2f pp (annualised, demeaned)", total)
  full_subtitle <- if (is.null(subtitle)) sub_total else
    paste(subtitle, sub_total, sep = " · ")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$median, y = .data$group)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_col(fill = bar_colour, alpha = 0.85, width = 0.7) +
    ggplot2::labs(title = title, subtitle = full_subtitle, caption = caption,
                  x = "Contribution to trend (pp, annualised)",
                  y = NULL) +
    mct_theme()
}

#' Horizontal-bar attribution of the trend revision at a single date.
#' `decomp` is the output of `decompose_revisions()`. Bars are sorted by
#' absolute contribution to the revision; positive bars (right of zero) pushed
#' the trend up vs the prior vintage, negative bars pulled it down.
plot_revision_attribution <- function(decomp, at_date,
                                      title = "Revision attribution",
                                      subtitle = NULL,
                                      caption = NULL,
                                      bar_colour = "#1f5582") {
  d <- dplyr::filter(decomp, .data$date == at_date)
  if (nrow(d) == 0) {
    stop(sprintf(
      "decompose_revisions has no rows at date %s; check overlap.",
      format(at_date)
    ))
  }
  d <- d |>
    dplyr::arrange(.data$delta_contrib) |>
    dplyr::mutate(group = factor(.data$group, levels = .data$group))
  total <- sum(d$delta_contrib)
  sub_total <- sprintf("Total revision at this date: %+.3f pp", total)
  full_subtitle <- if (is.null(subtitle)) sub_total else
    paste(subtitle, sub_total, sep = " · ")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$delta_contrib, y = .data$group)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_col(fill = bar_colour, alpha = 0.85, width = 0.7) +
    ggplot2::labs(title = title, subtitle = full_subtitle, caption = caption,
                  x = "Δ contribution to trend (pp, current − prior)",
                  y = NULL) +
    mct_theme()
}

#' Revisions panel: this-month vs last-month trend estimates at the dates
#' both fits cover. `df` is the output of `posterior_revisions()`.
plot_revisions <- function(df, title = "Trend revisions vs prior vintage",
                           subtitle = NULL, caption = NULL,
                           current_colour = "#1f5582",
                           prior_colour = "#a3a3a3") {
  ggplot2::ggplot(df, ggplot2::aes(.data$date)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$curr_lower, ymax = .data$curr_upper),
                         fill = current_colour, alpha = 0.18) +
    ggplot2::geom_line(ggplot2::aes(y = .data$prior_median),
                       colour = prior_colour, linewidth = 0.7,
                       linetype = "22") +
    ggplot2::geom_line(ggplot2::aes(y = .data$curr_median),
                       colour = current_colour, linewidth = 0.7) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = "Annualised % (demeaned)") +
    mct_theme()
}
