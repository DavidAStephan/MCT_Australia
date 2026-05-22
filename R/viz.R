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

#' Trend line with credible band. `summary` must be the output of
#' `summarise_path()` (columns date, lower, median, upper).
plot_trend_band <- function(summary, title = "Trend inflation",
                            subtitle = NULL, caption = NULL,
                            colour = "#1f5582", ylab = "Annualised % (demeaned)") {
  ggplot2::ggplot(summary, ggplot2::aes(.data$date)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = colour, alpha = 0.20) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median), colour = colour,
                       linewidth = 0.7) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = ylab) +
    mct_theme()
}

#' Common-share path — same shape as trend, scaled to [0, 1].
plot_common_share <- function(summary, title = "Common share of inflation variance",
                              subtitle = NULL, caption = NULL,
                              colour = "#7a4e7d") {
  ggplot2::ggplot(summary, ggplot2::aes(.data$date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = colour, alpha = 0.20) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median), colour = colour,
                       linewidth = 0.7) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = NULL) +
    mct_theme()
}

#' Compare two model paths (e.g., trend_A vs trend_B) on the same axes.
plot_compare_paths <- function(summary_A, summary_B, labels = c("A", "B"),
                               title = "Variant comparison",
                               subtitle = NULL, caption = NULL,
                               colours = c("#1f5582", "#b45b1f"),
                               ylab = "Annualised % (demeaned)") {
  df <- dplyr::bind_rows(
    dplyr::mutate(summary_A, variant = labels[1]),
    dplyr::mutate(summary_B, variant = labels[2])
  )
  ggplot2::ggplot(df, ggplot2::aes(.data$date, colour = .data$variant,
                                   fill = .data$variant)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median), linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = ylab) +
    mct_theme()
}

#' Stacked-area chart of sector contributions to trend inflation. `contribs`
#' is the long tibble from `sector_contributions()`.
plot_sector_contribs <- function(contribs,
                                 since = NULL,
                                 title = "Sector contributions to trend inflation",
                                 subtitle = NULL,
                                 caption = NULL) {
  if (!is.null(since)) {
    contribs <- dplyr::filter(contribs, .data$date >= since)
  }
  ggplot2::ggplot(contribs, ggplot2::aes(.data$date, .data$median,
                                         fill = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_area(position = "stack", alpha = 0.85,
                       colour = "white", linewidth = 0.1) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = "Percentage points (annualised)") +
    mct_theme() +
    ggplot2::theme(legend.position = "right",
                   legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)))
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
