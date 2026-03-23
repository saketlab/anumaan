# plots.R
# GBD-style burden plots for AMR analysis
# Replicates Figures 4 and 6 from Antimicrobial Resistance Collaborators,
# Lancet 2022, and provides hospital-level DALY comparison plots.

# Shared theme for burden plots
.burden_theme <- function(base_size = 16, grid_x = TRUE, legend_pos = "right") {
  th <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 2
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 2,
        color = "grey40"
      ),
      legend.position = legend_pos,
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = base_size - 1
      )
    )
  if (!grid_x) {
    th <- th + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  }
  th
}


#' Compute Hospital-Level DALY Breakdown
#'
#' Distributes pooled YLL/YLD totals to individual hospitals proportionally
#' by deaths and discharges, and computes per-1000 rates.
#'
#' @param hospital_counts Data frame with columns: \code{center_name},
#'   \code{deaths_h}, \code{discharged_h}, \code{cases_h}.
#' @param total_deaths Numeric. Pooled total deaths across all hospitals.
#' @param total_discharged Numeric. Pooled total discharges.
#' @param yll_base Numeric. Total baseline YLL.
#' @param yll_associated Numeric. Total YLL associated with resistance.
#' @param yll_attributable Numeric. Total YLL attributable to resistance.
#' @param yld_base Numeric. Total baseline YLD.
#' @param yld_associated Numeric. Total YLD associated.
#' @param yld_attributable Numeric. Total YLD attributable.
#'
#' @return Data frame with hospital-level YLL, YLD, DALY (absolute and
#'   per 1000 cases) columns.
#' @export
compute_hospital_daly <- function(hospital_counts,
                                  total_deaths,
                                  total_discharged,
                                  yll_base,
                                  yll_associated,
                                  yll_attributable,
                                  yld_base,
                                  yld_associated,
                                  yld_attributable) {
  hospital_counts %>%
    dplyr::mutate(
      .death_frac = deaths_h / total_deaths,
      .disch_frac = discharged_h / total_discharged,
      .per_1k = 1000 / cases_h,
      # Absolute
      YLL_base = yll_base * .data$.death_frac,
      YLD_base = yld_base * .data$.disch_frac,
      YLL_associated = yll_associated * .data$.death_frac,
      YLD_associated = yld_associated * .data$.disch_frac,
      DALY_associated = .data$YLL_associated + .data$YLD_associated,
      YLL_attributable = yll_attributable * .data$.death_frac,
      YLD_attributable = yld_attributable * .data$.disch_frac,
      DALY_attributable = .data$YLL_attributable + .data$YLD_attributable,
      # Per 1000 cases
      YLL_base_per_1000 = .data$YLL_base * .data$.per_1k,
      YLD_base_per_1000 = .data$YLD_base * .data$.per_1k,
      YLL_assoc_per_1000 = .data$YLL_associated * .data$.per_1k,
      YLD_assoc_per_1000 = .data$YLD_associated * .data$.per_1k,
      DALY_assoc_per_1000 = .data$DALY_associated * .data$.per_1k,
      YLL_attr_per_1000 = .data$YLL_attributable * .data$.per_1k,
      YLD_attr_per_1000 = .data$YLD_attributable * .data$.per_1k,
      DALY_attr_per_1000 = .data$DALY_attributable * .data$.per_1k
    ) %>%
    dplyr::select(-".death_frac", -".disch_frac", -".per_1k") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4)))
}


#' GBD Figure 4 — Deaths by Organism Group (Overlapping Bars)
#'
#' Overlapping bar chart showing deaths associated with vs attributable to AMR,
#' grouped by organism, replicating GBD Lancet 2022 Figure 4.
#'
#' @param deaths_df Data frame with columns: \code{org_group},
#'   \code{deaths_associated}, \code{deaths_attributable}.
#' @param total_deaths Numeric. Total deaths for subtitle.
#' @param title Character. Plot title. Default generates a standard title.
#'
#' @return A \code{ggplot} object.
#' @export
plot_gbd_fig4 <- function(deaths_df, total_deaths, title = NULL) {
  if (is.null(title)) {
    title <- "Deaths Attributable to and Associated with AMR by Pathogen Group"
  }

  plot_df <- deaths_df %>%
    dplyr::mutate(
      org_group = factor(.data$org_group,
        levels = .data$org_group[order(.data$deaths_associated,
          decreasing = TRUE
        )]
      ),
      label_associated = round_to_sum(
        .data$deaths_associated,
        round(sum(.data$deaths_associated))
      ),
      label_attributable = round_to_sum(
        .data$deaths_attributable,
        round(sum(.data$deaths_attributable))
      )
    )

  pal <- get_amr_palette("burden_gbd4")

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = plot_df,
      ggplot2::aes(
        x = .data$org_group, y = .data$deaths_associated,
        fill = "Associated with resistance"
      ),
      width = 0.7
    ) +
    ggplot2::geom_col(
      data = plot_df,
      ggplot2::aes(
        x = .data$org_group, y = .data$deaths_attributable,
        fill = "Attributable to resistance"
      ),
      width = 0.4
    ) +
    ggplot2::geom_text(
      data = plot_df,
      ggplot2::aes(
        x = .data$org_group, y = .data$deaths_associated,
        label = .data$label_associated
      ),
      vjust = -0.5, size = 5, fontface = "bold",
      color = unname(pal["Associated with resistance"])
    ) +
    ggplot2::geom_text(
      data = plot_df,
      ggplot2::aes(
        x = .data$org_group, y = .data$deaths_attributable,
        label = .data$label_attributable
      ),
      vjust = -0.5, size = 4.5, fontface = "bold",
      color = unname(pal["Attributable to resistance"])
    ) +
    ggplot2::scale_fill_manual(values = pal, name = "Resistance") +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = sprintf("Total deaths = %d", total_deaths),
      x        = "Pathogen group",
      y        = "Deaths (count)"
    ) +
    .burden_theme(legend_pos = c(0.98, 0.98)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 40, hjust = 1,
        vjust = 1
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.justification = c("right", "top"),
      legend.background = ggplot2::element_rect(
        color = "grey70",
        linewidth = 0.4
      ),
      legend.key.size = ggplot2::unit(0.6, "cm")
    )
}


#' GBD Figure 6 — Attributable Deaths Heatmap (Pathogen x Drug Class)
#'
#' Heatmap of deaths attributable to AMR by pathogen and antibiotic class,
#' replicating GBD Lancet 2022 Figure 6. Includes an "All pathogens" summary
#' row.
#'
#' @param profile_df Data frame with columns: \code{pathogen}, \code{profile}
#'   (drug class), \code{deaths_attributable}.
#' @param total_deaths Numeric. Total deaths for subtitle.
#' @param total_attributable Numeric. Total attributable deaths for subtitle.
#' @param title Character. Plot title.
#'
#' @return A \code{ggplot} object.
#' @export
plot_gbd_fig6 <- function(profile_df, total_deaths, total_attributable,
                          title = NULL) {
  if (is.null(title)) {
    title <- "Deaths Attributable to AMR by Pathogen \u2013 Drug Class"
  }

  fig6_raw <- profile_df %>%
    dplyr::mutate(drug_class = shorten_drug_class(.data$profile))

  all_row <- fig6_raw %>%
    dplyr::group_by(.data$drug_class) %>%
    dplyr::summarise(
      deaths_attributable = sum(.data$deaths_attributable,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(pathogen = "All pathogens")

  fig6_df <- dplyr::bind_rows(fig6_raw, all_row)

  # Ordering
  pathogen_order <- fig6_raw %>%
    dplyr::group_by(.data$pathogen) %>%
    dplyr::summarise(
      total = sum(.data$deaths_attributable, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$total)) %>%
    dplyr::pull(.data$pathogen)

  pathogen_levels <- c(rev(pathogen_order), "All pathogens")

  class_order <- all_row %>%
    dplyr::arrange(dplyr::desc(.data$deaths_attributable)) %>%
    dplyr::pull(.data$drug_class)

  fig6_df <- fig6_df %>%
    dplyr::mutate(
      pathogen = factor(.data$pathogen, levels = pathogen_levels),
      drug_class = factor(.data$drug_class, levels = class_order),
      .rounded = round(.data$deaths_attributable),
      cell_label = ifelse(.data$.rounded == 0, "",
        as.character(.data$.rounded)
      ),
      deaths_plot = ifelse(.data$.rounded == 0, NA_real_,
        .data$deaths_attributable
      )
    )

  fig6_with_value <- dplyr::filter(fig6_df, !is.na(.data$deaths_plot))
  max_deaths <- max(fig6_df$deaths_plot, na.rm = TRUE)

  ggplot2::ggplot(
    fig6_df,
    ggplot2::aes(x = .data$drug_class, y = .data$pathogen)
  ) +
    ggplot2::geom_tile(
      data = fig6_with_value,
      ggplot2::aes(fill = .data$deaths_plot),
      color = "grey80", linewidth = 0.35
    ) +
    ggplot2::geom_tile(
      data = dplyr::filter(
        fig6_with_value,
        .data$pathogen == "All pathogens"
      ),
      ggplot2::aes(fill = .data$deaths_plot),
      color = "grey30", linewidth = 0.9
    ) +
    ggplot2::geom_text(
      data = fig6_with_value,
      ggplot2::aes(
        label = .data$cell_label,
        color = .data$deaths_plot > max_deaths * 0.55
      ),
      size = 4.5, fontface = "bold"
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "white", "FALSE" = "black"), guide = "none"
    ) +
    ggplot2::scale_fill_gradientn(
      colours  = get_amr_palette("burden_gbd6"),
      na.value = "transparent",
      name     = "Attributable\nDeaths (count)",
      limits   = c(0, NA)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "Total deaths = %d  |  Attributable = %d",
        total_deaths, round(total_attributable)
      ),
      x = "Antibiotic class", y = NULL
    ) +
    .burden_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.border = ggplot2::element_rect(
        color = "grey50", fill = NA,
        linewidth = 0.6
      )
    )
}


#' Burden Comparison Bar Plot
#'
#' Horizontal grouped bar chart comparing associated vs attributable burden
#' (YLL, YLD, or DALY) across hospitals.
#'
#' @param hospital_daly Data frame from \code{compute_hospital_daly()}.
#' @param metric Character. One of \code{"YLL"}, \code{"YLD"}, or
#'   \code{"DALY"}.
#' @param label_map Optional named character vector mapping center names to
#'   display labels. Default \code{NULL} (use center names as-is).
#' @param syndrome_label Character. Label for the infectious syndrome used
#'   in axis/legend text (e.g., \code{"BSI patients"}, \code{"UTI cases"}).
#'   Default \code{"cases"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_burden_comparison <- function(hospital_daly, metric = "DALY",
                                   label_map = NULL,
                                   syndrome_label = "cases") {
  assoc_col <- paste0(metric, "_assoc_per_1000")
  attr_col <- paste0(metric, "_attr_per_1000")

  assoc_label <- paste(metric, "Associated")
  attr_label <- paste(metric, "Attributable")

  plot_df <- hospital_daly %>%
    dplyr::select(
      "center_name",
      !!assoc_label := dplyr::all_of(assoc_col),
      !!attr_label := dplyr::all_of(attr_col)
    ) %>%
    tidyr::pivot_longer(-"center_name",
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      center_name = reorder(.data$center_name, .data$value, max),
      metric = factor(.data$metric,
        levels = c(assoc_label, attr_label)
      )
    )

  pal_name <- paste0("burden_", tolower(metric))
  colors <- get_amr_palette(pal_name)
  names(colors) <- c(assoc_label, attr_label)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$center_name, y = .data$value, fill = .data$metric)
  ) +
    ggplot2::geom_col(
      position = "dodge", color = "grey20", linewidth = 0.3,
      width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data$value, 1)),
      position = ggplot2::position_dodge(width = 0.7),
      hjust = -0.15, size = 5, fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colors, name = "") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::labs(
      title = paste(metric, "Associated vs Attributable to AMR"),
      subtitle = paste("Per 1 000", syndrome_label),
      x = "", y = paste(metric, "per 1 000", syndrome_label)
    ) +
    .burden_theme(grid_x = FALSE, legend_pos = "bottom")

  if (!is.null(label_map)) {
    p <- p + ggplot2::scale_x_discrete(labels = label_map)
  }

  p
}


#' Burden Heatmap Across Hospitals
#'
#' Heatmap showing all burden metrics (YLL, YLD, DALY x associated,
#' attributable) per 1000 cases across hospitals.
#'
#' @param hospital_daly Data frame from \code{compute_hospital_daly()}.
#' @param label_map Optional named character vector for center name labels.
#' @param syndrome_label Character. Label for infectious syndrome in legends.
#'   Default \code{"cases"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_burden_heatmap <- function(hospital_daly, label_map = NULL,
                                syndrome_label = "cases") {
  metric_levels <- c(
    "YLL\nAssociated", "YLD\nAssociated", "DALY\nAssociated",
    "YLL\nAttributable", "YLD\nAttributable", "DALY\nAttributable"
  )

  heatmap_df <- hospital_daly %>%
    dplyr::select(
      "center_name",
      `YLL\nAssociated`    = "YLL_assoc_per_1000",
      `YLD\nAssociated`    = "YLD_assoc_per_1000",
      `DALY\nAssociated`   = "DALY_assoc_per_1000",
      `YLL\nAttributable`  = "YLL_attr_per_1000",
      `YLD\nAttributable`  = "YLD_attr_per_1000",
      `DALY\nAttributable` = "DALY_attr_per_1000"
    ) %>%
    tidyr::pivot_longer(-"center_name",
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      center_name = reorder(.data$center_name, .data$value, max),
      metric = factor(.data$metric, levels = metric_levels)
    )

  fill_mid <- max(heatmap_df$value, na.rm = TRUE) / 2

  p <- ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(x = .data$metric, y = .data$center_name, fill = .data$value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = round(.data$value, 1),
        color = .data$value > fill_mid
      ),
      size = 5, fontface = "bold"
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "white", "FALSE" = "black"), guide = "none"
    ) +
    ggplot2::scale_fill_distiller(
      palette = "Spectral", direction = -1, na.value = "grey92",
      name = paste0("Per 1 000\n", syndrome_label)
    ) +
    ggplot2::labs(
      title = "Burden Metrics by Hospital",
      subtitle = paste("All values normalised per 1 000", syndrome_label),
      x = "", y = ""
    ) +
    .burden_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "grey60", fill = NA,
        linewidth = 0.6
      )
    )

  if (!is.null(label_map)) {
    p <- p + ggplot2::scale_y_discrete(labels = label_map)
  }

  p
}


#' DALY Lollipop Plot
#'
#' Lollipop chart ranking hospitals by DALY (associated or attributable)
#' per 1000 cases, with point size proportional to case count.
#'
#' @param hospital_daly Data frame from \code{compute_hospital_daly()}.
#' @param type Character. \code{"attributable"} or \code{"associated"}.
#'   Default \code{"attributable"}.
#' @param label_map Optional named character vector for center name labels.
#' @param syndrome_label Character. Label for infectious syndrome in
#'   axis/legend text. Default \code{"cases"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_daly_lollipop <- function(hospital_daly, type = "attributable",
                               label_map = NULL,
                               syndrome_label = "cases") {
  y_col <- if (type == "attributable") {
    "DALY_attr_per_1000"
  } else {
    "DALY_assoc_per_1000"
  }

  pal_name <- paste0(
    "burden_lollipop_",
    if (type == "attributable") "attr" else "assoc"
  )
  colors <- get_amr_palette(pal_name)

  plot_df <- hospital_daly %>%
    dplyr::mutate(
      center_name = reorder(.data$center_name, .data[[y_col]]),
      y_val = .data[[y_col]]
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$center_name,
      y = .data$y_val
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$center_name, y = 0, yend = .data$y_val),
      color = "grey50", linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$cases_h, color = .data$y_val)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data$y_val, 1)),
      hjust = -0.5, size = 5, fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_gradient(
      low = colors[1], high = colors[2],
      name = paste0("DALY\n", type, "\nper 1 000")
    ) +
    ggplot2::scale_size_continuous(
      name = paste0("Total\n", syndrome_label),
      range = c(4, 10)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::labs(
      title = paste("DALY", tools::toTitleCase(type), "to AMR per Hospital"),
      subtitle = paste(
        "Per 1 000", syndrome_label,
        " |  Point size = total", syndrome_label
      ),
      x = "", y = paste("DALY", type, "per 1 000", syndrome_label)
    ) +
    .burden_theme(grid_x = FALSE)

  if (!is.null(label_map)) {
    p <- p + ggplot2::scale_x_discrete(labels = label_map)
  }

  p
}


#' Resistance Fraction Stacked Bar
#'
#' 100% stacked bar chart showing what percentage of DALY burden is
#' attributable to resistance vs baseline infectious burden, per hospital.
#'
#' @param hospital_daly Data frame from \code{compute_hospital_daly()}.
#' @param label_map Optional named character vector for center name labels.
#'
#' @return A \code{ggplot} object.
#' @export
plot_resistance_fraction <- function(hospital_daly, label_map = NULL) {
  plot_df <- hospital_daly %>%
    dplyr::mutate(
      pct_attributable = ifelse(.data$DALY_associated == 0, 0,
        .data$DALY_attributable / .data$DALY_associated * 100
      ),
      pct_baseline = 100 - .data$pct_attributable,
      center_name = reorder(.data$center_name, .data$pct_attributable)
    ) %>%
    dplyr::select("center_name", "pct_attributable", "pct_baseline") %>%
    tidyr::pivot_longer(-"center_name",
      names_to = "fraction",
      values_to = "pct"
    ) %>%
    dplyr::mutate(
      fraction = dplyr::recode(.data$fraction,
        pct_attributable = "Attributable to resistance",
        pct_baseline = "Baseline infectious burden"
      )
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$center_name, y = .data$pct, fill = .data$fraction)
  ) +
    ggplot2::geom_col(color = "white", linewidth = 0.3, width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(.data$pct, 1), "%")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 5, fontface = "bold", color = "white"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = get_amr_palette("burden_fraction"),
      name = ""
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(scale = 1)
    ) +
    ggplot2::labs(
      title = "Resistance Fraction of Total DALY Burden",
      subtitle = "What % of the infectious disease DALY burden is attributable to AMR?",
      x = "", y = "% of total DALY (associated)"
    ) +
    .burden_theme(grid_x = FALSE, legend_pos = "bottom")

  if (!is.null(label_map)) {
    p <- p + ggplot2::scale_x_discrete(labels = label_map)
  }

  p
}
