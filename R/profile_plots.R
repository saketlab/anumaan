# profile_plots.R
# Post-fit visualization for the joint Bernoulli resistance-profile
# pipeline: observed-vs-model calibration forest plots and profile/
# pairwise heatmaps. Sampler diagnostic plots (trace/rank/energy) live in
# profile_check_sampler() instead -- those are about whether the SAMPLER
# worked, not what the fitted profiles look like, and are built from raw
# draws rather than profile_ppc()/profile_aggregate_for_daly() output.
#
# Deliberately matches this package's existing plotting conventions
# (plot_fixed_effect_diagnostics.R's forest-plot style,
# plot_omega_diagnostics.R's heatmap style) rather than inventing a new
# visual language: same geom_errorbar/geom_point forest-plot construction,
# same geom_tile/geom_text heatmap construction, same
# class_display_label()/hospital_display_label()/class_pair_label() use.


#' Observed-vs-Model Calibration Forest Plot
#'
#' One CI forest plot per \code{\link{profile_ppc}} level: the model's
#' posterior mean + credible interval (horizontal error bar, blue) beside
#' the observed value (a red x at the same row) -- a reader sees at a
#' glance whether the empirical point falls inside the model's interval,
#' exactly what \code{interval_contains_observed} encodes numerically.
#' \code{status == "not_computable"} rows are dropped, not drawn as
#' empty/misleading rows. Faceted by the entity (class, class pair, or
#' state) so many hospitals compare at a glance within one panel.
#'
#' @param ppc_table One of \code{\link{profile_ppc}}'s three tables
#'   (\code{$marginal}, \code{$pairwise}, or \code{$complete_profile}).
#' @param level Character, one of \code{"marginal"}, \code{"pairwise"},
#'   \code{"complete_profile"} -- which table \code{ppc_table} is, since
#'   the observed/model/entity column names differ by level.
#' @param title_base Character. Prefixed to the plot title.
#' @param max_facets Integer or \code{NULL}. For \code{"complete_profile"}
#'   (up to \code{2^D} facets, unreadable fast for \code{D > 4}), keeps
#'   only the \code{max_facets} states with the largest mean model
#'   probability across groups. \code{NULL}: keep all. Ignored for
#'   \code{"marginal"}/\code{"pairwise"}. Default \code{12L}.
#'
#' @return A ggplot object, or \code{NULL} if nothing in \code{ppc_table}
#'   was evaluated.
#' @export
profile_plot_calibration_forest <- function(ppc_table,
                                             level = c("marginal", "pairwise", "complete_profile"),
                                             title_base = "",
                                             max_facets = 12L) {
  level <- match.arg(level)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }

  df <- ppc_table[ppc_table$status == "evaluated", , drop = FALSE]
  if (nrow(df) == 0L) {
    return(NULL)
  }
  n_not_computable <- sum(ppc_table$status == "not_computable")

  if (level == "marginal") {
    df$entity <- class_display_label(df$class)
    df$observed <- df$observed_resistance
  } else if (level == "pairwise") {
    df$entity <- class_pair_label(df$class_a, df$class_b)
    df$observed <- df$observed_RR_prop
    df$model_mean <- df$model_RR_mean
    df$model_lower <- df$model_RR_lower
    df$model_upper <- df$model_RR_upper
  } else {
    df$entity <- df$state
    df$observed <- df$empirical_frequency
    if (!is.null(max_facets) && length(unique(df$state)) > max_facets) {
      state_means <- tapply(df$model_mean, df$state, mean, na.rm = TRUE)
      top_states <- names(sort(state_means, decreasing = TRUE))[seq_len(max_facets)]
      df <- df[df$state %in% top_states, , drop = FALSE]
    }
  }

  df$group_display <- hospital_display_label(df$group)
  df$group_display <- factor(df$group_display, levels = sort(unique(df$group_display), decreasing = TRUE))

  ggplot2::ggplot(df, ggplot2::aes(y = .data$group_display)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$model_lower, xmax = .data$model_upper),
      height = 0, colour = "#4292C6", alpha = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(x = .data$model_mean, shape = "Model (posterior mean)"),
      colour = "#2166AC", size = 1.6
    ) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed, shape = "Observed"),
      colour = "#B2182B", size = 1.8
    ) +
    ggplot2::scale_shape_manual(name = NULL, values = c("Model (posterior mean)" = 16, "Observed" = 4)) +
    ggplot2::facet_wrap(~entity) +
    ggplot2::labs(
      title = paste(
        title_base, "--",
        switch(level, marginal = "Marginal", pairwise = "Pairwise", complete_profile = "Complete-Profile"),
        "Calibration"
      ),
      subtitle = sprintf(
        "Model: posterior mean + 95%% CI (blue). Observed: red x. %d not-computable cell(s) omitted.",
        n_not_computable
      ),
      x = "Probability", y = "Hospital"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11), strip.text = ggplot2::element_text(size = 8))
}


#' Hospital x Profile-State Heatmap
#'
#' Visualizes the estimated or observed resistance-profile distribution
#' across hospitals and states in one grid -- columns are hospitals, rows
#' are states (profiles), fill is the selected probability column.
#' Mirrors \code{\link{plot_omega_correlation_heatmap}}'s visual style.
#'
#' @param profile_table Either \code{\link{profile_ppc}}'s
#'   \code{$complete_profile} table or
#'   \code{\link{profile_aggregate_for_daly}}'s output.
#' @param value_col Character. Which column to plot as the fill --
#'   \code{"model_mean"} or \code{"empirical_frequency"} for a
#'   \code{profile_ppc()$complete_profile} table, or \code{"R_mean"} for a
#'   \code{\link{profile_aggregate_for_daly}} table.
#' @param title_base Character. Prefixed to the plot title.
#' @param max_states Integer or \code{NULL}. Keeps only the
#'   \code{max_states} states with the largest mean value across
#'   hospitals (for readability when \code{D} is large). Default
#'   \code{16L}.
#'
#' @return A ggplot object, or \code{NULL} if nothing to plot.
#' @export
profile_plot_profile_heatmap <- function(profile_table, value_col, title_base = "", max_states = 16L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  state_col <- if ("state" %in% names(profile_table)) "state" else "profile_delta"
  df <- profile_table[!is.na(profile_table[[value_col]]), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(NULL)
  }

  if (!is.null(max_states) && length(unique(df[[state_col]])) > max_states) {
    state_means <- tapply(df[[value_col]], df[[state_col]], mean, na.rm = TRUE)
    top_states <- names(sort(state_means, decreasing = TRUE))[seq_len(max_states)]
    df <- df[df[[state_col]] %in% top_states, , drop = FALSE]
  }

  df$group_display <- hospital_display_label(df$group)
  df$label <- sprintf("%.2f", df[[value_col]])

  ggplot2::ggplot(df, ggplot2::aes(x = .data$group_display, y = .data[[state_col]], fill = .data[[value_col]])) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.8) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1), name = value_col) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste(title_base, "-- Resistance Profile Distribution by Hospital"),
      x = NULL, y = "Profile state"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)
    )
}


#' Pairwise Co-Resistance Heatmap, One Hospital
#'
#' Class x class heatmap of the model's pairwise P(R,R) for one hospital,
#' mirroring \code{\link{plot_omega_correlation_heatmap}}'s exact visual
#' style. This is a joint-Bernoulli pairwise probability, not a latent
#' residual correlation matrix -- a different quantity that happens to
#' render the same way (classes x classes, symmetric, one cell per pair).
#'
#' @param ppc_pairwise \code{\link{profile_ppc}}'s \code{$pairwise} table.
#' @param group Character. Which hospital to plot.
#' @param class_cols Character vector, all classes in the fit's canonical
#'   order (for consistent row/column ordering).
#' @param title_base Character. Prefixed to the plot title.
#'
#' @return A ggplot object, or \code{NULL} if this group has no evaluated
#'   pairwise cells.
#' @export
profile_plot_pairwise_heatmap <- function(ppc_pairwise, group, class_cols, title_base = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  sub <- ppc_pairwise[ppc_pairwise$group == group & ppc_pairwise$status == "evaluated", , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(NULL)
  }

  mirrored <- sub
  mirrored$class_a <- sub$class_b
  mirrored$class_b <- sub$class_a
  diag_tbl <- data.frame(
    class_a = class_cols, class_b = class_cols, model_RR_mean = NA_real_,
    stringsAsFactors = FALSE
  )
  full <- rbind(
    sub[c("class_a", "class_b", "model_RR_mean")],
    mirrored[c("class_a", "class_b", "model_RR_mean")],
    diag_tbl
  )
  full$class_a <- factor(class_display_label(full$class_a), levels = class_display_label(class_cols))
  full$class_b <- factor(class_display_label(full$class_b), levels = rev(class_display_label(class_cols)))
  full$label <- ifelse(is.na(full$model_RR_mean), "", sprintf("%.2f", full$model_RR_mean))

  ggplot2::ggplot(full, ggplot2::aes(x = .data$class_a, y = .data$class_b, fill = .data$model_RR_mean)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "#B2182B", limits = c(0, 1), na.value = "grey90",
      name = "P(R,R)"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste(title_base, "--", hospital_display_label(group), "Pairwise Co-Resistance"),
      subtitle = "Model posterior-mean P(both resistant) -- exact joint Bernoulli marginalization, not simulated.",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      plot.subtitle = ggplot2::element_text(size = 8, face = "italic"),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)
    )
}
