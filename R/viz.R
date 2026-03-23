# viz.R
# Generic visualization functions for AMR data analysis

#' AMR Theme for ggplot2
#'
#' Universal theme for all AMR plots with consistent styling.
#'
#' @param base_size Base font size. Default 12.
#' @param legend_position Legend position. Default "top".
#'
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   amr_theme()
amr_theme <- function(base_size = 16, legend_position = "top") {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      legend.key = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.title.y = ggplot2::element_text(size = base_size + 2, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 4, face = "bold"),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 2,
        color = "grey40"
      ),
      plot.background = ggplot2::element_rect(fill = "white"),
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(
        size = base_size - 1,
        face = "bold"
      )
    )
}


#' Get Color Palette
#'
#' Returns predefined color palettes for AMR visualizations.
#'
#' @param palette_name Character. Name of palette: "resistance", "outcome", "default", "green", "spectral".
#'
#' @return Named vector of colors
#' @export
get_amr_palette <- function(palette_name = "default") {
  palettes <- list(
    resistance = c("R" = "#F44336", "S" = "#4CAF50", "I" = "#FFC107"),
    outcome = c(
      "Death" = "#E67E22", "Died" = "#E67E22",
      "Discharged" = "#5D6D7E",
      "LAMA" = "#BDC3C7",
      "Transferred to other hospital" = "#58D68D",
      "Referred" = "#58D68D"
    ),
    green = c("#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32"),
    spectral = RColorBrewer::brewer.pal(11, "Spectral"),
    default = RColorBrewer::brewer.pal(8, "Set2"),
    # GBD-style burden palettes
    burden_gbd4 = c(
      "Associated with resistance" = "#b2abd2",
      "Attributable to resistance" = "#4d2d8a"
    ),
    burden_gbd6 = c(
      "#2166ac", "#4393c3", "#92c5de", "#f7f7f7",
      "#fddbc7", "#f4a582", "#d6604d", "#b2182b"
    ),
    burden_yll = c(associated = "#d73027", attributable = "#f4a582"),
    burden_yld = c(associated = "#4393c3", attributable = "#92c5de"),
    burden_daly = c(associated = "#762a83", attributable = "#f1a340"),
    burden_lollipop_attr = c(low = "#fc8d59", high = "#7f0000"),
    burden_lollipop_assoc = c(low = "#9ecae1", high = "#08306b"),
    burden_fraction = c(
      "Attributable to resistance" = "#d73027",
      "Baseline infectious burden" = "#4393c3"
    )
  )

  if (palette_name %in% names(palettes)) {
    return(palettes[[palette_name]])
  } else {
    return(palettes$default)
  }
}


#' Generic Bar Plot
#'
#' Creates a bar plot with flexible column mapping.
#'
#' @param data Data frame
#' @param x Character. Column name for x-axis (categorical variable)
#' @param y Character. Column name for y-axis (numeric variable). If NULL, uses count.
#' @param fill Character. Optional column name for fill color
#' @param facet Character. Optional column name for faceting
#' @param title Character. Plot title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param palette Character or named vector. Color palette name or custom colors
#' @param show_labels Logical. Show value labels on bars. Default TRUE.
#' @param flip_coords Logical. Flip coordinates (horizontal bars). Default FALSE.
#' @param top_n Numeric. Show only top N categories. Default NULL (show all).
#' @param order Character. Order bars by: "desc" (descending), "asc" (ascending), or "none". Default "desc".
#' @param text_size Numeric. Label text size. NULL for auto-adjustment. Default NULL.
#' @param text_angle Numeric. X-axis text angle (0-90). NULL for auto-adjustment. Default NULL.
#' @param bar_width Numeric. Bar width (0-1). Default 0.8.
#'
#' @return A ggplot object
#' @export
plot_bar <- function(data, x, y = NULL, fill = NULL, facet = NULL,
                     title = NULL, xlab = NULL, ylab = NULL,
                     palette = "default", show_labels = TRUE,
                     flip_coords = FALSE, top_n = NULL,
                     order = "desc", text_size = NULL, text_angle = NULL,
                     bar_width = 0.8) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Filter top N if specified
  if (!is.null(top_n)) {
    if (is.null(y)) {
      # Count mode
      top_categories <- data %>%
        dplyr::count(!!rlang::sym(x), sort = TRUE) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(!!rlang::sym(x))
    } else {
      # Value mode
      top_categories <- data %>%
        dplyr::group_by(!!rlang::sym(x)) %>%
        dplyr::summarise(total = sum(!!rlang::sym(y), na.rm = TRUE)) %>%
        dplyr::slice_max(total, n = top_n) %>%
        dplyr::pull(!!rlang::sym(x))
    }
    data <- data %>% dplyr::filter(!!rlang::sym(x) %in% top_categories)
  }

  # Order categories if requested
  if (order != "none") {
    if (is.null(y)) {
      # Count mode - order by count
      category_order <- data %>%
        dplyr::count(!!rlang::sym(x), sort = TRUE) %>%
        dplyr::pull(!!rlang::sym(x))
    } else {
      # Value mode - order by sum of y
      category_order <- data %>%
        dplyr::group_by(!!rlang::sym(x)) %>%
        dplyr::summarise(total = sum(!!rlang::sym(y), na.rm = TRUE)) %>%
        dplyr::arrange(dplyr::desc(total)) %>%
        dplyr::pull(!!rlang::sym(x))
    }

    # Reverse for ascending
    if (order == "asc") {
      category_order <- rev(category_order)
    }

    # Reorder factor levels
    data[[x]] <- factor(data[[x]], levels = category_order)
  }

  # Auto-adjust text size based on number of categories
  if (is.null(text_size)) {
    n_categories <- length(unique(data[[x]]))
    text_size <- ifelse(n_categories > 20, 2.5,
      ifelse(n_categories > 10, 3, 3.5)
    )
  }

  # Build plot
  if (is.null(y)) {
    # Count mode
    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x)))
    if (!is.null(fill)) {
      p <- p + ggplot2::aes(fill = !!rlang::sym(fill))
    }
    p <- p + ggplot2::geom_bar(color = "black", width = bar_width)
  } else {
    # Value mode
    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(y)))
    if (!is.null(fill)) {
      p <- p + ggplot2::aes(fill = !!rlang::sym(fill))
    }
    p <- p + ggplot2::geom_col(color = "black", width = bar_width)
  }

  # Add labels
  if (show_labels) {
    if (is.null(y)) {
      p <- p + ggplot2::geom_text(
        stat = "count",
        ggplot2::aes(label = ggplot2::after_stat(count)),
        vjust = -0.3, fontface = "bold", size = text_size
      )
    } else {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = round(!!rlang::sym(y), 1)),
        vjust = -0.3, fontface = "bold", size = text_size
      )
    }
  }

  # Colors
  if (!is.null(fill)) {
    if (is.character(palette) && length(palette) == 1) {
      colors <- get_amr_palette(palette)

      # Check if we have enough colors
      n_categories <- length(unique(data[[fill]]))

      if (n_categories <= length(colors)) {
        # Use predefined palette
        p <- p + ggplot2::scale_fill_manual(values = colors)
      } else {
        # Use a color palette that can handle any number of categories
        message(sprintf(
          "Using automatic colors for %d categories (palette only has %d)",
          n_categories, length(colors)
        ))
        p <- p + ggplot2::scale_fill_hue()
      }
    } else if (length(palette) > 1) {
      p <- p + ggplot2::scale_fill_manual(values = palette)
    }
  }

  # Faceting
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free_y")
  }

  # Flip
  if (flip_coords) {
    p <- p + ggplot2::coord_flip()
  }

  # Auto-adjust text angle based on number of categories and flip status
  if (is.null(text_angle)) {
    n_categories <- length(unique(data[[x]]))
    if (!flip_coords) {
      # Vertical bars - angle x-axis text
      text_angle <- ifelse(n_categories > 10, 45,
        ifelse(n_categories > 5, 30, 0)
      )
    } else {
      # Horizontal bars - no angle needed
      text_angle <- 0
    }
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    amr_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = text_angle, hjust = ifelse(text_angle > 0, 1, 0.5))
    )

  return(p)
}


#' Stacked Bar Plot
#'
#' Creates a stacked bar plot showing composition.
#'
#' @param data Data frame
#' @param x Character. Column name for x-axis
#' @param fill Character. Column name for stacking/fill variable
#' @param y Character. Optional column name for y-axis. If NULL, uses count.
#' @param position Character. "stack" (default), "fill" (100 percent stacked), or "dodge" (grouped)
#' @param facet Character. Optional column name for faceting
#' @param title Character. Plot title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param palette Character or named vector. Color palette
#' @param show_labels Logical. Show value labels. Default TRUE.
#' @param show_percentages Logical. Show percentages in labels. Default TRUE.
#' @param flip_coords Logical. Flip coordinates. Default FALSE.
#' @param order Character. Order bars by: "desc" (descending), "asc" (ascending), or "none". Default "desc".
#' @param text_size Numeric. Label text size. NULL for auto-adjustment. Default NULL.
#' @param text_angle Numeric. X-axis text angle (0-90). NULL for auto-adjustment. Default NULL.
#' @param bar_width Numeric. Bar width (0-1). Default 0.7.
#'
#' @return A ggplot object
#' @export
plot_stacked_bar <- function(data, x, fill, y = NULL, position = "stack",
                             facet = NULL, title = NULL, xlab = NULL, ylab = NULL,
                             palette = "default", show_labels = TRUE,
                             show_percentages = TRUE, flip_coords = FALSE,
                             order = "desc", text_size = NULL, text_angle = NULL,
                             bar_width = 0.7) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Order categories if requested
  if (order != "none") {
    if (is.null(y)) {
      category_order <- data %>%
        dplyr::count(!!rlang::sym(x), sort = TRUE) %>%
        dplyr::pull(!!rlang::sym(x))
    } else {
      category_order <- data %>%
        dplyr::group_by(!!rlang::sym(x)) %>%
        dplyr::summarise(total = sum(!!rlang::sym(y), na.rm = TRUE)) %>%
        dplyr::arrange(dplyr::desc(total)) %>%
        dplyr::pull(!!rlang::sym(x))
    }

    # Reverse for ascending
    if (order == "asc") {
      category_order <- rev(category_order)
    }

    # Reorder factor levels
    data[[x]] <- factor(data[[x]], levels = category_order)
  }

  # Auto-adjust text size based on number of categories
  if (is.null(text_size)) {
    n_categories <- length(unique(data[[x]]))
    text_size <- ifelse(n_categories > 20, 2.5,
      ifelse(n_categories > 10, 3, 3.5)
    )
  }

  # Calculate proportions for labels
  if (show_percentages) {
    if (is.null(y)) {
      label_data <- data %>%
        dplyr::group_by(!!rlang::sym(x), !!rlang::sym(fill)) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(!!rlang::sym(x)) %>%
        dplyr::mutate(
          total = sum(count),
          percentage = round(count / total * 100, 1)
        )
    } else {
      label_data <- data %>%
        dplyr::group_by(!!rlang::sym(x), !!rlang::sym(fill)) %>%
        dplyr::summarise(value = sum(!!rlang::sym(y), na.rm = TRUE), .groups = "drop") %>%
        dplyr::group_by(!!rlang::sym(x)) %>%
        dplyr::mutate(
          total = sum(value),
          percentage = round(value / total * 100, 1)
        )
    }
  }

  # Build plot
  if (is.null(y)) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x), fill = !!rlang::sym(fill))) +
      ggplot2::geom_bar(position = position, color = "black", width = bar_width)
  } else {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(y), fill = !!rlang::sym(fill))) +
      ggplot2::geom_col(position = position, color = "black", width = bar_width)
  }

  # Add labels
  if (show_labels && show_percentages) {
    if (is.null(y)) {
      label_text <- paste0(label_data$count, " (", label_data$percentage, "%)")
    } else {
      label_text <- paste0(round(label_data$value, 1), " (", label_data$percentage, "%)")
    }

    p <- p + ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(label = label_text, y = if (is.null(y)) count else value),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = text_size
    )
  } else if (show_labels) {
    if (is.null(y)) {
      p <- p + ggplot2::geom_text(
        stat = "count",
        ggplot2::aes(label = ggplot2::after_stat(count)),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white", fontface = "bold", size = text_size
      )
    } else {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = round(!!rlang::sym(y), 1)),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white", fontface = "bold", size = text_size
      )
    }
  }

  # Colors
  if (is.character(palette) && length(palette) == 1) {
    colors <- get_amr_palette(palette)
    p <- p + ggplot2::scale_fill_manual(values = colors)
  } else if (length(palette) > 1) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }

  # Percentage axis for fill position
  if (position == "fill") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent_format())
  }

  # Faceting
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free")
  }

  # Flip
  if (flip_coords) {
    p <- p + ggplot2::coord_flip()
  }

  # Auto-adjust text angle based on number of categories and flip status
  if (is.null(text_angle)) {
    n_categories <- length(unique(data[[x]]))
    if (!flip_coords) {
      # Vertical bars - angle x-axis text
      text_angle <- ifelse(n_categories > 10, 45,
        ifelse(n_categories > 5, 30, 0)
      )
    } else {
      # Horizontal bars - no angle needed
      text_angle <- 0
    }
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    amr_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = text_angle, hjust = ifelse(text_angle > 0, 1, 0.5))
    )

  return(p)
}


#' Line Chart
#'
#' Creates a line chart for trends over time or continuous variables.
#'
#' @param data Data frame
#' @param x Character. Column name for x-axis (usually time/date)
#' @param y Character. Column name for y-axis
#' @param group Character. Optional column name for grouping lines
#' @param color Character. Optional column name for line color
#' @param facet Character. Optional column name for faceting
#' @param title Character. Plot title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param palette Character or named vector. Color palette
#' @param show_points Logical. Show points on lines. Default TRUE.
#' @param show_labels Logical. Show value labels. Default FALSE.
#'
#' @return A ggplot object
#' @export
plot_line <- function(data, x, y, group = NULL, color = NULL, facet = NULL,
                      title = NULL, xlab = NULL, ylab = NULL,
                      palette = "default", show_points = TRUE, show_labels = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Build plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(y)))

  # Add grouping
  if (!is.null(group)) {
    p <- p + ggplot2::aes(group = !!rlang::sym(group))
  }

  # Add color
  if (!is.null(color)) {
    p <- p + ggplot2::aes(color = !!rlang::sym(color))
  }

  # Line
  p <- p + ggplot2::geom_line(size = 1.2)

  # Points
  if (show_points) {
    p <- p + ggplot2::geom_point(size = 3)
  }

  # Labels
  if (show_labels) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = round(!!rlang::sym(y), 1)),
      vjust = -0.5, fontface = "bold"
    )
  }

  # Colors
  if (!is.null(color)) {
    if (is.character(palette) && length(palette) == 1) {
      colors <- get_amr_palette(palette)
      p <- p + ggplot2::scale_color_manual(values = colors)
    } else if (length(palette) > 1) {
      p <- p + ggplot2::scale_color_manual(values = palette)
    }
  }

  # Faceting
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free")
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    amr_theme()

  return(p)
}


#' Histogram
#'
#' Creates a histogram for continuous variables.
#'
#' @param data Data frame
#' @param x Character. Column name for continuous variable
#' @param fill Character. Optional column name for fill color
#' @param facet Character. Optional column name for faceting
#' @param bins Numeric. Number of bins. Default 30.
#' @param title Character. Plot title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param palette Character or named vector. Color palette
#'
#' @return A ggplot object
#' @export
plot_histogram <- function(data, x, fill = NULL, facet = NULL, bins = 30,
                           title = NULL, xlab = NULL, ylab = "Count",
                           palette = "default") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Build plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x)))

  if (!is.null(fill)) {
    p <- p + ggplot2::aes(fill = !!rlang::sym(fill)) +
      ggplot2::geom_histogram(bins = bins, color = "black", alpha = 0.7)
  } else {
    p <- p + ggplot2::geom_histogram(bins = bins, fill = "steelblue", color = "black")
  }

  # Colors
  if (!is.null(fill)) {
    if (is.character(palette) && length(palette) == 1) {
      colors <- get_amr_palette(palette)
      p <- p + ggplot2::scale_fill_manual(values = colors)
    } else if (length(palette) > 1) {
      p <- p + ggplot2::scale_fill_manual(values = palette)
    }
  }

  # Faceting
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free")
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    amr_theme()

  return(p)
}


#' Heatmap
#'
#' Creates a heatmap (tile plot) for matrix data.
#'
#' @param data Data frame
#' @param x Character. Column name for x-axis
#' @param y Character. Column name for y-axis
#' @param fill Character. Column name for fill value
#' @param title Character. Plot title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param palette Character or vector. Color palette. Default "green".
#' @param show_labels Logical. Show value labels in tiles. Default TRUE.
#'
#' @return A ggplot object
#' @export
plot_heatmap <- function(data, x, y, fill, title = NULL,
                         xlab = NULL, ylab = NULL,
                         palette = "green", show_labels = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Build plot
  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = !!rlang::sym(x),
    y = !!rlang::sym(y),
    fill = !!rlang::sym(fill)
  )) +
    ggplot2::geom_tile(color = "black")

  # Labels
  if (show_labels) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = round(!!rlang::sym(fill), 2)),
      color = "black", fontface = "bold", size = 3.5
    )
  }

  # Colors
  if (is.character(palette) && length(palette) == 1) {
    colors <- get_amr_palette(palette)
    p <- p + ggplot2::scale_fill_gradientn(colours = colors, na.value = "grey90")
  } else if (length(palette) > 1) {
    p <- p + ggplot2::scale_fill_gradientn(colours = palette, na.value = "grey90")
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  return(p)
}


#' Resistance Pattern Heatmap
#'
#' Creates a heatmap showing resistance patterns across isolates and
#' antibiotic classes. Shows S (susceptible), R (resistant), Partial (mixed),
#' or NT (not tested) for each isolate-class combination.
#'
#' @param data Data frame with isolate/event data
#' @param isolate_col Character. Column name for isolate/event ID
#' @param class_col Character. Column name for antibiotic class
#' @param result_col Character. Column name for resistance result (S/R/I/NT)
#' @param top_n Numeric. Show only top N isolates. Default NULL (show all).
#' @param title Character. Plot title
#' @param show_labels Logical. Show result labels in tiles. Default FALSE.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # Resistance heatmap
#' plot_resistance_heatmap(data,
#'   isolate_col = "event_id",
#'   class_col = "antibiotic_class",
#'   result_col = "resistance_status",
#'   top_n = 20
#' )
plot_resistance_heatmap <- function(data, isolate_col, class_col, result_col,
                                    top_n = NULL, title = NULL, show_labels = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Filter top N isolates if specified
  if (!is.null(top_n)) {
    top_isolates <- data %>%
      dplyr::distinct(!!rlang::sym(isolate_col)) %>%
      dplyr::slice_head(n = top_n)
    data <- data %>%
      dplyr::filter(!!rlang::sym(isolate_col) %in% top_isolates[[isolate_col]])
  }

  # Summarize per isolate x class
  class_summary <- data %>%
    dplyr::mutate(
      result_clean = toupper(trimws(!!rlang::sym(result_col)))
    ) %>%
    dplyr::group_by(!!rlang::sym(isolate_col), !!rlang::sym(class_col)) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_tested = sum(result_clean %in% c("S", "R", "I")),
      n_s = sum(result_clean == "S"),
      n_r = sum(result_clean == "R"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      class_result = dplyr::case_when(
        n_tested == 0 ~ "NT", # Nothing tested
        n_s == n_tested ~ "S", # All susceptible
        n_r == n_tested ~ "R", # All resistant
        TRUE ~ "Partial" # Mixed
      ),
      class_result = factor(class_result, levels = c("S", "R", "Partial", "NT"))
    )

  # Plot
  p <- ggplot2::ggplot(class_summary, ggplot2::aes(
    x = !!rlang::sym(class_col),
    y = factor(!!rlang::sym(isolate_col)),
    fill = class_result
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = c(
        "S" = "#2ECC71", # Green
        "R" = "#E74C3C", # Red
        "Partial" = "#F39C12", # Orange
        "NT" = "grey80" # Grey
      )
    ) +
    ggplot2::labs(
      x = "",
      y = "Isolate",
      fill = "Result",
      title = if (is.null(title)) "Resistance Pattern Heatmap" else title
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 9),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  # Add labels if requested
  if (show_labels) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = class_result),
      color = "black",
      size = 3
    )
  }

  return(p)
}


#' Proportion Bar Plot
#'
#' Creates a 100% stacked bar plot showing proportions with percentage labels.
#'
#' @param data Data frame
#' @param x Character. Column name for x-axis
#' @param fill Character. Column name for fill/stacking variable
#' @param facet Character. Optional column name for faceting
#' @param title Character. Plot title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param palette Character or named vector. Color palette
#' @param show_counts Logical. Show counts alongside percentages. Default TRUE.
#' @param show_totals Logical. Show total counts above bars. Default TRUE.
#' @param flip_coords Logical. Flip coordinates. Default FALSE.
#' @param order Character. Order bars by: "desc" (descending), "asc" (ascending), or "none". Default "desc".
#' @param text_size Numeric. Label text size. NULL for auto-adjustment. Default NULL.
#' @param text_angle Numeric. X-axis text angle (0-90). NULL for auto-adjustment. Default NULL.
#' @param bar_width Numeric. Bar width (0-1). Default 0.7.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # 100% proportion plot
#' plot_proportion(data,
#'   x = "organism", fill = "resistance_status",
#'   palette = "resistance"
#' )
plot_proportion <- function(data, x, fill, facet = NULL,
                            title = NULL, xlab = NULL, ylab = "Proportion",
                            palette = "default", show_counts = TRUE,
                            show_totals = TRUE, flip_coords = FALSE,
                            order = "desc", text_size = NULL, text_angle = NULL,
                            bar_width = 0.7) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Calculate proportions
  prop_data <- data %>%
    dplyr::group_by(!!rlang::sym(x), !!rlang::sym(fill)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(!!rlang::sym(x)) %>%
    dplyr::mutate(
      total = sum(count),
      proportion = count / total,
      percentage = round(proportion * 100, 1)
    ) %>%
    dplyr::ungroup()

  # Order categories if requested
  if (order != "none") {
    category_order <- prop_data %>%
      dplyr::group_by(!!rlang::sym(x)) %>%
      dplyr::summarise(total = sum(count)) %>%
      dplyr::arrange(dplyr::desc(total)) %>%
      dplyr::pull(!!rlang::sym(x))

    # Reverse for ascending
    if (order == "asc") {
      category_order <- rev(category_order)
    }

    # Reorder factor levels
    prop_data[[x]] <- factor(prop_data[[x]], levels = category_order)
  }

  # Auto-adjust text size based on number of categories
  if (is.null(text_size)) {
    n_categories <- length(unique(prop_data[[x]]))
    text_size <- ifelse(n_categories > 20, 2.5,
      ifelse(n_categories > 10, 3, 3.5)
    )
  }

  # Build plot
  p <- ggplot2::ggplot(prop_data, ggplot2::aes(
    x = !!rlang::sym(x),
    y = proportion,
    fill = !!rlang::sym(fill)
  )) +
    ggplot2::geom_col(width = bar_width, color = "black")

  # Add labels
  if (show_counts) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = paste0(count, "\n(", percentage, "%)")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = text_size,
      color = "white",
      fontface = "bold"
    )
  } else {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = paste0(percentage, "%")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = text_size,
      color = "white",
      fontface = "bold"
    )
  }

  # Add totals above bars
  if (show_totals) {
    total_labels <- prop_data %>%
      dplyr::distinct(!!rlang::sym(x), total)

    p <- p + ggplot2::geom_text(
      data = total_labels,
      ggplot2::aes(x = !!rlang::sym(x), y = 1.05, label = paste0("n=", total)),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    )
  }

  # Colors
  if (is.character(palette) && length(palette) == 1) {
    colors <- get_amr_palette(palette)
    p <- p + ggplot2::scale_fill_manual(values = colors)
  } else if (length(palette) > 1) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }

  # Percentage y-axis
  p <- p + ggplot2::scale_y_continuous(
    labels = scales::percent_format(),
    expand = ggplot2::expansion(mult = c(0, 0.1))
  )

  # Faceting
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free_x")
  }

  # Flip
  if (flip_coords) {
    p <- p + ggplot2::coord_flip()
  }

  # Auto-adjust text angle based on number of categories and flip status
  if (is.null(text_angle)) {
    n_categories <- length(unique(prop_data[[x]]))
    if (!flip_coords) {
      # Vertical bars - angle x-axis text
      text_angle <- ifelse(n_categories > 10, 45,
        ifelse(n_categories > 5, 30, 0)
      )
    } else {
      # Horizontal bars - no angle needed
      text_angle <- 0
    }
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    amr_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = text_angle, hjust = ifelse(text_angle > 0, 1, 0.5))
    )

  return(p)
}


#' Grouped (Dodged) Bar Plot
#'
#' Creates a grouped bar plot for side-by-side comparisons.
#'
#' @param data Data frame
#' @param x Character. Column name for x-axis
#' @param fill Character. Column name for grouping variable
#' @param y Character. Optional column name for y-axis. If NULL, uses count.
#' @param facet Character. Optional column name for faceting
#' @param title Character. Plot title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param palette Character or named vector. Color palette
#' @param show_labels Logical. Show value labels. Default TRUE.
#' @param flip_coords Logical. Flip coordinates. Default FALSE.
#' @param order Character. Order bars by: "desc" (descending), "asc" (ascending), or "none". Default "desc".
#' @param text_size Numeric. Label text size. NULL for auto-adjustment. Default NULL.
#' @param text_angle Numeric. X-axis text angle (0-90). NULL for auto-adjustment. Default NULL.
#' @param bar_width Numeric. Bar width (0-1). Default 0.7.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # Grouped bar for comparing MRSA vs MDR
#' plot_grouped_bar(data, x = "hospital", fill = "mrsa_status")
plot_grouped_bar <- function(data, x, fill, y = NULL, facet = NULL,
                             title = NULL, xlab = NULL, ylab = NULL,
                             palette = "default", show_labels = TRUE,
                             flip_coords = FALSE,
                             order = "desc", text_size = NULL, text_angle = NULL,
                             bar_width = 0.7) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Order categories if requested
  if (order != "none") {
    if (is.null(y)) {
      category_order <- data %>%
        dplyr::count(!!rlang::sym(x), sort = TRUE) %>%
        dplyr::pull(!!rlang::sym(x))
    } else {
      category_order <- data %>%
        dplyr::group_by(!!rlang::sym(x)) %>%
        dplyr::summarise(total = sum(!!rlang::sym(y), na.rm = TRUE)) %>%
        dplyr::arrange(dplyr::desc(total)) %>%
        dplyr::pull(!!rlang::sym(x))
    }

    # Reverse for ascending
    if (order == "asc") {
      category_order <- rev(category_order)
    }

    # Reorder factor levels
    data[[x]] <- factor(data[[x]], levels = category_order)
  }

  # Auto-adjust text size based on number of categories
  if (is.null(text_size)) {
    n_categories <- length(unique(data[[x]]))
    text_size <- ifelse(n_categories > 20, 2.5,
      ifelse(n_categories > 10, 3, 3.5)
    )
  }

  # Build plot
  if (is.null(y)) {
    # Count mode
    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x), fill = !!rlang::sym(fill))) +
      ggplot2::geom_bar(
        position = ggplot2::position_dodge(width = 0.8),
        width = bar_width,
        color = "black"
      )
  } else {
    # Value mode
    p <- ggplot2::ggplot(data, ggplot2::aes(
      x = !!rlang::sym(x),
      y = !!rlang::sym(y),
      fill = !!rlang::sym(fill)
    )) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = 0.8),
        width = bar_width,
        color = "black"
      )
  }

  # Add labels
  if (show_labels) {
    if (is.null(y)) {
      p <- p + ggplot2::geom_text(
        stat = "count",
        ggplot2::aes(label = ggplot2::after_stat(count)),
        position = ggplot2::position_dodge(width = 0.8),
        vjust = -0.3,
        fontface = "bold",
        size = text_size
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = round(!!rlang::sym(y), 1)),
        position = ggplot2::position_dodge(width = 0.8),
        vjust = -0.3,
        fontface = "bold",
        size = text_size
      )
    }
  }

  # Colors
  if (is.character(palette) && length(palette) == 1) {
    colors <- get_amr_palette(palette)
    p <- p + ggplot2::scale_fill_manual(values = colors)
  } else if (length(palette) > 1) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }

  # Faceting
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free")
  }

  # Flip
  if (flip_coords) {
    p <- p + ggplot2::coord_flip()
  }

  # Auto-adjust text angle based on number of categories and flip status
  if (is.null(text_angle)) {
    n_categories <- length(unique(data[[x]]))
    if (!flip_coords) {
      # Vertical bars - angle x-axis text
      text_angle <- ifelse(n_categories > 10, 45,
        ifelse(n_categories > 5, 30, 0)
      )
    } else {
      # Horizontal bars - no angle needed
      text_angle <- 0
    }
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    amr_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = text_angle, hjust = ifelse(text_angle > 0, 1, 0.5))
    )

  return(p)
}
