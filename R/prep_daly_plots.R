# prep_daly_plots.R
# DALY burden computation and plotting functions (YLL / YLD)


# COMPUTATION

#' Compute Hospital-Level DALY Breakdown
#'
#' Distributes pooled YLL/YLD totals to individual hospitals proportionally
#' by deaths and discharges, and computes per-1 000 rates.
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
#'   per 1 000 cases) columns.
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
      .death_frac       = deaths_h / total_deaths,
      .disch_frac       = discharged_h / total_discharged,
      .per_1k           = 1000 / cases_h,
      # Absolute
      YLL_base          = yll_base          * .data$.death_frac,
      YLD_base          = yld_base          * .data$.disch_frac,
      YLL_associated    = yll_associated    * .data$.death_frac,
      YLD_associated    = yld_associated    * .data$.disch_frac,
      DALY_associated   = .data$YLL_associated + .data$YLD_associated,
      YLL_attributable  = yll_attributable  * .data$.death_frac,
      YLD_attributable  = yld_attributable  * .data$.disch_frac,
      DALY_attributable = .data$YLL_attributable + .data$YLD_attributable,
      # Per 1 000 cases
      YLL_base_per_1000   = .data$YLL_base          * .data$.per_1k,
      YLD_base_per_1000   = .data$YLD_base          * .data$.per_1k,
      YLL_assoc_per_1000  = .data$YLL_associated    * .data$.per_1k,
      YLD_assoc_per_1000  = .data$YLD_associated    * .data$.per_1k,
      DALY_assoc_per_1000 = .data$DALY_associated   * .data$.per_1k,
      YLL_attr_per_1000   = .data$YLL_attributable  * .data$.per_1k,
      YLD_attr_per_1000   = .data$YLD_attributable  * .data$.per_1k,
      DALY_attr_per_1000  = .data$DALY_attributable * .data$.per_1k
    ) %>%
    dplyr::select(-".death_frac", -".disch_frac", -".per_1k") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4)))
}


# YLL / YLD ASSOCIATED vs ATTRIBUTABLE \u2014 GROUPED BAR BY HOSPITAL


#' Plot YLL or YLD Associated vs Attributable per Hospital
#'
#' Produces a horizontal grouped bar chart comparing the associated and
#' attributable burden (YLL or YLD) per 1 000 patients across hospitals.
#' Bars are sorted by the associated value (largest at top by default).
#' Each bar is labelled with its rounded value.
#'
#' The input \code{data} should be a hospital-level summary data frame
#' (e.g. the output of \code{compute_hospital_daly()}) that already contains
#' the per-1 000 columns. If the summary was computed across multiple
#' syndromes, use \code{syndrome_col} and \code{syndrome_name} to restrict
#' the plot to one syndrome before pivoting.
#'
#' @param data         Data frame. One row per hospital (or one row per
#'   hospital x syndrome). Must contain the columns named by
#'   \code{center_col}, \code{assoc_col}, and \code{attr_col}.
#' @param metric       Character. \code{"YLL"} or \code{"YLD"}. Controls
#'   which default column names and axis labels are used. Ignored when
#'   \code{assoc_col} and \code{attr_col} are supplied explicitly.
#'   Default \code{"YLL"}.
#' @param center_col   Character. Column containing hospital / centre names.
#'   Default \code{"center_name"}.
#' @param assoc_col    Character or \code{NULL}. Column for the associated
#'   per-1 000 value. If \code{NULL} (default), auto-set to
#'   \code{"YLL_assoc_per_1000"} or \code{"YLD_assoc_per_1000"} based on
#'   \code{metric}.
#' @param attr_col     Character or \code{NULL}. Column for the attributable
#'   per-1 000 value. If \code{NULL} (default), auto-set to
#'   \code{"YLL_attr_per_1000"} or \code{"YLD_attr_per_1000"} based on
#'   \code{metric}.
#' @param syndrome_col  Character or \code{NULL}. Column name containing
#'   syndrome labels (e.g. \code{"infectious_syndrome"}). If \code{NULL}
#'   (default), no syndrome filtering is applied.
#' @param syndrome_name Character or \code{NULL}. The syndrome value to
#'   retain (e.g. \code{"Bloodstream infection"}). Requires
#'   \code{syndrome_col} to be set. If \code{NULL} (default), all rows
#'   are used.
#' @param sort_by      Character. \code{"associated"} (default),
#'   \code{"attributable"}, or \code{"none"}. Controls hospital ordering
#'   on the y-axis.
#' @param label_digits Integer. Decimal places in bar value labels. Default 1.
#' @param colours      Named character vector with exactly two elements:
#'   \code{"associated"} and \code{"attributable"}. If \code{NULL}
#'   (default), a colour pair is chosen automatically based on
#'   \code{metric} (red tones for YLL, blue tones for YLD).
#' @param base_size    Numeric. Base font size. Default 14.
#' @param title        Character. Custom plot title. Auto-generated if
#'   \code{NULL}.
#' @param subtitle     Character. Custom subtitle. Auto-generated if
#'   \code{NULL}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
plot_burden_by_hospital <- function(data,
                                    metric        = c("YLL", "YLD"),
                                    center_col    = "center_name",
                                    assoc_col     = NULL,
                                    attr_col      = NULL,
                                    syndrome_col  = NULL,
                                    syndrome_name = NULL,
                                    sort_by       = c("associated", "attributable", "none"),
                                    label_digits  = 1,
                                    colours       = NULL,
                                    base_size     = 14,
                                    title         = NULL,
                                    subtitle      = NULL) {

  # -- 0. match args -----------------------------------------------------------
  metric  <- match.arg(metric)
  sort_by <- match.arg(sort_by)

  # -- 1. resolve column names -------------------------------------------------
  if (is.null(assoc_col))
    assoc_col <- if (metric == "YLL") "YLL_assoc_per_1000" else "YLD_assoc_per_1000"
  if (is.null(attr_col))
    attr_col  <- if (metric == "YLL") "YLL_attr_per_1000"  else "YLD_attr_per_1000"

  # -- 2. validate columns -----------------------------------------------------
  required_cols <- c(center_col, assoc_col, attr_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s",
                 paste(missing_cols, collapse = ", ")))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- 3. syndrome filter ------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 4. default colours ------------------------------------------------------
  if (is.null(colours)) {
    colours <- if (metric == "YLL") {
      c("associated" = "#d73027", "attributable" = "#f4a582")
    } else {
      c("associated" = "#4393c3", "attributable" = "#92c5de")
    }
  }
  if (!all(c("associated", "attributable") %in% names(colours)))
    stop("'colours' must be a named vector with elements 'associated' and 'attributable'.")

  # -- 5. legend label strings -------------------------------------------------
  assoc_label <- paste(metric, "Associated")
  attr_label  <- paste(metric, "Attributable")
  y_label     <- sprintf("%s per 1 000 admissions", metric)

  # -- 6. sort hospitals -------------------------------------------------------
  ctr_sym   <- rlang::sym(center_col)
  sort_col  <- if (sort_by == "attributable") attr_col else assoc_col

  if (sort_by != "none") {
    ordered_levels <- data %>%
      dplyr::arrange(!!rlang::sym(sort_col)) %>%
      dplyr::pull(!!ctr_sym)
  } else {
    ordered_levels <- data %>% dplyr::pull(!!ctr_sym)
  }

  # -- 7. reshape to long format -----------------------------------------------
  plot_df <- data %>%
    dplyr::select(
      center       = !!ctr_sym,
      associated   = !!rlang::sym(assoc_col),
      attributable = !!rlang::sym(attr_col)
    ) %>%
    tidyr::pivot_longer(
      cols      = c(associated, attributable),
      names_to  = "metric_type",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      center      = factor(center, levels = ordered_levels),
      metric_type = factor(
        dplyr::case_when(
          metric_type == "associated"   ~ assoc_label,
          metric_type == "attributable" ~ attr_label
        ),
        levels = c(assoc_label, attr_label)
      )
    )

  # -- 8. auto titles ----------------------------------------------------------
  syndrome_tag  <- if (!is.null(syndrome_name)) sprintf(" \u2014 %s", syndrome_name) else ""
  auto_title    <- title    %||% sprintf("%s Associated vs Attributable per Hospital%s", metric, syndrome_tag)
  auto_subtitle <- subtitle %||% "Per 1 000 admissions  |  Attributable = burden directly due to resistance"

  # -- 9. colour vector keyed to legend labels ---------------------------------
  fill_values <- stats::setNames(
    c(colours["associated"], colours["attributable"]),
    c(assoc_label, attr_label)
  )

  # -- 10. build plot ----------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = center, y = value, fill = metric_type)
  ) +
    ggplot2::geom_col(
      position  = ggplot2::position_dodge(width = 0.75),
      color     = "grey20",
      linewidth = 0.3,
      width     = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(value, label_digits)),
      position = ggplot2::position_dodge(width = 0.75),
      hjust    = -0.2,
      size     = 3.8,
      fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill_values, name = "") +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = auto_subtitle,
      x        = "",
      y        = y_label
    ) +
    eda_theme(base_size = base_size, legend_position = "bottom") +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank()
    )

  return(p)
}


# YLL / YLD ASSOCIATED vs ATTRIBUTABLE \u2014 GROUPED BAR BY ORGANISM


#' Plot YLL or YLD Associated vs Attributable by Organism
#'
#' Produces a horizontal grouped bar chart showing the top \code{n} organisms
#' ranked by their associated burden (YLL or YLD), with a paired attributable
#' bar alongside each. Organisms are ranked by the associated value so the
#' highest-burden pathogen appears at the top.
#'
#' The input \code{data} should be a pathogen-level summary data frame with
#' one row per organism, already containing both the associated and attributable
#' burden columns. For YLD this is typically the output of
#' \code{compute_yld_associated()} left-joined with
#' \code{compute_yld_attributable()}. For YLL it is the
#' \code{by_pathogen} slot of \code{compute_yll_associated()} left-joined with
#' the equivalent slot from \code{compute_yll_attributable()}.
#'
#' @param data         Data frame. One row per organism. Must contain the
#'   columns named by \code{organism_col}, \code{assoc_col}, and
#'   \code{attr_col}.
#' @param metric       Character. \code{"YLL"} or \code{"YLD"}. Controls
#'   default column names, colours, and axis labels. Ignored when
#'   \code{assoc_col} and \code{attr_col} are supplied explicitly.
#'   Default \code{"YLL"}.
#' @param n            Integer. Number of top organisms to display. Default 10.
#' @param organism_col Character or \code{NULL}. Column containing organism /
#'   pathogen names. If \code{NULL} (default), auto-set to
#'   \code{"organism_name"} for YLL (matching \code{daly_calc_yll_*} examples)
#'   and \code{"pathogen"} for YLD (matching \code{daly_calc_yld_*} defaults).
#' @param assoc_col    Character or \code{NULL}. Column for the associated
#'   burden value. If \code{NULL} (default), auto-set to
#'   \code{"YLL_associated_k"} (YLL) or \code{"YLD_associated"} (YLD) --
#'   matching the output column names of \code{compute_yll_associated()} and
#'   \code{compute_yld_associated()} respectively.
#' @param attr_col     Character or \code{NULL}. Column for the attributable
#'   burden value. If \code{NULL} (default), auto-set to
#'   \code{"YLL_attributable_k"} (YLL) or \code{"YLD_attributable"} (YLD) --
#'   matching the output column names of \code{compute_yll_attributable()} and
#'   \code{compute_yld_attributable()} respectively.
#' @param n_admissions Numeric or \code{NULL}. When supplied, the associated
#'   and attributable values are divided by \code{n_admissions} and multiplied
#'   by 1 000 before plotting. Default \code{NULL} (no normalisation; values
#'   are plotted in thousands as returned by the YLL/YLD functions).
#' @param syndrome_col  Character or \code{NULL}. Column name containing
#'   syndrome labels. If \code{NULL} (default), no filtering is applied.
#' @param syndrome_name Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}. Default \code{NULL}.
#' @param label_digits Integer. Decimal places in bar value labels. Default 1.
#' @param colours      Named character vector with elements \code{"associated"}
#'   and \code{"attributable"}. If \code{NULL} (default), auto-chosen based
#'   on \code{metric}.
#' @param base_size    Numeric. Base font size. Default 14.
#' @param title        Character. Custom plot title. Auto-generated if
#'   \code{NULL}.
#' @param subtitle     Character. Custom subtitle. Auto-generated if
#'   \code{NULL}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
plot_burden_by_organism <- function(data,
                                    metric        = c("YLL", "YLD"),
                                    n             = 10,
                                    organism_col  = NULL,
                                    assoc_col     = NULL,
                                    attr_col      = NULL,
                                    n_admissions  = NULL,
                                    syndrome_col  = NULL,
                                    syndrome_name = NULL,
                                    label_digits  = 1,
                                    colours       = NULL,
                                    base_size     = 14,
                                    title         = NULL,
                                    subtitle      = NULL) {

  # -- 0. match args -----------------------------------------------------------
  metric <- match.arg(metric)

  # -- 1. resolve column names -------------------------------------------------
  # YLL: pathogen column is "organism_name" (daly_calc_yll_* examples)
  # YLD: pathogen column is "pathogen"      (daly_calc_yld_* defaults)
  if (is.null(organism_col))
    organism_col <- if (metric == "YLL") "organism_name" else "pathogen"
  if (is.null(assoc_col))
    assoc_col <- if (metric == "YLL") "YLL_associated_k" else "YLD_associated"
  if (is.null(attr_col))
    attr_col  <- if (metric == "YLL") "YLL_attributable_k" else "YLD_attributable"

  # -- 2. validate columns -----------------------------------------------------
  required_cols <- c(organism_col, assoc_col, attr_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s",
                 paste(missing_cols, collapse = ", ")))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- 3. syndrome filter ------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 4. normalise if n_admissions provided -----------------------------------
  if (!is.null(n_admissions)) {
    if (!is.numeric(n_admissions) || n_admissions <= 0)
      stop("'n_admissions' must be a positive number.")
    data[[assoc_col]] <- data[[assoc_col]] / n_admissions * 1000
    data[[attr_col]]  <- data[[attr_col]]  / n_admissions * 1000
  }

  # -- 5. default colours ------------------------------------------------------
  if (is.null(colours)) {
    colours <- if (metric == "YLL") {
      c("associated" = "#d73027", "attributable" = "#f4a582")
    } else {
      c("associated" = "#4393c3", "attributable" = "#92c5de")
    }
  }
  if (!all(c("associated", "attributable") %in% names(colours)))
    stop("'colours' must be a named vector with elements 'associated' and 'attributable'.")

  # -- 6. legend label strings -------------------------------------------------
  assoc_label <- paste(metric, "Associated")
  attr_label  <- paste(metric, "Attributable")
  y_label     <- if (!is.null(n_admissions)) {
    sprintf("%s per 1 000 admissions", metric)
  } else {
    sprintf("%s (thousands)", metric)
  }

  # -- 7. select top n organisms by associated value ---------------------------
  org_sym   <- rlang::sym(organism_col)
  assoc_sym <- rlang::sym(assoc_col)
  attr_sym  <- rlang::sym(attr_col)

  top_orgs <- data %>%
    dplyr::filter(!is.na(!!assoc_sym)) %>%
    dplyr::slice_max(!!assoc_sym, n = n, with_ties = FALSE) %>%
    dplyr::arrange(!!assoc_sym) %>%
    dplyr::pull(!!org_sym)

  plot_df <- data %>%
    dplyr::filter(!!org_sym %in% top_orgs) %>%
    dplyr::select(
      organism     = !!org_sym,
      associated   = !!assoc_sym,
      attributable = !!attr_sym
    ) %>%
    tidyr::pivot_longer(
      cols      = c(associated, attributable),
      names_to  = "metric_type",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      organism    = factor(organism, levels = top_orgs),
      metric_type = factor(
        dplyr::case_when(
          metric_type == "associated"   ~ assoc_label,
          metric_type == "attributable" ~ attr_label
        ),
        levels = c(assoc_label, attr_label)
      )
    )

  # -- 7. auto titles ----------------------------------------------------------
  syndrome_tag  <- if (!is.null(syndrome_name)) sprintf(" \u2014 %s", syndrome_name) else ""
  auto_title    <- title    %||% sprintf("Top %d Organisms by %s Burden%s", n, metric, syndrome_tag)
  auto_subtitle <- subtitle %||% if (!is.null(n_admissions)) {
    "Per 1 000 admissions  |  Attributable = burden directly due to resistance"
  } else {
    "Attributable = burden directly due to resistance"
  }

  # -- 8. colour vector keyed to legend labels ---------------------------------
  fill_values <- stats::setNames(
    c(colours["associated"], colours["attributable"]),
    c(assoc_label, attr_label)
  )

  # -- 9. build plot -----------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = organism, y = value, fill = metric_type)
  ) +
    ggplot2::geom_col(
      position  = ggplot2::position_dodge(width = 0.75),
      color     = "grey20",
      linewidth = 0.3,
      width     = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(value, label_digits)),
      position = ggplot2::position_dodge(width = 0.75),
      hjust    = -0.2,
      size     = 3.8,
      fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill_values, name = "") +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = auto_subtitle,
      x        = "",
      y        = y_label
    ) +
    eda_theme(base_size = base_size, legend_position = "bottom") +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank()
    )

  return(p)
}


# YLL HEATMAP \u2014 RESISTANCE CLASS x PATHOGEN GROUP

#' Heatmap of YLL per 1 000 Admissions by Resistance Class and Pathogen Group
#'
#' Produces a tile heatmap where the x-axis shows resistance classes (dominant
#' AMR profiles), the y-axis shows pathogen groups (or individual pathogens),
#' and each cell's fill reflects the associated or attributable YLL burden per
#' 1 000 admissions.
#'
#' The input \code{data} should be the \code{$by_class} slot (for associated)
#' or the \code{$by_profile} slot (for attributable) returned by
#' \code{daly_calc_yll_associated()} / \code{daly_calc_yll_attributable()}.
#' Both slots contain columns \code{pathogen}, \code{profile}, and
#' \code{YLL_Kdelta}.
#'
#' If an organism-group column has been pre-joined onto \code{data} (e.g.
#' \code{org_group}), supply its name via \code{group_col} and the function
#' will aggregate to the group x class level. When \code{group_col} is
#' \code{NULL} or absent, individual pathogen names are used on the y-axis.
#'
#' @param data         Data frame. Must contain the columns named by
#'   \code{pathogen_col}, \code{class_col}, and \code{value_col}.
#' @param type         Character. \code{"associated"} (default) or
#'   \code{"attributable"}. Used only for axis and title labels.
#' @param n_admissions Numeric or \code{NULL}. Total admissions used as the
#'   denominator for per-1 000 normalisation. If \code{NULL} (default),
#'   values are plotted as supplied (assumed pre-normalised).
#' @param pathogen_col Character. Column containing pathogen names.
#'   Default \code{"pathogen"}.
#' @param class_col    Character. Column containing the resistance class /
#'   dominant profile. Default \code{"profile"}.
#' @param value_col    Character. Required. Column containing the YLL value to
#'   plot (e.g. \code{"YLL_class"} after distributing per-pathogen YLL across
#'   resistance classes, or \code{"YLL_attributable_k"} if your data already
#'   has one row per pathogen x class).
#' @param group_col    Character or \code{NULL}. Column containing organism
#'   group labels (e.g. \code{"org_group"}). If \code{NULL} (default) or not
#'   present in \code{data}, pathogen names are used on the y-axis.
#' @param syndrome_col  Character or \code{NULL}. Syndrome filter column.
#'   Default \code{NULL}.
#' @param syndrome_name Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}. Default \code{NULL}.
#' @param show_values  Logical. Print rounded values inside each cell.
#'   Default \code{TRUE}.
#' @param value_digits Integer. Decimal places for cell labels. Default 2.
#' @param palette      Character. RColorBrewer palette for the fill scale.
#'   Default \code{"Spectral"}.
#' @param base_size    Numeric. Base font size. Default 14.
#' @param title        Character. Custom plot title. Auto-generated if
#'   \code{NULL}.
#' @param subtitle     Character. Custom subtitle. Auto-generated if
#'   \code{NULL}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_yll_heatmap <- function(data,
                             type          = c("associated", "attributable"),
                             value_col,
                             n_admissions  = NULL,
                             pathogen_col  = "pathogen",
                             class_col     = "profile",
                             group_col     = NULL,
                             syndrome_col  = NULL,
                             syndrome_name = NULL,
                             show_values   = TRUE,
                             value_digits  = 2,
                             palette       = "Spectral",
                             base_size     = 14,
                             title         = NULL,
                             subtitle      = NULL) {

  # -- 0. match args -----------------------------------------------------------
  type <- match.arg(type)

  # -- 1. validate required args -----------------------------------------------
  if (missing(value_col) || is.null(value_col))
    stop("'value_col' is required. Supply the column name containing your YLL values (e.g. \"YLL_class\" or \"YLL_attributable_k\").")
  if (!is.null(n_admissions) && (!is.numeric(n_admissions) || n_admissions <= 0))
    stop("'n_admissions' must be a positive number.")

  # -- 2. determine y-axis column (group if available, else pathogen) -----------
  use_group <- !is.null(group_col) && group_col %in% names(data)
  y_col     <- if (use_group) group_col else pathogen_col

  # -- 3. validate columns -----------------------------------------------------
  required_cols <- unique(c(pathogen_col, class_col, value_col,
                             if (use_group) group_col else NULL))
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s", paste(missing_cols, collapse = ", ")))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- 4. syndrome filter ------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) & data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 5. aggregate and normalise ----------------------------------------------
  y_sym     <- rlang::sym(y_col)
  class_sym <- rlang::sym(class_col)
  value_sym <- rlang::sym(value_col)

  heat_df <- data %>%
    dplyr::group_by(!!y_sym, !!class_sym) %>%
    dplyr::summarise(yll_sum = sum(!!value_sym, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      yll_per_1000 = if (!is.null(n_admissions)) yll_sum / n_admissions * 1000 else yll_sum,
      yll_per_1000 = ifelse(yll_per_1000 < 0, NA_real_, yll_per_1000)
    )

  # -- 6. axis ordering by total burden ----------------------------------------
  class_order <- heat_df %>%
    dplyr::group_by(!!class_sym) %>%
    dplyr::summarise(tot = sum(yll_per_1000, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(tot)) %>%
    dplyr::pull(!!class_sym)

  group_order <- heat_df %>%
    dplyr::group_by(!!y_sym) %>%
    dplyr::summarise(tot = sum(yll_per_1000, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(tot) %>%
    dplyr::pull(!!y_sym)

  heat_df <- heat_df %>%
    dplyr::mutate(
      !!class_sym := factor(!!class_sym, levels = class_order),
      !!y_sym      := factor(!!y_sym,    levels = group_order)
    )

  # -- 7. midpoint for dynamic text colour -------------------------------------
  mid_val <- stats::median(heat_df$yll_per_1000, na.rm = TRUE)

  # -- 8. auto titles ----------------------------------------------------------
  type_label    <- if (type == "associated") "Associated" else "Attributable"
  syndrome_tag  <- if (!is.null(syndrome_name)) sprintf(" \u2014 %s", syndrome_name) else ""
  auto_title    <- title    %||% sprintf("YLL %s per 1 000 Admissions by Resistance Class%s",
                                         type_label, syndrome_tag)
  auto_subtitle <- subtitle %||% "YLL delta (thousands) normalised to total admissions; negative values masked"

  # -- 9. build base plot ------------------------------------------------------
  p <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = !!class_sym, y = !!y_sym, fill = yll_per_1000)
  ) +
    ggplot2::geom_tile(colour = "grey80", linewidth = 0.4) +
    ggplot2::scale_fill_distiller(
      palette   = palette,
      direction = -1,
      na.value  = "white",
      name      = "YLL\nper 1 000",
      guide     = ggplot2::guide_colourbar(barwidth = 1, barheight = 8)
    ) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = auto_subtitle,
      x        = "Resistance class",
      y        = ""
    ) +
    eda_theme(base_size = base_size, legend_position = "right") +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 40, hjust = 1),
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "grey70", fill = NA)
    )

  # -- 10. optional cell labels with dynamic text colour -----------------------
  if (show_values) {
    dark_df  <- heat_df %>%
      dplyr::filter(!is.na(yll_per_1000) & yll_per_1000 >= mid_val) %>%
      dplyr::mutate(cell_lbl = format(round(yll_per_1000, value_digits), nsmall = value_digits))
    light_df <- heat_df %>%
      dplyr::filter(!is.na(yll_per_1000) & yll_per_1000 < mid_val) %>%
      dplyr::mutate(cell_lbl = format(round(yll_per_1000, value_digits), nsmall = value_digits))

    if (nrow(dark_df) > 0)
      p <- p + ggplot2::geom_text(
        data = dark_df,
        ggplot2::aes(label = cell_lbl),
        colour = "white", size = 3.2, fontface = "bold"
      )
    if (nrow(light_df) > 0)
      p <- p + ggplot2::geom_text(
        data = light_df,
        ggplot2::aes(label = cell_lbl),
        colour = "grey20", size = 3.2
      )
  }

  return(p)
}


# YLD HEATMAP \u2014 ORGANISM GROUP x ASSOCIATED / ATTRIBUTABLE

#' Heatmap of YLD per 1 000 Admissions by Organism Group
#'
#' Produces a two-column tile heatmap where the x-axis shows YLD Associated
#' and YLD Attributable, the y-axis shows organism groups (or individual
#' pathogens), and each cell's fill reflects the YLD burden per 1 000
#' admissions.
#'
#' The input \code{data} should contain one row per pathogen with
#' \code{YLD_associated} and \code{YLD_attributable} columns -- typically
#' the output of \code{daly_calc_yld_attributable()} (which augments the
#' data frame with both columns). If an organism-group column has been
#' pre-joined (e.g. \code{org_group}), supply its name via \code{group_col}
#' to aggregate cells to the group level.
#'
#' @param data         Data frame. Must contain the columns named by
#'   \code{pathogen_col}, \code{assoc_col}, and \code{attr_col}.
#' @param n_admissions Numeric or \code{NULL}. Total admissions used as the
#'   denominator for per-1 000 normalisation. If \code{NULL} (default),
#'   values are plotted as supplied (assumed pre-normalised).
#' @param pathogen_col Character. Column containing pathogen names.
#'   Default \code{"pathogen"}.
#' @param assoc_col    Character. Column containing YLD associated values.
#'   Default \code{"YLD_associated"}.
#' @param attr_col     Character. Column containing YLD attributable values.
#'   Default \code{"YLD_attributable"}.
#' @param group_col    Character or \code{NULL}. Column containing organism
#'   group labels (e.g. \code{"org_group"}). If \code{NULL} (default) or not
#'   present in \code{data}, pathogen names are used on the y-axis.
#' @param syndrome_col  Character or \code{NULL}. Syndrome filter column.
#'   Default \code{NULL}.
#' @param syndrome_name Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}. Default \code{NULL}.
#' @param show_values  Logical. Print rounded values inside each cell.
#'   Default \code{TRUE}.
#' @param value_digits Integer. Decimal places for cell labels. Default 2.
#' @param palette      Character. RColorBrewer palette for the fill scale.
#'   Default \code{"Blues"}.
#' @param coord_fixed  Logical. Whether to apply \code{ggplot2::coord_fixed()}
#'   to produce square tiles. Default \code{TRUE}.
#' @param base_size    Numeric. Base font size. Default 14.
#' @param title        Character. Custom plot title. Auto-generated if
#'   \code{NULL}.
#' @param subtitle     Character. Custom subtitle. Auto-generated if
#'   \code{NULL}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_yld_heatmap <- function(data,
                             n_admissions  = NULL,
                             pathogen_col  = "pathogen",
                             assoc_col     = "YLD_associated",
                             attr_col      = "YLD_attributable",
                             group_col     = NULL,
                             syndrome_col  = NULL,
                             syndrome_name = NULL,
                             show_values   = TRUE,
                             value_digits  = 2,
                             palette       = "Blues",
                             coord_fixed   = TRUE,
                             base_size     = 14,
                             title         = NULL,
                             subtitle      = NULL) {

  # -- 1. validate n_admissions if provided ------------------------------------
  if (!is.null(n_admissions) && (!is.numeric(n_admissions) || n_admissions <= 0))
    stop("'n_admissions' must be a positive number.")

  # -- 2. determine y-axis column ----------------------------------------------
  use_group <- !is.null(group_col) && group_col %in% names(data)
  y_col     <- if (use_group) group_col else pathogen_col

  # -- 3. validate columns -----------------------------------------------------
  required_cols <- unique(c(pathogen_col, assoc_col, attr_col,
                             if (use_group) group_col else NULL))
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s", paste(missing_cols, collapse = ", ")))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- 4. syndrome filter ------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) & data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 5. aggregate by y-axis group --------------------------------------------
  y_sym     <- rlang::sym(y_col)
  assoc_sym <- rlang::sym(assoc_col)
  attr_sym  <- rlang::sym(attr_col)

  agg_df <- data %>%
    dplyr::group_by(!!y_sym) %>%
    dplyr::summarise(
      associated   = sum(!!assoc_sym, na.rm = TRUE),
      attributable = sum(!!attr_sym,  na.rm = TRUE),
      .groups = "drop"
    )

  # -- 6. normalise (if n_admissions provided) and pivot to long ---------------
  heat_df <- agg_df %>%
    dplyr::mutate(
      associated   = if (!is.null(n_admissions)) associated   / n_admissions * 1000 else associated,
      attributable = if (!is.null(n_admissions)) attributable / n_admissions * 1000 else attributable
    ) %>%
    tidyr::pivot_longer(
      cols      = c(associated, attributable),
      names_to  = "metric_type",
      values_to = "yld_per_1000"
    ) %>%
    dplyr::mutate(
      yld_per_1000 = ifelse(yld_per_1000 < 0, NA_real_, yld_per_1000),
      metric_type  = dplyr::recode(metric_type,
        associated   = "YLD Associated",
        attributable = "YLD Attributable"
      ),
      metric_type = factor(metric_type, levels = c("YLD Associated", "YLD Attributable"))
    )

  # -- 7. y-axis ordering by total YLD (highest at top) -----------------------
  group_order <- heat_df %>%
    dplyr::group_by(!!y_sym) %>%
    dplyr::summarise(tot = sum(yld_per_1000, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(tot) %>%
    dplyr::pull(!!y_sym)

  heat_df <- heat_df %>%
    dplyr::mutate(!!y_sym := factor(!!y_sym, levels = group_order))

  # -- 8. midpoint for dynamic text colour ------------------------------------
  mid_val <- stats::median(heat_df$yld_per_1000, na.rm = TRUE)

  # -- 9. auto titles ----------------------------------------------------------
  syndrome_tag  <- if (!is.null(syndrome_name)) sprintf(" \u2014 %s", syndrome_name) else ""
  auto_title    <- title    %||% sprintf("YLD per 1 000 Admissions by Organism Group%s", syndrome_tag)
  auto_subtitle <- subtitle %||% "Attributable = YLD directly due to resistance"

  # -- 10. build base plot -----------------------------------------------------
  p <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = metric_type, y = !!y_sym, fill = yld_per_1000)
  ) +
    ggplot2::geom_tile(colour = "grey80", linewidth = 0.4) +
    ggplot2::scale_fill_distiller(
      palette   = palette,
      direction = 1,
      na.value  = "grey92",
      name      = "YLD\nper 1 000",
      guide     = ggplot2::guide_colourbar(barwidth = 1, barheight = 8)
    ) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = auto_subtitle,
      x        = "",
      y        = "Pathogen group"
    ) +
    eda_theme(base_size = base_size, legend_position = "right") +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "grey80", fill = NA, linewidth = 0.5)
    )

  if (coord_fixed)
    p <- p + ggplot2::coord_fixed()

  # -- 11. optional cell labels ------------------------------------------------
  if (show_values) {
    label_df <- heat_df %>%
      dplyr::filter(!is.na(yld_per_1000)) %>%
      dplyr::mutate(cell_lbl = format(round(yld_per_1000, value_digits), nsmall = value_digits))

    p <- p + ggplot2::geom_text(
      data     = label_df,
      ggplot2::aes(label = cell_lbl),
      colour   = "black",
      size     = 3.5,
      fontface = "bold"
    )
  }

  return(p)
}
