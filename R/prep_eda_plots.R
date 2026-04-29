# eda_plots.R
# Reusable EDA plotting functions for AMR stewardship datasets


# THEME & PALETTE HELPERS


#' EDA ggplot2 Theme
#'
#' Consistent theme used across all EDA plot functions in this file.
#' Built on top of \code{ggpubr::theme_pubr()}.
#'
#' @param base_size Numeric. Base font size. Default 14.
#' @param legend_position Character. Legend position ("top", "bottom",
#'   "left", "right", "none"). Default "top".
#'
#' @return A ggplot2 theme object.
#'
eda_theme <- function(base_size = 14, legend_position = "top") {
  ggpubr::theme_pubr(base_size = base_size) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "black",
                                               linewidth = 0.8),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA,
                                               linewidth = 0.8),
      panel.spacing    = ggplot2::unit(1.2, "lines"),
      legend.position  = legend_position,
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}



# TOP PATHOGENS


#' Plot Top Pathogens by Unique Patients
#'
#' Produces a horizontal bar chart of the top \code{n} organisms ranked by
#' unique patient count. Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} -- one panel per centre, top \emph{n} organisms
#'     within each centre independently. Facet strips show patient and organism
#'     counts for that centre.
#'   \item \strong{"overall"} -- all centres pooled; single chart of the top
#'     \emph{n} organisms across the entire dataset.
#'   \item \strong{"single"} -- one specific centre; pass the centre name via
#'     \code{center}.
#' }
#'
#' Labels on each bar show absolute count and percentage of that centre's
#' (or overall) unique patient total.
#'
#' @param data         Data frame. Long-format AMR dataset (one row per
#'   isolate/antibiotic record).
#' @param n            Integer. Number of top organisms to show. Default 5.
#' @param mode         Character. One of \code{"faceted"}, \code{"overall"},
#'   or \code{"single"}. Default \code{"faceted"}.
#' @param center       Character. Required when \code{mode = "single"}.
#'   Exact name of the centre to plot (must match a value in
#'   \code{center_col}).
#' @param patient_col  Character. Column name for patient identifiers.
#'   Default \code{"PatientInformation_id"}.
#' @param organism_col Character. Column name for organism names.
#'   Default \code{"organism_name"}.
#' @param center_col   Character. Column name for centre/facility names.
#'   Default \code{"center_name"}.
#' @param fill_colour  Character. Bar fill colour (any R colour string or hex).
#'   Default \code{"steelblue"}.
#' @param ncol         Integer. Number of columns in the facet grid
#'   (\code{mode = "faceted"} only). Default 2.
#' @param base_size    Numeric. Base font size passed to \code{eda_theme()}.
#'   Default 14.
#' @param title        Character. Custom plot title. If \code{NULL} a title is
#'   generated automatically from \code{mode} and \code{n}. Default \code{NULL}.
#' @param syndrome_col Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name Character or \code{NULL}. The syndrome value to retain
#'   (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col} to be set.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object. Print it or pass to \code{ggsave()}.
#' @export
#'

plot_top_organisms <- function(data,
                               n             = 5,
                               mode          = c("faceted", "overall", "single"),
                               center        = NULL,
                               patient_col   = "PatientInformation_id",
                               organism_col  = "organism_name",
                               center_col    = "center_name",
                               fill_colour   = "steelblue",
                               ncol          = 2,
                               base_size     = 14,
                               title         = NULL,
                               syndrome_col  = NULL,
                               syndrome_name = NULL) {

  # -- 0. match mode argument ------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate required columns -----------------------------------------
  required_cols <- c(patient_col, organism_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center)) {
      stop("'center' must be provided when mode = 'single'.")
    }
    available <- unique(data[[center_col]])
    if (!center %in% available) {
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col,
        paste(sort(available), collapse = ", ")
      ))
    }
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  org_sym <- rlang::sym(organism_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. cleaning -- remove blank / NA organism names -----------------------
  org_data <- data %>%
    dplyr::filter(
      !is.na(!!org_sym),
      trimws(as.character(!!org_sym)) != ""
    )

  # -- 4. optional: filter to a single centre -------------------------------
  if (mode == "single") {
    org_data <- org_data %>%
      dplyr::filter(!!ctr_sym == center)
  }

  # -- 5. branch by mode -----------------------------------------------------

  # -- 5a. OVERALL -----------------------------------------------------------
  if (mode == "overall") {

    total_patients <- dplyr::n_distinct(org_data[[patient_col]])

    top_org <- org_data %>%
      dplyr::group_by(!!org_sym) %>%
      dplyr::summarise(
        patients = dplyr::n_distinct(!!pt_sym),
        .groups  = "drop"
      ) %>%
      dplyr::slice_max(patients, n = n) %>%
      dplyr::mutate(
        percent = 100 * patients / total_patients,
        label   = paste0(patients, " (", round(percent, 1), "%)"),
        !!organism_col := stats::reorder(!!org_sym, patients)
      )

    auto_title <- title %||% sprintf("Top %d Organisms by Unique Patients \u2014 All Centres Pooled", n)

    p <- ggplot2::ggplot(
      top_org,
      ggplot2::aes(
        x = !!org_sym,
        y = patients
      )
    ) +
      ggplot2::geom_col(fill = fill_colour, width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        hjust     = -0.1,
        size      = 3.8,
        fontface  = "bold"
      ) +
      ggplot2::coord_flip() +
      ggplot2::expand_limits(y = max(top_org$patients) * 1.25) +
      ggplot2::labs(
        x     = "Organism",
        y     = "Number of Unique Patients",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none")

    return(p)
  }

  # -- 5b. FACETED or SINGLE (share the same per-centre logic) --------------

  # Centre summary for facet strip labels
  center_summary <- org_data %>%
    dplyr::group_by(!!ctr_sym) %>%
    dplyr::summarise(
      unique_patients  = dplyr::n_distinct(!!pt_sym),
      unique_organisms = dplyr::n_distinct(!!org_sym),
      .groups          = "drop"
    ) %>%
    dplyr::mutate(
      facet_label = paste0(
        !!ctr_sym,
        "\nPatients: ",  unique_patients,
        " | Organisms: ", unique_organisms
      )
    )

  # Total unique patients per centre (denominator for %)
  total_per_center <- org_data %>%
    dplyr::group_by(!!ctr_sym) %>%
    dplyr::summarise(
      total_patients = dplyr::n_distinct(!!pt_sym),
      .groups        = "drop"
    )

  # Top n organisms per centre
  top_org <- org_data %>%
    dplyr::group_by(!!ctr_sym, !!org_sym) %>%
    dplyr::summarise(
      patients = dplyr::n_distinct(!!pt_sym),
      .groups  = "drop"
    ) %>%
    dplyr::left_join(total_per_center, by = center_col) %>%
    dplyr::mutate(
      percent = 100 * patients / total_patients,
      label   = paste0(patients, " (", round(percent, 1), "%)")
    ) %>%
    dplyr::group_by(!!ctr_sym) %>%
    dplyr::slice_max(patients, n = n) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      center_summary[, c(center_col, "facet_label")],
      by = center_col
    )

  # -- 5c. SINGLE: drop facet_wrap, use just the centre name in the title ----
  if (mode == "single") {

    auto_title <- title %||% sprintf(
      "Top %d Organisms by Unique Patients \u2014 %s", n, center
    )

    top_org <- top_org %>%
      dplyr::mutate(
        !!organism_col := stats::reorder(!!org_sym, patients)
      )

    p <- ggplot2::ggplot(
      top_org,
      ggplot2::aes(
        x = !!org_sym,
        y = patients
      )
    ) +
      ggplot2::geom_col(fill = fill_colour, width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        hjust    = -0.1,
        size     = 3.8,
        fontface = "bold"
      ) +
      ggplot2::coord_flip() +
      ggplot2::expand_limits(y = max(top_org$patients) * 1.25) +
      ggplot2::labs(
        x     = "Organism",
        y     = "Number of Unique Patients",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none")

    return(p)
  }

  # -- 5d. FACETED -----------------------------------------------------------

  auto_title <- title %||% sprintf(
    "Top %d Organisms by Unique Patients \u2014 All Centres", n
  )

  p <- ggplot2::ggplot(
    top_org,
    ggplot2::aes(
      x = tidytext::reorder_within(!!org_sym, patients, !!ctr_sym),
      y = patients
    )
  ) +
    ggplot2::geom_col(fill = fill_colour, width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      hjust    = -0.1,
      size     = 3.8,
      fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    tidytext::scale_x_reordered() +
    ggplot2::facet_wrap(~facet_label, ncol = ncol, scales = "free_y") +
    ggplot2::expand_limits(y = max(top_org$patients) * 1.25) +
    ggplot2::labs(
      x     = "Organism",
      y     = "Number of Unique Patients",
      title = auto_title
    ) +
    eda_theme(base_size = base_size, legend_position = "none")

  return(p)
}



# ANTIBIOTIC SUSCEPTIBILITY PATTERN


#' Plot Antibiotic Susceptibility Pattern (Stacked R/S Bars)
#'
#' Produces a horizontal stacked bar chart showing the number of admissions
#' tested as Resistant (R) or Susceptible (S) for the top \code{n} antibiotics.
#' Percentage labels are shown inside each bar segment.
#'
#' \strong{Worst-phenotype rule:} when a patient has multiple records for the
#' same organism-antibiotic combination, any single R result marks the whole
#' episode as R. This prevents double-counting repeat cultures from the same
#' infection episode.
#'
#' \strong{Counting unit:} bars show unique patients (not isolate rows), so a
#' patient tested against the same antibiotic twice is counted once.
#'
#' Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} -- top \emph{n} antibiotics per centre
#'     independently, one panel per centre.
#'   \item \strong{"overall"} -- all centres pooled; top \emph{n} antibiotics
#'     globally.
#'   \item \strong{"single"} -- one specific centre; pass the centre name via
#'     \code{center}.
#' }
#'
#' @param data           Data frame. Long-format AMR dataset.
#' @param n              Integer. Number of top antibiotics to show (ranked by
#'   total patients tested). Default 5.
#' @param mode           Character. One of \code{"faceted"}, \code{"overall"},
#'   or \code{"single"}. Default \code{"faceted"}.
#' @param center         Character. Required when \code{mode = "single"}.
#'   Exact name of the centre to plot.
#' @param patient_col    Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param antibiotic_col Character. Antibiotic name column.
#'   Default \code{"antibiotic_name"}.
#' @param value_col      Character. Susceptibility result column containing
#'   "R" and "S" values. Default \code{"antibiotic_value"}.
#' @param organism_col   Character. Organism name column (used in
#'   worst-phenotype deduplication). Default \code{"organism_name"}.
#' @param center_col     Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param colours        Named character vector. Fill colours for R and S.
#'   Default \code{c("R" = "#D73027", "S" = "#1A9850")}.
#' @param ncol           Integer. Facet columns (\code{mode = "faceted"} only).
#'   Default 2.
#' @param base_size      Numeric. Base font size. Default 14.
#' @param title          Character. Custom plot title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col   Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name  Character or \code{NULL}. The syndrome value to retain
#'   (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col} to be set.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
#'

plot_abx_susceptibility <- function(data,
                                    n              = 5,
                                    mode           = c("faceted", "overall", "single"),
                                    center         = NULL,
                                    patient_col    = "PatientInformation_id",
                                    antibiotic_col = "antibiotic_name",
                                    value_col      = "antibiotic_value",
                                    organism_col   = "organism_name",
                                    center_col     = "center_name",
                                    colours        = c("R" = "#D73027", "S" = "#1A9850"),
                                    ncol           = 2,
                                    base_size      = 14,
                                    title          = NULL,
                                    syndrome_col   = NULL,
                                    syndrome_name  = NULL) {

  # -- 0. match mode ---------------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, antibiotic_col, value_col, organism_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center)) {
      stop("'center' must be provided when mode = 'single'.")
    }
    available <- unique(data[[center_col]])
    if (!center %in% available) {
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col,
        paste(sort(available), collapse = ", ")
      ))
    }
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  abx_sym <- rlang::sym(antibiotic_col)
  val_sym <- rlang::sym(value_col)
  org_sym <- rlang::sym(organism_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. clean: keep only valid R/S rows, drop blank names -----------------
  abx_clean <- data %>%
    dplyr::filter(
      !!val_sym %in% c("R", "S"),
      !is.na(!!abx_sym), trimws(as.character(!!abx_sym)) != "",
      !is.na(!!org_sym),  trimws(as.character(!!org_sym))  != ""
    ) %>%
    dplyr::distinct(
      !!ctr_sym, !!pt_sym, !!org_sym, !!abx_sym, !!val_sym
    )

  # -- 4. optional: filter to a single centre -------------------------------
  if (mode == "single") {
    abx_clean <- abx_clean %>%
      dplyr::filter(!!ctr_sym == center)
  }

  # -- 5. worst-phenotype rule -----------------------------------------------
  # Per patient x organism x antibiotic: any R -> episode = R
  if (mode == "overall") {
    abx_episode <- abx_clean %>%
      dplyr::group_by(!!pt_sym, !!org_sym, !!abx_sym) %>%
      dplyr::summarise(
        final_abx = ifelse(any(!!val_sym == "R"), "R", "S"),
        .groups   = "drop"
      )
  } else {
    abx_episode <- abx_clean %>%
      dplyr::group_by(!!ctr_sym, !!pt_sym, !!org_sym, !!abx_sym) %>%
      dplyr::summarise(
        final_abx = ifelse(any(!!val_sym == "R"), "R", "S"),
        .groups   = "drop"
      )
  }

  # -- 6. summarise: count unique patients per antibiotic x R/S -------------
  if (mode == "overall") {

    abx_summary <- abx_episode %>%
      dplyr::group_by(!!abx_sym, final_abx) %>%
      dplyr::summarise(
        n       = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!abx_sym) %>%
      dplyr::mutate(
        total_tested = sum(n),
        percent      = round(100 * n / total_tested, 1),
        label        = paste0(percent, "%")
      ) %>%
      dplyr::ungroup()

    # Top n antibiotics
    top_abx <- abx_summary %>%
      dplyr::distinct(!!abx_sym, total_tested) %>%
      dplyr::slice_max(total_tested, n = n)

    abx_summary <- abx_summary %>%
      dplyr::inner_join(top_abx, by = c(antibiotic_col, "total_tested")) %>%
      dplyr::mutate(
        !!antibiotic_col := stats::reorder(!!abx_sym, total_tested)
      )

    auto_title <- title %||% sprintf(
      "Antibiotic Susceptibility Pattern (Top %d) \u2014 All Centres Pooled", n
    )

    p <- ggplot2::ggplot(
      abx_summary,
      ggplot2::aes(
        x    = !!abx_sym,
        y    = n,
        fill = final_abx
      )
    ) +
      ggplot2::geom_col(width = 0.7, color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        color    = "white",
        size     = 4,
        fontface = "bold"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(
        values = colours,
        labels = c("R" = "Resistant", "S" = "Susceptible")
      ) +
      ggplot2::labs(
        x     = "Antibiotic",
        y     = "Number of Admissions",
        fill  = "Susceptibility",
        title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 7. per-centre summary (faceted + single) ------------------------------
  abx_summary <- abx_episode %>%
    dplyr::group_by(!!ctr_sym, !!abx_sym, final_abx) %>%
    dplyr::summarise(
      n       = dplyr::n_distinct(!!pt_sym),
      .groups = "drop"
    ) %>%
    dplyr::group_by(!!ctr_sym, !!abx_sym) %>%
    dplyr::mutate(
      total_tested = sum(n),
      percent      = round(100 * n / total_tested, 1),
      label        = paste0(percent, "%")
    ) %>%
    dplyr::ungroup()

  # Top n antibiotics per centre
  top_abx <- abx_summary %>%
    dplyr::distinct(!!ctr_sym, !!abx_sym, total_tested) %>%
    dplyr::group_by(!!ctr_sym) %>%
    dplyr::slice_max(total_tested, n = n) %>%
    dplyr::ungroup()

  abx_summary <- abx_summary %>%
    dplyr::inner_join(top_abx, by = c(center_col, antibiotic_col, "total_tested"))

  # -- 8. SINGLE: no facet ---------------------------------------------------
  if (mode == "single") {

    auto_title <- title %||% sprintf(
      "Antibiotic Susceptibility Pattern (Top %d) \u2014 %s", n, center
    )

    abx_summary <- abx_summary %>%
      dplyr::mutate(
        !!antibiotic_col := stats::reorder(!!abx_sym, total_tested)
      )

    p <- ggplot2::ggplot(
      abx_summary,
      ggplot2::aes(
        x    = !!abx_sym,
        y    = n,
        fill = final_abx
      )
    ) +
      ggplot2::geom_col(width = 0.7, color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        color    = "white",
        size     = 4,
        fontface = "bold"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(
        values = colours,
        labels = c("R" = "Resistant", "S" = "Susceptible")
      ) +
      ggplot2::labs(
        x     = "Antibiotic",
        y     = "Number of Admissions",
        fill  = "Susceptibility",
        title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 9. FACETED ------------------------------------------------------------
  auto_title <- title %||% sprintf(
    "Antibiotic Susceptibility Pattern (Top %d) \u2014 All Centres", n
  )

  p <- ggplot2::ggplot(
    abx_summary,
    ggplot2::aes(
      x    = tidytext::reorder_within(!!abx_sym, total_tested, !!ctr_sym),
      y    = n,
      fill = final_abx
    )
  ) +
    ggplot2::geom_col(width = 0.7, color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      color    = "white",
      size     = 4,
      fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    tidytext::scale_x_reordered() +
    ggplot2::facet_wrap(stats::as.formula(paste("~", center_col)), ncol = ncol, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = colours,
      labels = c("R" = "Resistant", "S" = "Susceptible")
    ) +
    ggplot2::labs(
      x     = "Antibiotic",
      y     = "Number of Admissions",
      fill  = "Susceptibility",
      title = auto_title
    ) +
    eda_theme(base_size = base_size)

  return(p)
}



# ANTIBIOTIC RESISTANCE HEATMAP


#' Plot Antibiotic Resistance Heatmap
#'
#' Produces a ggplot2 tile heatmap where:
#' \itemize{
#'   \item \strong{x-axis} -- antibiotics (ordered by total patients tested,
#'     most tested on the left)
#'   \item \strong{y-axis} -- centres
#'   \item \strong{fill colour} -- proportion of patients who were Resistant
#'     for that antibiotic at that centre (0 = fully susceptible, 1 = fully
#'     resistant)
#' }
#'
#' Only R and S results are used. The same \strong{worst-phenotype rule}
#' as \code{plot_abx_susceptibility()} is applied: per patient x organism x
#' antibiotic, any single R marks the episode as R.
#'
#' Cells with no data (antibiotic not tested at that centre) are shown in
#' light grey.
#'
#' Supports two display modes:
#' \itemize{
#'   \item \strong{"all"} (default) -- all centres on the y-axis.
#'   \item \strong{"single"} -- filter to one centre before plotting (produces
#'     a single-row heatmap useful for per-centre reports).
#' }
#'
#' @param data           Data frame. Long-format AMR dataset.
#' @param n              Integer or \code{NULL}. If an integer, only the top
#'   \emph{n} antibiotics by total patients tested are shown. If \code{NULL},
#'   all antibiotics are shown. Default \code{NULL}.
#' @param mode           Character. \code{"all"} or \code{"single"}.
#'   Default \code{"all"}.
#' @param center         Character. Required when \code{mode = "single"}.
#'   Exact centre name.
#' @param patient_col    Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param antibiotic_col Character. Antibiotic name column.
#'   Default \code{"antibiotic_name"}.
#' @param value_col      Character. Susceptibility result column (R/S values).
#'   Default \code{"antibiotic_value"}.
#' @param organism_col   Character. Organism name column.
#'   Default \code{"organism_name"}.
#' @param center_col     Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param show_values    Logical. Print proportion values inside tiles.
#'   Default \code{TRUE}.
#' @param midpoint       Numeric. Midpoint of the colour gradient (0-1).
#'   Default \code{0.5}.
#' @param low_colour     Character. Colour for proportion = 0 (fully
#'   susceptible). Default \code{"#1a9850"} (green).
#' @param mid_colour     Character. Colour for the midpoint.
#'   Default \code{"#fee08b"} (yellow).
#' @param high_colour    Character. Colour for proportion = 1 (fully
#'   resistant). Default \code{"#d73027"} (red).
#' @param base_size      Numeric. Base font size. Default 14.
#' @param title          Character. Custom plot title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col   Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name  Character or \code{NULL}. The syndrome value to retain
#'   (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col} to be set.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
#'

plot_abx_heatmap <- function(data,
                             n              = NULL,
                             mode           = c("all", "single"),
                             center         = NULL,
                             patient_col    = "PatientInformation_id",
                             antibiotic_col = "antibiotic_name",
                             value_col      = "antibiotic_value",
                             organism_col   = "organism_name",
                             center_col     = "center_name",
                             show_values    = TRUE,
                             midpoint       = 0.5,
                             low_colour     = "#1a9850",
                             mid_colour     = "#fee08b",
                             high_colour    = "#d73027",
                             base_size      = 14,
                             title          = NULL,
                             syndrome_col   = NULL,
                             syndrome_name  = NULL) {

  # -- 0. match mode ---------------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, antibiotic_col, value_col,
                     organism_col, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center)) {
      stop("'center' must be provided when mode = 'single'.")
    }
    available <- unique(data[[center_col]])
    if (!center %in% available) {
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col,
        paste(sort(available), collapse = ", ")
      ))
    }
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  abx_sym <- rlang::sym(antibiotic_col)
  val_sym <- rlang::sym(value_col)
  org_sym <- rlang::sym(organism_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. clean: keep only valid R/S rows, drop blank names -----------------
  abx_clean <- data %>%
    dplyr::filter(
      !!val_sym %in% c("R", "S"),
      !is.na(!!abx_sym), trimws(as.character(!!abx_sym)) != "",
      !is.na(!!org_sym),  trimws(as.character(!!org_sym))  != ""
    ) %>%
    dplyr::distinct(
      !!ctr_sym, !!pt_sym, !!org_sym, !!abx_sym, !!val_sym
    )

  # -- 4. optional: filter to a single centre -------------------------------
  if (mode == "single") {
    abx_clean <- abx_clean %>%
      dplyr::filter(!!ctr_sym == center)
  }

  # -- 5. worst-phenotype rule -----------------------------------------------
  abx_episode <- abx_clean %>%
    dplyr::group_by(!!ctr_sym, !!pt_sym, !!org_sym, !!abx_sym) %>%
    dplyr::summarise(
      final_abx = ifelse(any(!!val_sym == "R"), "R", "S"),
      .groups   = "drop"
    )

  # -- 6. resistance proportion per centre x antibiotic ---------------------
  # Count R and total, then compute proportion R
  abx_counts <- abx_episode %>%
    dplyr::group_by(!!ctr_sym, !!abx_sym, final_abx) %>%
    dplyr::summarise(
      n       = dplyr::n_distinct(!!pt_sym),
      .groups = "drop"
    ) %>%
    dplyr::group_by(!!ctr_sym, !!abx_sym) %>%
    dplyr::mutate(total_tested = sum(n)) %>%
    dplyr::ungroup()

  # Keep only R rows to get resistance proportion
  heat_data <- abx_counts %>%
    dplyr::filter(final_abx == "R") %>%
    dplyr::mutate(
      proportion = n / total_tested,
      value_label = sprintf("%.2f", proportion)
    )

  # -- 7. optional: keep top n antibiotics by total patients tested ----------
  if (!is.null(n)) {
    if (!is.numeric(n) || n < 1) {
      stop("'n' must be a positive integer or NULL.")
    }

    top_abx <- abx_counts %>%
      dplyr::group_by(!!abx_sym) %>%
      dplyr::summarise(grand_total = sum(total_tested), .groups = "drop") %>%
      dplyr::slice_max(grand_total, n = n) %>%
      dplyr::pull(!!abx_sym)

    heat_data <- heat_data %>%
      dplyr::filter(!!abx_sym %in% top_abx)
  }

  # -- 8. order antibiotics by total patients tested (most tested = left) ----
  abx_order <- abx_counts %>%
    dplyr::group_by(!!abx_sym) %>%
    dplyr::summarise(grand_total = sum(total_tested), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(grand_total)) %>%
    dplyr::pull(!!abx_sym)

  # restrict order to antibiotics actually in heat_data after filtering
  abx_order <- abx_order[abx_order %in% unique(heat_data[[antibiotic_col]])]

  heat_data <- heat_data %>%
    dplyr::mutate(
      !!antibiotic_col := factor(!!abx_sym, levels = abx_order)
    )

  # -- 9. build title --------------------------------------------------------
  if (mode == "single") {
    auto_title <- title %||% sprintf(
      "Antibiotic Resistance Proportions \u2014 %s", center
    )
  } else {
    auto_title <- title %||% "Antibiotic Resistance Proportions Across Centres"
  }

  # -- 10. build plot --------------------------------------------------------
  p <- ggplot2::ggplot(
    heat_data,
    ggplot2::aes(
      x    = !!abx_sym,
      y    = !!ctr_sym,
      fill = proportion
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradient2(
      low      = low_colour,
      mid      = mid_colour,
      high     = high_colour,
      midpoint = midpoint,
      limits   = c(0, 1),
      na.value = "grey90",
      name     = "Resistance\nProportion"
    ) +
    ggplot2::scale_x_discrete(
      # grey out missing antibiotics at a centre would be gaps; label only
    ) +
    ggplot2::labs(
      x     = "Antibiotic",
      y     = "Centre",
      title = auto_title
    ) +
    eda_theme(base_size = base_size, legend_position = "right") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # -- 11. optional: value labels inside tiles -------------------------------
  if (show_values) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = value_label),
        size     = 3.2,
        fontface = "bold",
        color    = "white"
      )
  }

  return(p)
}



# FINAL OUTCOME DISTRIBUTION


#' Plot Distribution of Final Outcomes
#'
#' Produces a bar chart of final outcome counts per centre, with each bar
#' labelled with the count and percentage of that centre's (or overall) unique
#' patients. Bars are coloured by outcome category.
#'
#' \strong{Deduplication:} a patient may appear multiple times in a long-format
#' dataset (one row per isolate or antibiotic). This function takes the
#' \emph{first} recorded outcome per patient per centre so each patient is
#' counted once.
#'
#' \strong{Referred merging:} by default, \code{"Transferred to other hospital"}
#' is recoded to \code{"Referred"} since they represent the same clinical event
#' across different centres. Set \code{merge_referred = FALSE} to keep them
#' separate.
#'
#' Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} (default) -- one panel per centre.
#'   \item \strong{"overall"} -- all centres pooled into one chart.
#'   \item \strong{"single"} -- one specific centre; pass the name via
#'     \code{center}.
#' }
#'
#' @param data            Data frame. Long-format AMR dataset.
#' @param mode            Character. One of \code{"faceted"}, \code{"overall"},
#'   or \code{"single"}. Default \code{"faceted"}.
#' @param center          Character. Required when \code{mode = "single"}.
#'   Exact centre name.
#' @param patient_col     Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param outcome_col     Character. Final outcome column.
#'   Default \code{"final_outcome"}.
#' @param center_col      Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param merge_referred  Logical. Recode
#'   \code{"Transferred to other hospital"} to \code{"Referred"}.
#'   Default \code{TRUE}.
#' @param palette         Named character vector mapping outcome values to
#'   colours. Any outcome not in the vector gets a default ggplot colour.
#'   Default uses a preset palette for common outcome categories.
#' @param ncol            Integer. Facet columns (\code{mode = "faceted"} only).
#'   Default 2.
#' @param base_size       Numeric. Base font size. Default 14.
#' @param title           Character. Custom plot title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col    Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name   Character or \code{NULL}. The syndrome value to retain
#'   (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col} to be set.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
#'

plot_outcome_distribution <- function(data,
                                      mode           = c("faceted", "overall", "single"),
                                      center         = NULL,
                                      patient_col    = "PatientInformation_id",
                                      outcome_col    = "final_outcome",
                                      center_col     = "center_name",
                                      merge_referred = TRUE,
                                      palette        = c(
                                        "Death"                        = "#E74C3C",
                                        "Died"                         = "#E74C3C",
                                        "Discharged"                   = "#2ECC71",
                                        "LAMA"                         = "#95A5A6",
                                        "Referred"                     = "#3498DB",
                                        "Transferred to other hospital"= "#3498DB"
                                      ),
                                      ncol           = 2,
                                      base_size      = 14,
                                      title          = NULL,
                                      syndrome_col   = NULL,
                                      syndrome_name  = NULL) {

  # -- 0. match mode ---------------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, outcome_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center)) {
      stop("'center' must be provided when mode = 'single'.")
    }
    available <- unique(data[[center_col]])
    if (!center %in% available) {
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col,
        paste(sort(available), collapse = ", ")
      ))
    }
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  out_sym <- rlang::sym(outcome_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. clean: drop missing outcomes --------------------------------------
  clean <- data %>%
    dplyr::filter(
      !is.na(!!out_sym),
      trimws(as.character(!!out_sym)) != ""
    )

  # -- 4. optionally merge referred categories -------------------------------
  if (merge_referred) {
    clean <- clean %>%
      dplyr::mutate(
        !!outcome_col := dplyr::case_when(
          !!out_sym == "Transferred to other hospital" ~ "Referred",
          TRUE ~ as.character(!!out_sym)
        )
      )
    # refresh symbol after mutate
    out_sym <- rlang::sym(outcome_col)
  }

  # -- 5. optional: filter to one centre ------------------------------------
  if (mode == "single") {
    clean <- clean %>%
      dplyr::filter(!!ctr_sym == center)
  }

  # -- 6. deduplicate: one outcome per patient (per centre) -----------------
  # Takes the first recorded outcome to handle long-format duplicates
  if (mode == "overall") {
    outcome_unique <- clean %>%
      dplyr::group_by(!!pt_sym) %>%
      dplyr::summarise(
        !!outcome_col := dplyr::first(!!out_sym),
        .groups = "drop"
      )
  } else {
    outcome_unique <- clean %>%
      dplyr::group_by(!!ctr_sym, !!pt_sym) %>%
      dplyr::summarise(
        !!outcome_col := dplyr::first(!!out_sym),
        .groups = "drop"
      )
  }

  # -- 7. summarise counts and percentages ----------------------------------
  if (mode == "overall") {

    summary_out <- outcome_unique %>%
      dplyr::group_by(!!out_sym) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        total   = sum(n),
        percent = 100 * n / total,
        label   = paste0(n, " (", round(percent, 1), "%)"),
        !!outcome_col := stats::reorder(!!out_sym, n)
      )

    auto_title <- title %||%
      "Distribution of Final Outcomes \u2014 All Centres Pooled"

    p <- ggplot2::ggplot(
      summary_out,
      ggplot2::aes(x = !!out_sym, y = n, fill = !!out_sym)
    ) +
      ggplot2::geom_col(width = 0.7, color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        vjust    = -0.3,
        size     = 3.6,
        fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::scale_fill_manual(values = palette, na.value = "grey70") +
      ggplot2::labs(
        x     = "Final Outcome",
        y     = "Number of Patients",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none")

    return(p)
  }

  # -- 8. per-centre summary (faceted + single) ------------------------------
  summary_out <- outcome_unique %>%
    dplyr::group_by(!!ctr_sym, !!out_sym) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(!!ctr_sym) %>%
    dplyr::mutate(
      total   = sum(n),
      percent = 100 * n / total,
      label   = paste0(n, " (", round(percent, 1), "%)")
    ) %>%
    dplyr::ungroup()

  # -- 9. SINGLE -------------------------------------------------------------
  if (mode == "single") {

    auto_title <- title %||% sprintf(
      "Distribution of Final Outcomes \u2014 %s", center
    )

    summary_out <- summary_out %>%
      dplyr::mutate(
        !!outcome_col := stats::reorder(!!out_sym, n)
      )

    p <- ggplot2::ggplot(
      summary_out,
      ggplot2::aes(x = !!out_sym, y = n, fill = !!out_sym)
    ) +
      ggplot2::geom_col(width = 0.7, color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        vjust    = -0.3,
        size     = 3.6,
        fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::scale_fill_manual(values = palette, na.value = "grey70") +
      ggplot2::labs(
        x     = "Final Outcome",
        y     = "Number of Patients",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none")

    return(p)
  }

  # -- 10. FACETED -----------------------------------------------------------
  auto_title <- title %||%
    "Distribution of Final Outcomes \u2014 All Centres"

  p <- ggplot2::ggplot(
    summary_out,
    ggplot2::aes(
      x    = tidytext::reorder_within(!!out_sym, n, !!ctr_sym),
      y    = n,
      fill = !!out_sym
    )
  ) +
    ggplot2::geom_col(width = 0.7, color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      vjust    = -0.3,
      size     = 3.6,
      fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    tidytext::scale_x_reordered() +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", center_col)),
      ncol   = ncol,
      scales = "free"
    ) +
    ggplot2::scale_fill_manual(values = palette, na.value = "grey70") +
    ggplot2::labs(
      x     = "Final Outcome",
      y     = "Number of Patients",
      fill  = "Outcome",
      title = auto_title
    ) +
    eda_theme(base_size = base_size, legend_position = "none") +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(2.5, "lines"),
      plot.margin   = ggplot2::margin(15, 20, 15, 20)
    )

  return(p)
}



# FINAL OUTCOME PROPORTIONS BY RESISTANCE STATUS


#' Plot Final Outcome Proportions for Resistant or Susceptible Patients
#'
#' Produces a 100\% horizontal stacked bar chart showing the distribution of
#' final outcomes (Death, Discharged, LAMA, Referred) for the top \code{n}
#' organisms, split by whether the patient had a Resistant or Susceptible
#' result.
#'
#' Each bar represents one organism. The fill shows outcome proportions.
#' A \code{n=} label is printed above each bar showing the total number of
#' unique patients in that organism x resistance group.
#'
#' \strong{Resistance grouping:} \code{resistance_filter = "R"} shows patients
#' who had \emph{at least one} R result for the given organism. \code{"S"} shows
#' patients with only S results. Because these groups overlap (a patient
#' resistant to one antibiotic may be susceptible to another), calling the
#' function twice -- once for "R" and once for "S" -- gives the most complete
#' picture.
#'
#' \strong{Top N organisms} are ranked by total unique patients per centre
#' (or globally for \code{mode = "overall"}), calculated before the R/S
#' split so the same organisms appear in both R and S plots.
#'
#' Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} (default) -- one panel per centre.
#'   \item \strong{"overall"} -- all centres pooled.
#'   \item \strong{"single"} -- one specific centre.
#' }
#'
#' @param data              Data frame. Long-format AMR dataset.
#' @param n                 Integer. Number of top organisms to show.
#'   Default 5.
#' @param resistance_filter Character. \code{"R"} or \code{"S"}. Which patient
#'   group to display. Default \code{"R"}.
#' @param mode              Character. One of \code{"faceted"},
#'   \code{"overall"}, or \code{"single"}. Default \code{"faceted"}.
#' @param center            Character. Required when \code{mode = "single"}.
#'   Exact centre name.
#' @param patient_col       Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param organism_col      Character. Organism name column.
#'   Default \code{"organism_name"}.
#' @param antibiotic_col    Character. Antibiotic name column.
#'   Default \code{"antibiotic_name"}.
#' @param value_col         Character. Susceptibility result column (R/S).
#'   Default \code{"antibiotic_value"}.
#' @param outcome_col       Character. Final outcome column.
#'   Default \code{"final_outcome"}.
#' @param center_col        Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param merge_referred    Logical. Recode
#'   \code{"Transferred to other hospital"} to \code{"Referred"}.
#'   Default \code{TRUE}.
#' @param palette           Named character vector mapping outcome values to
#'   colours. Unmatched outcomes get grey automatically.
#' @param ncol              Integer. Facet columns
#'   (\code{mode = "faceted"} only). Default 2.
#' @param base_size         Numeric. Base font size. Default 14.
#' @param title             Character. Custom title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col      Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name     Character or \code{NULL}. The syndrome value to
#'   retain (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col}.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
#'

plot_outcome_by_organism <- function(data,
                                     n                 = 5,
                                     resistance_filter = c("R", "S"),
                                     mode              = c("faceted", "overall", "single"),
                                     center            = NULL,
                                     patient_col       = "PatientInformation_id",
                                     organism_col      = "organism_name",
                                     antibiotic_col    = "antibiotic_name",
                                     value_col         = "antibiotic_value",
                                     outcome_col       = "final_outcome",
                                     center_col        = "center_name",
                                     merge_referred    = TRUE,
                                     palette           = c(
                                       "Death"                         = "#E74C3C",
                                       "Died"                          = "#E74C3C",
                                       "Discharged"                    = "#2ECC71",
                                       "LAMA"                          = "#95A5A6",
                                       "Referred"                      = "#3498DB",
                                       "Transferred to other hospital" = "#3498DB"
                                     ),
                                     ncol              = 2,
                                     base_size         = 14,
                                     title             = NULL,
                                     syndrome_col      = NULL,
                                     syndrome_name     = NULL) {

  # -- 0. match args ---------------------------------------------------------
  resistance_filter <- match.arg(resistance_filter)
  mode              <- match.arg(mode)

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, organism_col, antibiotic_col,
                     value_col, outcome_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center)) {
      stop("'center' must be provided when mode = 'single'.")
    }
    available <- unique(data[[center_col]])
    if (!center %in% available) {
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col,
        paste(sort(available), collapse = ", ")
      ))
    }
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  org_sym <- rlang::sym(organism_col)
  abx_sym <- rlang::sym(antibiotic_col)
  val_sym <- rlang::sym(value_col)
  out_sym <- rlang::sym(outcome_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. clean: drop blank organisms and outcomes ---------------------------
  clean <- data %>%
    dplyr::filter(
      !is.na(!!org_sym), trimws(as.character(!!org_sym)) != "",
      !is.na(!!out_sym), trimws(as.character(!!out_sym)) != ""
    )

  # -- 4. merge referred categories -----------------------------------------
  if (merge_referred) {
    clean <- clean %>%
      dplyr::mutate(
        !!outcome_col := dplyr::case_when(
          !!out_sym == "Transferred to other hospital" ~ "Referred",
          TRUE ~ as.character(!!out_sym)
        )
      )
    out_sym <- rlang::sym(outcome_col)
  }

  # -- 5. filter to one centre if single mode --------------------------------
  if (mode == "single") {
    clean <- clean %>% dplyr::filter(!!ctr_sym == center)
  }

  # -- 6. identify top N organisms (before R/S split) -----------------------
  # Ranked by total unique patients so both R and S plots show the same organisms
  if (mode == "overall") {
    top_org <- clean %>%
      dplyr::filter(!is.na(!!org_sym), trimws(as.character(!!org_sym)) != "") %>%
      dplyr::group_by(!!org_sym) %>%
      dplyr::summarise(
        unique_patients = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::slice_max(unique_patients, n = n)
  } else {
    top_org <- clean %>%
      dplyr::filter(!is.na(!!org_sym), trimws(as.character(!!org_sym)) != "") %>%
      dplyr::group_by(!!ctr_sym, !!org_sym) %>%
      dplyr::summarise(
        unique_patients = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!ctr_sym) %>%
      dplyr::slice_max(unique_patients, n = n) %>%
      dplyr::ungroup()
  }

  # -- 7. filter to top N organisms + valid R/S rows -------------------------
  join_cols <- if (mode == "overall") organism_col else c(center_col, organism_col)

  df_filtered <- clean %>%
    dplyr::semi_join(top_org, by = join_cols) %>%
    dplyr::filter(
      !!val_sym %in% c("R", "S"),
      !is.na(!!out_sym), trimws(as.character(!!out_sym)) != ""
    ) %>%
    dplyr::distinct(
      !!ctr_sym, !!pt_sym, !!org_sym,
      !!abx_sym, !!val_sym, !!out_sym
    )

  # -- 8. summarise: unique patients per organism x R/S x outcome -----------
  if (mode == "overall") {
    summary_data <- df_filtered %>%
      dplyr::group_by(!!org_sym, !!val_sym, !!out_sym) %>%
      dplyr::summarise(
        unique_patients = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!org_sym, !!val_sym) %>%
      dplyr::mutate(
        total_patients = sum(unique_patients),
        proportion     = unique_patients / total_patients,
        percent        = round(100 * proportion, 1),
        label          = paste0(percent, "%")
      ) %>%
      dplyr::ungroup()

    label_data <- summary_data %>%
      dplyr::distinct(!!org_sym, !!val_sym, total_patients)

  } else {
    summary_data <- df_filtered %>%
      dplyr::group_by(!!ctr_sym, !!org_sym, !!val_sym, !!out_sym) %>%
      dplyr::summarise(
        unique_patients = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!ctr_sym, !!org_sym, !!val_sym) %>%
      dplyr::mutate(
        total_patients = sum(unique_patients),
        proportion     = unique_patients / total_patients,
        percent        = round(100 * proportion, 1),
        label          = paste0(percent, "%")
      ) %>%
      dplyr::ungroup()

    label_data <- summary_data %>%
      dplyr::distinct(!!ctr_sym, !!org_sym, !!val_sym, total_patients)
  }

  # -- 9. filter both tables to the chosen resistance group ------------------
  plot_data  <- summary_data %>% dplyr::filter(!!val_sym == resistance_filter)
  label_data <- label_data   %>% dplyr::filter(!!val_sym == resistance_filter)

  if (nrow(plot_data) == 0) {
    stop(sprintf(
      "No data remaining after filtering to resistance_filter = '%s'. ",
      resistance_filter
    ))
  }

  # -- 10. build auto title --------------------------------------------------
  rs_label <- if (resistance_filter == "R") "Resistant" else "Susceptible"

  if (mode == "overall") {
    auto_title <- title %||% sprintf(
      "Final Outcome Proportions (%s Patients) \u2014 Top %d Organisms, All Centres Pooled",
      rs_label, n
    )
  } else if (mode == "single") {
    auto_title <- title %||% sprintf(
      "Final Outcome Proportions (%s Patients) \u2014 Top %d Organisms, %s",
      rs_label, n, center
    )
  } else {
    auto_title <- title %||% sprintf(
      "Final Outcome Proportions (%s Patients) \u2014 Top %d Organisms",
      rs_label, n
    )
  }

  # -- 11. OVERALL plot ------------------------------------------------------
  if (mode == "overall") {

    plot_data  <- plot_data  %>% dplyr::mutate(!!organism_col := stats::reorder(!!org_sym, total_patients))
    label_data <- label_data %>% dplyr::mutate(!!organism_col := stats::reorder(!!org_sym, total_patients))

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = !!org_sym, y = proportion, fill = !!out_sym)
    ) +
      ggplot2::geom_col(color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white", size = 3.4, fontface = "bold"
      ) +
      ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(x = !!org_sym, y = 1.05,
                     label = paste0("n=", total_patients)),
        inherit.aes = FALSE, size = 3.4, fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        expand = ggplot2::expansion(mult = c(0, 0.12))
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = palette, na.value = "grey70") +
      ggplot2::labs(
        x = "Organism", y = "Proportion of Patients",
        fill = "Final Outcome", title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 12. FACETED + SINGLE share reorder_within logic ----------------------
  plot_data <- plot_data %>%
    dplyr::mutate(
      org_ordered = tidytext::reorder_within(!!org_sym, total_patients, !!ctr_sym)
    )

  label_data <- label_data %>%
    dplyr::mutate(
      org_ordered = tidytext::reorder_within(!!org_sym, total_patients, !!ctr_sym)
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = org_ordered, y = proportion, fill = !!out_sym)
  ) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white", size = 3.4, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = org_ordered, y = 1.05,
                   label = paste0("n=", total_patients)),
      inherit.aes = FALSE, size = 3.4, fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    tidytext::scale_x_reordered() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = palette, na.value = "grey70") +
    ggplot2::labs(
      x = "Organism", y = "Proportion of Patients",
      fill = "Final Outcome", title = auto_title
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.5, "lines"))

  # add facet only for faceted mode
  if (mode == "faceted") {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", center_col)),
      ncol = ncol, scales = "free_y"
    )
  }

  return(p)
}



# DEATH VS DISCHARGED -- SIDE-BY-SIDE BAR CHART


#' Plot Death vs Discharged Counts for Top Pathogens
#'
#' Produces a side-by-side (dodged) bar chart comparing Death and Discharged
#' patient counts for the top \code{n} organisms. LAMA and Referred are
#' excluded -- only the two definitive outcomes are shown.
#'
#' Each bar is labelled with the count and the percentage that outcome
#' represents out of all Death + Discharged patients for that organism
#' (e.g. \code{45 (32\%)} on a Death bar means 45 patients died, which is
#' 32\% of the organism's known outcomes). This makes the mortality rate
#' immediately readable without a separate calculation.
#'
#' \strong{Deduplication:} uses
#' \code{distinct(centre, patient, organism, outcome)} -- one record per
#' patient per organism per outcome -- consistent with the original EDA script.
#'
#' \strong{Top N organisms} are ranked by total unique patients with a known
#' outcome (Death or Discharged). For \code{mode = "faceted"} the ranking is
#' done independently per centre; for \code{"overall"} it is done globally.
#'
#' Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} (default) -- one panel per centre.
#'   \item \strong{"overall"} -- all centres pooled into one chart.
#'   \item \strong{"single"} -- one specific centre.
#' }
#'
#' @param data              Data frame. Long-format AMR dataset.
#' @param n                 Integer. Number of top organisms to show.
#'   Default 5.
#' @param mode              Character. One of \code{"faceted"},
#'   \code{"overall"}, or \code{"single"}. Default \code{"faceted"}.
#' @param center            Character. Required when \code{mode = "single"}.
#'   Exact centre name.
#' @param resistance_filter Character or \code{NULL}. When \code{"R"}, only
#'   patients who had at least one Resistant result for the organism are
#'   included (worst-phenotype rule: any R = R). When \code{"S"}, only
#'   patients with exclusively Susceptible results are included. \code{NULL}
#'   (default) includes all patients regardless of resistance status.
#' @param patient_col       Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param organism_col      Character. Organism name column.
#'   Default \code{"organism_name"}.
#' @param outcome_col       Character. Final outcome column.
#'   Default \code{"final_outcome"}.
#' @param center_col        Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param antibiotic_col    Character. Antibiotic name column. Only used when
#'   \code{resistance_filter} is not \code{NULL}.
#'   Default \code{"antibiotic_name"}.
#' @param value_col         Character. Susceptibility result column (R/S).
#'   Only used when \code{resistance_filter} is not \code{NULL}.
#'   Default \code{"antibiotic_value"}.
#' @param death_label       Character. Exact string used for death in the
#'   outcome column. Default \code{"Death"}.
#' @param discharged_label  Character. Exact string used for discharged.
#'   Default \code{"Discharged"}.
#' @param colours           Named character vector with keys matching
#'   \code{death_label} and \code{discharged_label}.
#'   Default \code{c("Death" = "#E74C3C", "Discharged" = "#2ECC71")}.
#' @param bar_width         Numeric. Width of each individual bar (0-1).
#'   Default \code{0.65}.
#' @param ncol              Integer. Facet columns
#'   (\code{mode = "faceted"} only). Default 2.
#' @param base_size         Numeric. Base font size. Default 14.
#' @param title             Character. Custom title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col      Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name     Character or \code{NULL}. The syndrome value to
#'   retain (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col}.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
#'

plot_death_discharged <- function(data,
                                  n                = 5,
                                  mode             = c("faceted", "overall", "single"),
                                  center           = NULL,
                                  resistance_filter = NULL,
                                  patient_col      = "PatientInformation_id",
                                  organism_col     = "organism_name",
                                  outcome_col      = "final_outcome",
                                  center_col       = "center_name",
                                  antibiotic_col   = "antibiotic_name",
                                  value_col        = "antibiotic_value",
                                  death_label      = "Death",
                                  discharged_label = "Discharged",
                                  colours          = NULL,
                                  bar_width        = 0.65,
                                  ncol             = 2,
                                  base_size        = 14,
                                  title            = NULL,
                                  syndrome_col     = NULL,
                                  syndrome_name    = NULL) {

  # -- 0. match mode ---------------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, organism_col, outcome_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center)) stop("'center' must be provided when mode = 'single'.")
    available <- unique(data[[center_col]])
    if (!center %in% available) {
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col, paste(sort(available), collapse = ", ")
      ))
    }
  }

  if (!is.null(resistance_filter)) {
    if (!resistance_filter %in% c("R", "S")) {
      stop("'resistance_filter' must be \"R\", \"S\", or NULL.")
    }
    missing_abx <- setdiff(c(antibiotic_col, value_col), names(data))
    if (length(missing_abx) > 0) {
      stop(sprintf(
        "resistance_filter requires column(s) not found in data: %s",
        paste(missing_abx, collapse = ", ")
      ))
    }
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. resolve colours ----------------------------------------------------
  if (is.null(colours)) {
    colours <- stats::setNames(
      c("#E74C3C", "#2ECC71"),
      c(death_label, discharged_label)
    )
  }

  # -- 3. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  org_sym <- rlang::sym(organism_col)
  out_sym <- rlang::sym(outcome_col)
  ctr_sym <- rlang::sym(center_col)
  abx_sym <- rlang::sym(antibiotic_col)
  val_sym <- rlang::sym(value_col)

  # -- 4. keep only Death and Discharged, drop blank organisms --------------
  clean <- data %>%
    dplyr::filter(
      !!out_sym %in% c(death_label, discharged_label),
      !is.na(!!org_sym), trimws(as.character(!!org_sym)) != ""
    )

  if (nrow(clean) == 0) {
    stop(sprintf(
      "No rows matched death_label = '%s' or discharged_label = '%s'. ",
      death_label, discharged_label
    ))
  }

  # -- 4b. resistance pre-filter (worst-phenotype at organism level) ---------
  # Classify each patient x organism as R (any R result) or S (all S results)
  # from the FULL dataset, then keep only the desired group.
  if (!is.null(resistance_filter)) {
    resist_class <- data %>%
      dplyr::filter(
        !!val_sym %in% c("R", "S"),
        !is.na(!!org_sym), trimws(as.character(!!org_sym)) != ""
      ) %>%
      dplyr::distinct(!!ctr_sym, !!pt_sym, !!org_sym, !!abx_sym, !!val_sym) %>%
      dplyr::group_by(!!ctr_sym, !!pt_sym, !!org_sym) %>%
      dplyr::summarise(
        organism_resist = ifelse(any(!!val_sym == "R"), "R", "S"),
        .groups = "drop"
      ) %>%
      dplyr::filter(organism_resist == resistance_filter)

    clean <- clean %>%
      dplyr::semi_join(
        resist_class,
        by = c(center_col, patient_col, organism_col)
      )

    if (nrow(clean) == 0) {
      stop(sprintf(
        "No patients remaining after applying resistance_filter = '%s'.",
        resistance_filter
      ))
    }
  }

  # -- 5. filter to one centre if single mode --------------------------------
  if (mode == "single") {
    clean <- clean %>% dplyr::filter(!!ctr_sym == center)
  }

  # -- 6. deduplicate: one row per patient x organism x outcome (x centre) --
  if (mode == "overall") {
    clean <- clean %>%
      dplyr::distinct(!!pt_sym, !!org_sym, !!out_sym)
  } else {
    clean <- clean %>%
      dplyr::distinct(!!ctr_sym, !!pt_sym, !!org_sym, !!out_sym)
  }

  # -- 7. identify top N organisms by total known-outcome patients -----------
  if (mode == "overall") {
    top_org <- clean %>%
      dplyr::group_by(!!org_sym) %>%
      dplyr::summarise(
        total_patients = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::slice_max(total_patients, n = n)
  } else {
    top_org <- clean %>%
      dplyr::group_by(!!ctr_sym, !!org_sym) %>%
      dplyr::summarise(
        total_patients = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!ctr_sym) %>%
      dplyr::slice_max(total_patients, n = n) %>%
      dplyr::ungroup()
  }

  join_cols <- if (mode == "overall") organism_col else c(center_col, organism_col)

  clean <- clean %>%
    dplyr::semi_join(top_org, by = join_cols)

  # -- 8. count unique patients per organism x outcome (x centre) -----------
  if (mode == "overall") {
    summary_data <- clean %>%
      dplyr::group_by(!!org_sym, !!out_sym) %>%
      dplyr::summarise(
        n       = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!org_sym) %>%
      dplyr::mutate(
        total   = sum(n),
        percent = round(100 * n / total, 1),
        label   = paste0(n, "\n(", percent, "%)")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        !!organism_col := stats::reorder(!!org_sym, total)
      )
  } else {
    summary_data <- clean %>%
      dplyr::group_by(!!ctr_sym, !!org_sym, !!out_sym) %>%
      dplyr::summarise(
        n       = dplyr::n_distinct(!!pt_sym),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!ctr_sym, !!org_sym) %>%
      dplyr::mutate(
        total   = sum(n),
        percent = round(100 * n / total, 1),
        label   = paste0(n, "\n(", percent, "%)")
      ) %>%
      dplyr::ungroup()
  }

  # -- 9. build auto title ---------------------------------------------------
  rs_tag <- if (!is.null(resistance_filter)) {
    if (resistance_filter == "R") " (Resistant Patients)" else " (Susceptible Patients)"
  } else ""

  if (mode == "overall") {
    auto_title <- title %||% sprintf(
      "Death vs Discharged%s \u2014 Top %d Organisms, All Centres Pooled",
      rs_tag, n
    )
  } else if (mode == "single") {
    auto_title <- title %||% sprintf(
      "Death vs Discharged%s \u2014 Top %d Organisms, %s",
      rs_tag, n, center
    )
  } else {
    auto_title <- title %||% sprintf(
      "Death vs Discharged%s \u2014 Top %d Organisms",
      rs_tag, n
    )
  }

  # -- 10. OVERALL plot ------------------------------------------------------
  if (mode == "overall") {

    p <- ggplot2::ggplot(
      summary_data,
      ggplot2::aes(
        x    = !!org_sym,
        y    = n,
        fill = !!out_sym
      )
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = bar_width + 0.05),
        width    = bar_width,
        color    = "black"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_dodge(width = bar_width + 0.05),
        vjust    = -0.3,
        size     = 3.2,
        fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.18))
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = colours) +
      ggplot2::labs(
        x     = "Organism",
        y     = "Number of Patients",
        fill  = "Outcome",
        title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 11. FACETED + SINGLE share reorder_within -----------------------------
  summary_data <- summary_data %>%
    dplyr::mutate(
      org_ordered = tidytext::reorder_within(!!org_sym, total, !!ctr_sym)
    )

  p <- ggplot2::ggplot(
    summary_data,
    ggplot2::aes(
      x    = org_ordered,
      y    = n,
      fill = !!out_sym
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = bar_width + 0.05),
      width    = bar_width,
      color    = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_dodge(width = bar_width + 0.05),
      vjust    = -0.3,
      size     = 3.2,
      fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.22))
    ) +
    tidytext::scale_x_reordered() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::labs(
      x     = "Organism",
      y     = "Number of Patients",
      fill  = "Outcome",
      title = auto_title
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.5, "lines"))

  if (mode == "faceted") {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", center_col)),
      ncol   = ncol,
      scales = "free"
    )
  }

  return(p)
}



# ANTIBIOTIC RESISTANCE DISTRIBUTION ACROSS SAMPLE TYPES


#' Plot Distribution of Antibiotic Resistance Across Sample Types
#'
#' Produces a horizontal 100\% stacked bar chart showing the proportion of
#' Resistant (R) and Susceptible (S) deduplicated AST results for the top
#' \code{n} sample types. Each bar sums to 100\% and is labelled with the R
#' and S percentages.
#'
#' \strong{Deduplication (worst-phenotype):} before counting, rows are
#' collapsed by patient x centre x sample type x organism x antibiotic. If
#' the same antibiotic appears more than once for a patient-organism-sample
#' combination (data entry error or repeat test), any single R result marks
#' the episode as R. This prevents duplicate rows from inflating counts.
#'
#' \strong{Counting unit:} deduplicated AST test results -- one record per
#' patient x organism x antibiotic x sample type combination.
#'
#' \strong{Top N ranking:} sample types are ranked by total deduplicated
#' test count, so the most-tested specimen types appear first.
#'
#' \strong{Default sample column:} \code{"sample_type"} -- the raw column.
#' Switch to \code{"specimen_normalized"} for standardised values.
#'
#' Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} (default) -- top \emph{n} sample types per
#'     centre independently, one panel per centre.
#'   \item \strong{"overall"} -- all centres pooled; top \emph{n} globally.
#'   \item \strong{"single"} -- one specific centre.
#' }
#'
#' @param data           Data frame. Long-format AMR dataset.
#' @param n              Integer. Number of top sample types to show (ranked
#'   by deduplicated test count). Default 10.
#' @param mode           Character. One of \code{"faceted"}, \code{"overall"},
#'   or \code{"single"}. Default \code{"faceted"}.
#' @param center         Character. Required when \code{mode = "single"}.
#'   Exact centre name.
#' @param patient_col    Character. Patient ID column. Used for deduplication.
#'   Default \code{"PatientInformation_id"}.
#' @param sample_col     Character. Sample/specimen type column.
#'   Default \code{"sample_type"}.
#' @param organism_col   Character. Organism name column. Used for
#'   deduplication. Default \code{"organism_name"}.
#' @param antibiotic_col Character. Antibiotic name column.
#'   Default \code{"antibiotic_name"}.
#' @param value_col      Character. Susceptibility result column (R/S values).
#'   Default \code{"antibiotic_value"}.
#' @param center_col     Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param colours        Named character vector for R and S bars.
#'   Default \code{c("R" = "#D73027", "S" = "#1A9850")}.
#' @param ncol           Integer. Facet columns (\code{mode = "faceted"} only).
#'   Default 2.
#' @param base_size      Numeric. Base font size. Default 14.
#' @param title          Character. Custom title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col   Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name  Character or \code{NULL}. The syndrome value to retain
#'   (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col} to be set.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
#'

plot_resistance_by_sample <- function(data,
                                      n              = 10,
                                      mode           = c("faceted", "overall", "single"),
                                      center         = NULL,
                                      patient_col    = "PatientInformation_id",
                                      sample_col     = "sample_type",
                                      organism_col   = "organism_name",
                                      antibiotic_col = "antibiotic_name",
                                      value_col      = "antibiotic_value",
                                      center_col     = "center_name",
                                      colours        = c("R" = "#D73027",
                                                         "S" = "#1A9850"),
                                      ncol           = 2,
                                      base_size      = 14,
                                      title          = NULL,
                                      syndrome_col   = NULL,
                                      syndrome_name  = NULL) {

  # -- 0. match mode ---------------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, sample_col, organism_col, antibiotic_col, value_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center)) {
      stop("'center' must be provided when mode = 'single'.")
    }
    available <- unique(data[[center_col]])
    if (!center %in% available) {
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col,
        paste(sort(available), collapse = ", ")
      ))
    }
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pat_sym <- rlang::sym(patient_col)
  smp_sym <- rlang::sym(sample_col)
  org_sym <- rlang::sym(organism_col)
  abx_sym <- rlang::sym(antibiotic_col)
  val_sym <- rlang::sym(value_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. clean: valid R/S rows, drop blank sample types and antibiotics -----
  clean <- data %>%
    dplyr::filter(
      !!val_sym %in% c("R", "S"),
      !is.na(!!smp_sym), trimws(as.character(!!smp_sym)) != "",
      !is.na(!!abx_sym), trimws(as.character(!!abx_sym)) != ""
    )

  # -- 4. filter to one centre if single mode --------------------------------
  if (mode == "single") {
    clean <- clean %>% dplyr::filter(!!ctr_sym == center)
  }

  # -- 5. deduplicate: worst-phenotype per patient x centre x sample x
  #       organism x antibiotic -- collapses duplicate/conflicting rows -------
  dedup_cols <- c(center_col, patient_col, sample_col, organism_col, antibiotic_col)
  if (mode == "overall") dedup_cols <- setdiff(dedup_cols, center_col)

  clean <- clean %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(dedup_cols))) %>%
    dplyr::summarise(
      !!value_col := if (any(!!val_sym == "R")) "R" else "S",
      .groups = "drop"
    )

  # -- 6. identify top N sample types (by deduplicated test count) ----------
  if (mode == "overall") {
    top_smp <- clean %>%
      dplyr::group_by(!!smp_sym) %>%
      dplyr::summarise(tests = dplyr::n(), .groups = "drop") %>%
      dplyr::slice_max(tests, n = n)
  } else {
    top_smp <- clean %>%
      dplyr::group_by(!!ctr_sym, !!smp_sym) %>%
      dplyr::summarise(tests = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(!!ctr_sym) %>%
      dplyr::slice_max(tests, n = n) %>%
      dplyr::ungroup()
  }

  join_cols <- if (mode == "overall") sample_col else c(center_col, sample_col)
  clean <- clean %>% dplyr::semi_join(top_smp, by = join_cols)

  # -- 7. summarise: deduplicated counts per sample type x R/S (x centre) ----
  if (mode == "overall") {
    summary_data <- clean %>%
      dplyr::group_by(!!smp_sym, !!val_sym) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(!!smp_sym) %>%
      dplyr::mutate(
        total_tests = sum(n),
        proportion  = n / total_tests,
        label       = paste0(round(100 * proportion, 1), "%")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(!!sample_col := stats::reorder(!!smp_sym, total_tests))

    auto_title <- title %||% sprintf(
      "Antibiotic Resistance by Sample Type (Top %d) \u2014 All Centres Pooled", n
    )

    p <- ggplot2::ggplot(
      summary_data,
      ggplot2::aes(x = !!smp_sym, y = proportion, fill = !!val_sym)
    ) +
      ggplot2::geom_col(width = 0.7, color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white", size = 3.4, fontface = "bold"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      ggplot2::scale_fill_manual(
        values = colours,
        labels = c("R" = "Resistant", "S" = "Susceptible")
      ) +
      ggplot2::labs(
        x     = "Sample Type",
        y     = "Proportion of AST Results",
        fill  = "Susceptibility",
        title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 8. per-centre summary (faceted + single) ------------------------------
  summary_data <- clean %>%
    dplyr::group_by(!!ctr_sym, !!smp_sym, !!val_sym) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(!!ctr_sym, !!smp_sym) %>%
    dplyr::mutate(
      total_tests = sum(n),
      proportion  = n / total_tests,
      label       = paste0(round(100 * proportion, 1), "%")
    ) %>%
    dplyr::ungroup()

  top_smp_with_total <- summary_data %>%
    dplyr::distinct(!!ctr_sym, !!smp_sym, total_tests) %>%
    dplyr::group_by(!!ctr_sym) %>%
    dplyr::slice_max(total_tests, n = n) %>%
    dplyr::ungroup()

  summary_data <- summary_data %>%
    dplyr::inner_join(top_smp_with_total,
                      by = c(center_col, sample_col, "total_tests"))

  # -- 9. SINGLE -------------------------------------------------------------
  if (mode == "single") {
    auto_title <- title %||% sprintf(
      "Antibiotic Resistance by Sample Type (Top %d) \u2014 %s", n, center
    )

    summary_data <- summary_data %>%
      dplyr::mutate(!!sample_col := stats::reorder(!!smp_sym, total_tests))

    p <- ggplot2::ggplot(
      summary_data,
      ggplot2::aes(x = !!smp_sym, y = proportion, fill = !!val_sym)
    ) +
      ggplot2::geom_col(width = 0.7, color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white", size = 3.4, fontface = "bold"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      ggplot2::scale_fill_manual(
        values = colours,
        labels = c("R" = "Resistant", "S" = "Susceptible")
      ) +
      ggplot2::labs(
        x     = "Sample Type",
        y     = "Proportion of AST Results",
        fill  = "Susceptibility",
        title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 10. FACETED -----------------------------------------------------------
  auto_title <- title %||% sprintf(
    "Antibiotic Resistance by Sample Type (Top %d) \u2014 All Centres", n
  )

  p <- ggplot2::ggplot(
    summary_data,
    ggplot2::aes(
      x    = tidytext::reorder_within(!!smp_sym, total_tests, !!ctr_sym),
      y    = proportion,
      fill = !!val_sym
    )
  ) +
    ggplot2::geom_col(width = 0.7, color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white", size = 3.4, fontface = "bold"
    ) +
    ggplot2::coord_flip() +
    tidytext::scale_x_reordered() +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", center_col)),
      ncol   = ncol,
      scales = "free_y"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_fill_manual(
      values = colours,
      labels = c("R" = "Resistant", "S" = "Susceptible")
    ) +
    ggplot2::labs(
      x     = "Sample Type",
      y     = "Proportion of AST Results",
      fill  = "Susceptibility",
      title = auto_title
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.5, "lines"))

  return(p)
}



# FINAL OUTCOME BY AGE BIN


#' Plot Final Outcome Proportions by Age Bin
#'
#' Produces a horizontal 100\% stacked bar chart showing how final outcomes
#' (Death, Discharged, LAMA, Referred) are distributed across age groups.
#' Each bar represents one age bin and sums to 100\%. An \code{n=} label above
#' each bar shows the total patients in that age group.
#'
#' \strong{Deduplication:} one record per patient per centre -- first recorded
#' outcome and age bin per patient per centre. Each patient contributes exactly
#' one bar segment per centre.
#'
#' \strong{Age bin ordering:} bins are ordered numerically (\code{<1},
#' \code{1-5}, \code{5-10}, \ldots, \code{85+}) automatically. Supply
#' \code{age_levels} to override with a custom order.
#'
#' Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} (default) -- one panel per centre.
#'   \item \strong{"overall"} -- all centres pooled; one chart.
#'   \item \strong{"single"} -- one specific centre.
#' }
#'
#' @param data           Data frame. Long-format AMR dataset.
#' @param mode           Character. One of \code{"faceted"}, \code{"overall"},
#'   or \code{"single"}. Default \code{"faceted"}.
#' @param center         Character. Required when \code{mode = "single"}.
#'   Exact centre name.
#' @param patient_col    Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param agebin_col     Character. Age bin column. Default \code{"Age_bin"}.
#' @param outcome_col    Character. Final outcome column.
#'   Default \code{"final_outcome"}.
#' @param center_col     Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param merge_referred Logical. Recode
#'   \code{"Transferred to other hospital"} to \code{"Referred"}.
#'   Default \code{TRUE}.
#' @param age_levels     Character vector or \code{NULL}. Custom ordering of
#'   age bin labels. If \code{NULL} (default), bins are sorted automatically
#'   (\code{<1} first, numeric ranges in order, open-ended bins last).
#' @param palette        Named character vector. Outcome -> fill colour mapping.
#'   Unmatched outcomes receive grey. Default uses red/green/grey/blue preset.
#' @param ncol           Integer. Facet columns (\code{mode = "faceted"} only).
#'   Default 2.
#' @param base_size      Numeric. Base font size. Default 14.
#' @param title          Character. Custom title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col   Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name  Character or \code{NULL}. The syndrome value to retain
#'   (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col} to be set.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
#'

plot_outcome_by_agebin <- function(data,
                                   mode           = c("faceted", "overall", "single"),
                                   center         = NULL,
                                   patient_col    = "PatientInformation_id",
                                   agebin_col     = "Age_bin",
                                   outcome_col    = "final_outcome",
                                   center_col     = "center_name",
                                   merge_referred = TRUE,
                                   age_levels     = NULL,
                                   palette        = c(
                                     "Death"                         = "#E74C3C",
                                     "Died"                          = "#E74C3C",
                                     "Discharged"                    = "#2ECC71",
                                     "LAMA"                          = "#95A5A6",
                                     "Referred"                      = "#3498DB",
                                     "Transferred to other hospital" = "#3498DB"
                                   ),
                                   ncol           = 2,
                                   base_size      = 14,
                                   title          = NULL,
                                   syndrome_col   = NULL,
                                   syndrome_name  = NULL) {

  # -- 0. match mode ---------------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, agebin_col, outcome_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (mode == "single") {
    if (is.null(center))
      stop("'center' must be provided when mode = 'single'.")
    available <- unique(data[[center_col]])
    if (!center %in% available)
      stop(sprintf(
        "'%s' not found in column '%s'. Available values: %s",
        center, center_col, paste(sort(available), collapse = ", ")
      ))
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  age_sym <- rlang::sym(agebin_col)
  out_sym <- rlang::sym(outcome_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. clean: drop missing age bins and outcomes --------------------------
  clean <- data %>%
    dplyr::filter(
      !is.na(!!age_sym), trimws(as.character(!!age_sym)) != "",
      !is.na(!!out_sym), trimws(as.character(!!out_sym)) != ""
    )

  # -- 4. merge referred categories -----------------------------------------
  if (merge_referred) {
    clean <- clean %>%
      dplyr::mutate(
        !!outcome_col := dplyr::case_when(
          !!out_sym == "Transferred to other hospital" ~ "Referred",
          TRUE ~ as.character(!!out_sym)
        )
      )
    out_sym <- rlang::sym(outcome_col)
  }

  # -- 5. filter to single centre --------------------------------------------
  if (mode == "single") {
    clean <- clean %>% dplyr::filter(!!ctr_sym == center)
  }

  # -- 6. deduplicate: one row per patient (x centre) -----------------------
  if (mode == "overall") {
    clean <- clean %>%
      dplyr::group_by(!!pt_sym) %>%
      dplyr::summarise(
        !!agebin_col  := dplyr::first(!!age_sym),
        !!outcome_col := dplyr::first(!!out_sym),
        .groups = "drop"
      )
  } else {
    clean <- clean %>%
      dplyr::group_by(!!ctr_sym, !!pt_sym) %>%
      dplyr::summarise(
        !!agebin_col  := dplyr::first(!!age_sym),
        !!outcome_col := dplyr::first(!!out_sym),
        .groups = "drop"
      )
    age_sym <- rlang::sym(agebin_col)
    out_sym <- rlang::sym(outcome_col)
  }

  # -- 7. order age bins -----------------------------------------------------
  if (!is.null(age_levels)) {
    lvls <- age_levels
  } else {
    all_bins <- unique(clean[[agebin_col]])
    lt_bins   <- all_bins[grepl("^<", all_bins)]
    plus_bins <- all_bins[grepl("\\+$", all_bins)]
    range_bins <- setdiff(all_bins, c(lt_bins, plus_bins))
    lower_num  <- suppressWarnings(as.numeric(sub("[--].*", "", range_bins)))
    range_bins <- range_bins[order(lower_num)]
    lvls <- c(sort(lt_bins), range_bins, sort(plus_bins))
  }

  clean <- clean %>%
    dplyr::mutate(!!agebin_col := factor(!!age_sym, levels = lvls))

  # -- 8. summarise proportions ----------------------------------------------
  if (mode == "overall") {
    summary_data <- clean %>%
      dplyr::group_by(!!age_sym, !!out_sym) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(!!age_sym) %>%
      dplyr::mutate(
        total      = sum(count),
        proportion = count / total,
        label      = scales::percent(proportion, accuracy = 1)
      ) %>%
      dplyr::ungroup()

    n_labels <- summary_data %>%
      dplyr::distinct(!!age_sym, total)

  } else {
    summary_data <- clean %>%
      dplyr::group_by(!!ctr_sym, !!age_sym, !!out_sym) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(!!ctr_sym, !!age_sym) %>%
      dplyr::mutate(
        total      = sum(count),
        proportion = count / total,
        label      = scales::percent(proportion, accuracy = 1)
      ) %>%
      dplyr::ungroup()

    n_labels <- summary_data %>%
      dplyr::distinct(!!ctr_sym, !!age_sym, total)
  }

  # -- 9. build fill colours -------------------------------------------------
  outcomes_in_data <- unique(as.character(summary_data[[outcome_col]]))
  fill_vals <- ifelse(
    outcomes_in_data %in% names(palette),
    palette[outcomes_in_data],
    "#AAAAAA"
  )
  names(fill_vals) <- outcomes_in_data

  # -- 10. OVERALL -----------------------------------------------------------
  if (mode == "overall") {
    auto_title <- title %||%
      "Proportion of Final Outcomes by Age Bin \u2014 All Centres Pooled"

    p <- ggplot2::ggplot(
      summary_data,
      ggplot2::aes(x = !!age_sym, y = proportion, fill = !!out_sym)
    ) +
      ggplot2::geom_col(color = "black", width = 0.8) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3.2, color = "white", fontface = "bold"
      ) +
      ggplot2::geom_text(
        data = n_labels,
        ggplot2::aes(x = !!age_sym, y = 1.05,
                     label = paste0("n=", total)),
        inherit.aes = FALSE,
        size = 3.4, fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        expand = ggplot2::expansion(mult = c(0, 0.12))
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = fill_vals) +
      ggplot2::labs(
        x     = "Age Bin (years)",
        y     = "Proportion of Patients",
        fill  = "Final Outcome",
        title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 11. SINGLE ------------------------------------------------------------
  if (mode == "single") {
    auto_title <- title %||% sprintf(
      "Proportion of Final Outcomes by Age Bin \u2014 %s", center
    )

    p <- ggplot2::ggplot(
      summary_data,
      ggplot2::aes(x = !!age_sym, y = proportion, fill = !!out_sym)
    ) +
      ggplot2::geom_col(color = "black", width = 0.8) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3.2, color = "white", fontface = "bold"
      ) +
      ggplot2::geom_text(
        data = n_labels,
        ggplot2::aes(x = !!age_sym, y = 1.05,
                     label = paste0("n=", total)),
        inherit.aes = FALSE,
        size = 3.4, fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        expand = ggplot2::expansion(mult = c(0, 0.12))
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = fill_vals) +
      ggplot2::labs(
        x     = "Age Bin (years)",
        y     = "Proportion of Patients",
        fill  = "Final Outcome",
        title = auto_title
      ) +
      eda_theme(base_size = base_size)

    return(p)
  }

  # -- 12. FACETED -----------------------------------------------------------
  auto_title <- title %||%
    "Proportion of Final Outcomes by Age Bin \u2014 All Centres"

  p <- ggplot2::ggplot(
    summary_data,
    ggplot2::aes(x = !!age_sym, y = proportion, fill = !!out_sym)
  ) +
    ggplot2::geom_col(color = "black", width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3.2, color = "white", fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = n_labels,
      ggplot2::aes(x = !!age_sym, y = 1.05,
                   label = paste0("n=", total)),
      inherit.aes = FALSE,
      size = 3.4, fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", center_col)),
      ncol   = ncol,
      scales = "free_y"
    ) +
    ggplot2::scale_fill_manual(values = fill_vals) +
    ggplot2::labs(
      x     = "Age Bin (years)",
      y     = "Proportion of Patients",
      fill  = "Final Outcome",
      title = auto_title
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.5, "lines"))

  return(p)
}



# MONO VS POLYMICROBIAL BY FACILITY


#' Plot Mono vs Polymicrobial Infections by Facility
#'
#' Produces a dodged bar chart comparing the number of unique patients with
#' Monomicrobial vs Polymicrobial infections at each centre. Hospitals appear
#' on the x-axis; the two bars per hospital are coloured by infection type.
#' Count labels are printed above each bar.
#'
#' \strong{Deduplication:} \code{distinct(patient, centre, infection type)} --
#' each unique patient-centre-infection type combination is counted once. A
#' patient who has both mono and polymicrobial episodes at the same centre
#' contributes to both bars, matching the original EDA script.
#'
#' \strong{Infection type:} the \code{poly_col} column is expected to contain
#' binary values (\code{0} = Monomicrobial, \code{1} = Polymicrobial). The
#' labels shown in the legend and plot are controlled by \code{mono_label}
#' and \code{poly_label}.
#'
#' @param data         Data frame. Long-format AMR dataset.
#' @param patient_col  Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param poly_col     Character. Column with binary infection type
#'   (\code{0} = mono, \code{1} = poly). Default \code{"is_polymicrobial"}.
#' @param center_col   Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param mono_label   Character. Label for monomicrobial bars.
#'   Default \code{"Monomicrobial"}.
#' @param poly_label   Character. Label for polymicrobial bars.
#'   Default \code{"Polymicrobial"}.
#' @param colours      Named character vector. Fill colours keyed by
#'   \code{mono_label} and \code{poly_label}.
#'   Default blue/red.
#' @param bar_width    Numeric. Width of each individual bar. Default \code{0.68}.
#' @param base_size    Numeric. Base font size. Default \code{14}.
#' @param title        Character. Custom title. Auto-generated if \code{NULL}.
#' @param syndrome_col  Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels (e.g. \code{"infectious_syndrome"}).
#'   If \code{NULL} (default), no syndrome filtering is applied.
#' @param syndrome_name Character or \code{NULL}. The syndrome value to retain
#'   (e.g. \code{"BSI"}, \code{"VAP"}). Requires \code{syndrome_col} to be set.
#'   If \code{NULL} (default), all syndromes are included.
#'
#' @return A \code{ggplot} object.
#' @export
plot_mono_poly_by_facility <- function(data,
                                       patient_col  = "PatientInformation_id",
                                       poly_col     = "is_polymicrobial",
                                       center_col   = "center_name",
                                       mono_label   = "Monomicrobial",
                                       poly_label   = "Polymicrobial",
                                       colours      = NULL,
                                       bar_width    = 0.68,
                                       base_size    = 14,
                                       title        = NULL,
                                       syndrome_col  = NULL,
                                       syndrome_name = NULL) {

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, poly_col, center_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s",
                 paste(missing_cols, collapse = ", ")))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  pol_sym <- rlang::sym(poly_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. recode 0/1 -> infection_type factor --------------------------------
  clean <- data %>%
    dplyr::filter(!is.na(!!pol_sym)) %>%
    dplyr::mutate(
      infection_type = factor(
        dplyr::if_else(!!pol_sym == 1L | !!pol_sym == TRUE,
                       poly_label, mono_label),
        levels = c(mono_label, poly_label)
      )
    )

  # -- 4. deduplicate: one row per patient x centre x infection type ---------
  plot_data <- clean %>%
    dplyr::distinct(!!pt_sym, !!ctr_sym, infection_type) %>%
    dplyr::count(!!ctr_sym, infection_type, name = "n_patients")

  # -- 5. resolve colours ----------------------------------------------------
  if (is.null(colours)) {
    colours <- stats::setNames(
      c("#2196F3", "#F44336"),
      c(mono_label, poly_label)
    )
  }

  # -- 6. build plot ---------------------------------------------------------
  auto_title <- title %||%
    "Mono vs Polymicrobial Infections by Facility"

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = !!ctr_sym, y = n_patients, fill = infection_type)
  ) +
    ggplot2::geom_col(
      position  = ggplot2::position_dodge(width = bar_width + 0.04),
      width     = bar_width,
      alpha     = 0.88
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_patients),
      position  = ggplot2::position_dodge(width = bar_width + 0.04),
      vjust     = -0.45,
      size      = 3.4,
      fontface  = "bold"
    ) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      x     = "Facility",
      y     = "Number of Unique Patients",
      fill  = "Infection Type",
      title = auto_title
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1,
                                                  face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )

  return(p)
}



# HAI vs CAI BY FACILITY


#' Plot HAI vs CAI Infection Distribution by Facility
#'
#' Derives Healthcare-Associated Infection (HAI) vs Community-Associated
#' Infection (CAI) status from dates, then produces a stacked bar chart
#' showing unique patient counts per centre. Each bar carries count and
#' percentage labels inside each segment, and an \code{n=} total label
#' above.
#'
#' \strong{Classification rule:} for each patient x centre, the earliest
#' culture date is used. If the time from admission to that culture is
#' \eqn{\geq} \code{hai_threshold_hours} (default 48 h), the patient is
#' classified as HAI; otherwise CAI. The existing \code{type_of_infection}
#' column in the dataset is \emph{ignored} -- classification is derived
#' entirely from dates.
#'
#' \strong{Deduplication:} one classification per patient per centre
#' (earliest culture date).
#'
#' @param data                  Data frame. Long-format AMR dataset.
#' @param patient_col           Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param center_col            Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param admission_col         Character. Admission date column.
#'   Default \code{"admission_date"}.
#' @param culture_col           Character. Culture/specimen date column.
#'   Default \code{"culture_date"}.
#' @param hai_threshold_hours   Numeric. Hours from admission to first
#'   culture at or above which a case is classified as HAI. Default \code{48}.
#' @param hai_label             Character. Label for HAI bars.
#'   Default \code{"HAI"}.
#' @param cai_label             Character. Label for CAI bars.
#'   Default \code{"CAI"}.
#' @param colours               Named character vector. Fill colours keyed
#'   by \code{hai_label} and \code{cai_label}. Default red/green.
#' @param base_size             Numeric. Base font size. Default \code{14}.
#' @param title                 Character. Custom title. Auto-generated if
#'   \code{NULL}.
#' @param syndrome_col          Character or \code{NULL}. Column name
#'   containing syndrome/infection category labels. \code{NULL} = no filter.
#' @param syndrome_name         Character or \code{NULL}. Syndrome value to
#'   retain. Requires \code{syndrome_col}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_hai_cai_by_facility <- function(data,
                                     patient_col         = "PatientInformation_id",
                                     center_col          = "center_name",
                                     admission_col       = "admission_date",
                                     culture_col         = "culture_date",
                                     hai_threshold_hours = 48,
                                     hai_label           = "HAI",
                                     cai_label           = "CAI",
                                     colours             = NULL,
                                     base_size           = 14,
                                     title               = NULL,
                                     syndrome_col        = NULL,
                                     syndrome_name       = NULL) {

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, center_col, admission_col, culture_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s",
                 paste(missing_cols, collapse = ", ")))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  ctr_sym <- rlang::sym(center_col)
  adm_sym <- rlang::sym(admission_col)
  cul_sym <- rlang::sym(culture_col)

  # -- 3. classify HAI / CAI -- one row per patient x centre -----------------
  classified <- data %>%
    dplyr::filter(!is.na(!!adm_sym), !is.na(!!cul_sym)) %>%
    dplyr::group_by(!!pt_sym, !!ctr_sym) %>%
    dplyr::slice_min(!!cul_sym, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      infection_type_derived = factor(
        dplyr::if_else(
          as.numeric(difftime(!!cul_sym, !!adm_sym, units = "hours")) >=
            hai_threshold_hours,
          hai_label, cai_label
        ),
        levels = c(cai_label, hai_label)
      )
    )

  # -- 4. count unique patients per centre x infection type -----------------
  plot_data <- classified %>%
    dplyr::count(!!ctr_sym, infection_type_derived, name = "n_patients") %>%
    dplyr::group_by(!!ctr_sym) %>%
    dplyr::mutate(
      total      = sum(n_patients),
      proportion = n_patients / total,
      label      = paste0(n_patients, "\n(", round(100 * proportion), "%)")
    ) %>%
    dplyr::ungroup()

  n_labels <- plot_data %>%
    dplyr::distinct(!!ctr_sym, total)

  # -- 5. resolve colours ----------------------------------------------------
  if (is.null(colours)) {
    colours <- stats::setNames(
      c("#2ECC71", "#E74C3C"),
      c(cai_label, hai_label)
    )
  }

  # -- 6. build plot ---------------------------------------------------------
  auto_title <- title %||% sprintf(
    "HAI vs CAI by Facility (threshold \u2265 %d h from admission)",
    hai_threshold_hours
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = !!ctr_sym, y = n_patients, fill = infection_type_derived)
  ) +
    ggplot2::geom_col(color = "black", width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      size     = 3.2,
      color    = "white",
      fontface = "bold",
      lineheight = 0.9
    ) +
    ggplot2::geom_text(
      data = n_labels,
      ggplot2::aes(x = !!ctr_sym, y = total,
                   label = paste0("n=", total)),
      inherit.aes = FALSE,
      vjust    = -0.45,
      size     = 3.4,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.14))
    ) +
    ggplot2::labs(
      x     = "Facility",
      y     = "Number of Unique Patients",
      fill  = "Infection Type",
      title = auto_title
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1,
                                                  face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )

  return(p)
}



# PATIENT DISTRIBUTION BY LOCATION TYPE ACROSS FACILITIES


#' Plot Patient Distribution by Location Type Across Facilities
#'
#' Produces a dodged bar chart showing the number of unique patients in each
#' location type (ICU, Ward, Other) per facility. Raw location strings are
#' normalised automatically via case-insensitive pattern matching, so
#' centre-specific ICU names (e.g. \code{"Respiratory ICU (RICU)"}) all map
#' to \code{"ICU"}.
#'
#' \strong{Normalisation:} a location value is classified as:
#' \itemize{
#'   \item \strong{ICU} -- if the lower-cased string contains
#'     \code{icu_pattern} (default \code{"icu"})
#'   \item \strong{Ward} -- else if it contains \code{ward_pattern}
#'     (default \code{"ward"})
#'   \item \strong{Other} -- everything else
#' }
#'
#' \strong{Deduplication:} \code{distinct(patient, centre, location type)} --
#' each unique patient-centre-location combination is counted once.
#'
#' @param data           Data frame. Long-format AMR dataset.
#' @param patient_col    Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param location_col   Character. Location/unit column.
#'   Default \code{"location"}.
#' @param center_col     Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param icu_pattern    Character. Case-insensitive regex pattern that
#'   identifies ICU values. Default \code{"icu"}.
#' @param ward_pattern   Character. Case-insensitive regex pattern that
#'   identifies Ward values. Default \code{"ward"}.
#' @param icu_label      Character. Display label for ICU bars.
#'   Default \code{"ICU"}.
#' @param ward_label     Character. Display label for Ward bars.
#'   Default \code{"Ward"}.
#' @param other_label    Character. Display label for all other locations.
#'   Default \code{"Other"}.
#' @param colours        Named character vector. Fill colours keyed by
#'   \code{icu_label}, \code{ward_label}, and \code{other_label}.
#'   Default pink/green/grey.
#' @param bar_width      Numeric. Width of each individual bar.
#'   Default \code{0.68}.
#' @param base_size      Numeric. Base font size. Default \code{14}.
#' @param title          Character. Custom title. Auto-generated if
#'   \code{NULL}.
#' @param syndrome_col   Character or \code{NULL}. Column name containing
#'   syndrome/infection category labels. \code{NULL} = no filter.
#' @param syndrome_name  Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_location_by_facility <- function(data,
                                      patient_col   = "PatientInformation_id",
                                      location_col  = "location",
                                      center_col    = "center_name",
                                      icu_pattern   = "icu",
                                      ward_pattern  = "ward",
                                      icu_label     = "ICU",
                                      ward_label    = "Ward",
                                      other_label   = "Other",
                                      colours       = NULL,
                                      bar_width     = 0.68,
                                      base_size     = 14,
                                      title         = NULL,
                                      syndrome_col  = NULL,
                                      syndrome_name = NULL) {

  # -- 1. validate columns ---------------------------------------------------
  required_cols <- c(patient_col, location_col, center_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s",
                 paste(missing_cols, collapse = ", ")))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter ---------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols --------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  loc_sym <- rlang::sym(location_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. normalise location -> ICU / Ward / Other ---------------------------
  clean <- data %>%
    dplyr::filter(!is.na(!!loc_sym)) %>%
    dplyr::mutate(
      location_type = dplyr::case_when(
        grepl(icu_pattern,  tolower(as.character(!!loc_sym))) ~ icu_label,
        grepl(ward_pattern, tolower(as.character(!!loc_sym))) ~ ward_label,
        TRUE                                                  ~ other_label
      ),
      location_type = factor(location_type,
                             levels = c(icu_label, ward_label, other_label))
    )

  # -- 4. deduplicate: one row per patient x centre x location type ----------
  plot_data <- clean %>%
    dplyr::distinct(!!pt_sym, !!ctr_sym, location_type) %>%
    dplyr::count(!!ctr_sym, location_type, name = "n_patients")

  # -- 5. resolve colours ----------------------------------------------------
  if (is.null(colours)) {
    colours <- stats::setNames(
      c("#E91E63", "#4CAF50", "#9E9E9E"),
      c(icu_label, ward_label, other_label)
    )
  }

  # -- 6. build plot ---------------------------------------------------------
  auto_title <- title %||% "Patient Distribution by Location Type across Facilities"

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = !!ctr_sym, y = n_patients, fill = location_type)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = bar_width + 0.04),
      width    = bar_width,
      alpha    = 0.88
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_patients),
      position = ggplot2::position_dodge(width = bar_width + 0.04),
      vjust    = -0.45,
      size     = 3.4,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      x     = "Facility",
      y     = "Number of Unique Patients",
      fill  = "Location",
      title = auto_title
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1,
                                                  face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )

  return(p)
}



# LOS RIDGE PLOT


#' Ridge / Density Plot of Length of Stay (LOS)
#'
#' Visualises the distribution of length of stay (LOS) on a log10 x-axis using
#' either a stacked ridge plot (\code{mode = "all"}) or a single-centre density
#' curve (\code{mode = "single"}).
#'
#' LOS is computed as
#' \code{difftime(discharge_col, admission_col, units = "days")}.
#' One row per unique patient x centre x admission_date x discharge_date is
#' retained before plotting. Rows with \code{LOS < min_los} or
#' \code{LOS > max_los} are dropped.
#'
#' @param data A data frame.
#' @param mode Character. One of \code{"all"} (ridge plot, all centres stacked
#'   on y-axis, ordered by median LOS ascending) or \code{"single"} (density
#'   curve for one centre specified by \code{center}).
#' @param center Character. Centre name required when \code{mode = "single"}.
#' @param patient_col Character. Column of patient identifiers.
#' @param center_col Character. Column of centre identifiers.
#' @param admission_col Character. Column containing admission date/datetime.
#'   Default \code{"admission_date"}.
#' @param discharge_col Character. Column containing discharge / final outcome
#'   date. Default \code{"final_outcome_date"}.
#' @param filter_outcome Character or \code{NULL}. If not \code{NULL}, only
#'   rows whose outcome column (\code{outcome_col}) equals this value are used
#'   (e.g. \code{"Discharged"}). Set \code{NULL} to include all outcomes.
#'   Default \code{NULL}.
#' @param outcome_col Character. Column of patient outcomes, used only when
#'   \code{filter_outcome} is not \code{NULL}. Default \code{"final_outcome"}.
#' @param min_los Numeric. Minimum LOS (days) to include. Rows below this
#'   are removed (prevents \code{log(0)} issues). Default \code{1}.
#' @param max_los Numeric. Maximum LOS (days) to include. Default \code{365}.
#' @param colours Named character vector of colours for centres. If \code{NULL}
#'   (default), the \pkg{RColorBrewer} "Set2" palette is used.
#' @param scale Numeric. \code{scale} argument passed to
#'   \code{ggridges::geom_density_ridges()}. Controls overlap between ridges.
#'   Default \code{1.2}.
#' @param alpha Numeric. Fill transparency (0-1). Default \code{0.7}.
#' @param base_size Numeric. Base font size. Default \code{14}.
#' @param title Character or \code{NULL}. Plot title. Auto-generated if
#'   \code{NULL}.
#' @param syndrome_col Character or \code{NULL}. Column of syndrome labels for
#'   pre-filtering. Default \code{NULL}.
#' @param syndrome_name Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}. Default \code{NULL}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_los_ridge <- function(
    data,
    mode              = c("all", "single"),
    center            = NULL,
    patient_col       = "patient_id",
    center_col        = "center_name",
    admission_col     = "admission_date",
    discharge_col     = "final_outcome_date",
    filter_outcome    = NULL,
    outcome_col       = "final_outcome",
    min_los           = 1,
    max_los           = 365,
    colours           = NULL,
    scale             = 1.2,
    alpha             = 0.7,
    base_size         = 14,
    title             = NULL,
    syndrome_col      = NULL,
    syndrome_name     = NULL
) {

  # -- step 0: validate arguments -----------------------------------------------
  mode <- match.arg(mode)
  if (mode == "single" && is.null(center))
    stop("'center' must be provided when mode = 'single'.")

  # -- step 1: column validation ------------------------------------------------
  required_cols <- c(patient_col, center_col, admission_col, discharge_col)
  if (!is.null(filter_outcome)) required_cols <- c(required_cols, outcome_col)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0)
    stop("The following required columns are missing from data: ",
         paste(missing_cols, collapse = ", "))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter -------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- step 2: tidy-eval symbols -------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  ctr_sym <- rlang::sym(center_col)
  adm_sym <- rlang::sym(admission_col)
  dis_sym <- rlang::sym(discharge_col)

  # -- step 3: optional outcome filter ------------------------------------------
  if (!is.null(filter_outcome)) {
    out_sym <- rlang::sym(outcome_col)
    data <- data %>%
      dplyr::filter(!is.na(!!out_sym), !!out_sym == filter_outcome)
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", outcome_col, filter_outcome))
  }

  # -- step 4: compute LOS -------------------------------------------------------
  clean <- data %>%
    dplyr::filter(!is.na(!!adm_sym), !is.na(!!dis_sym)) %>%
    dplyr::distinct(!!pt_sym, !!ctr_sym, !!adm_sym, !!dis_sym) %>%
    dplyr::mutate(
      LOS_days = as.numeric(difftime(!!dis_sym, !!adm_sym, units = "days"))
    ) %>%
    dplyr::filter(LOS_days >= min_los, LOS_days <= max_los)

  if (nrow(clean) == 0)
    stop("No valid LOS values remain after filtering.")

  # -- step 5: mode-specific prep ------------------------------------------------
  x_breaks <- c(1, 7, 14, 30, 60, 90, 180, 365)
  x_labels <- c("1", "7", "14", "30", "60", "90", "180", "365")

  if (mode == "all") {

    # Order centres by median LOS (ascending)
    center_order <- clean %>%
      dplyr::group_by(!!ctr_sym) %>%
      dplyr::summarise(med_los = stats::median(LOS_days), .groups = "drop") %>%
      dplyr::arrange(med_los) %>%
      dplyr::pull(!!ctr_sym)

    clean <- clean %>%
      dplyr::mutate(
        !!center_col := factor(!!ctr_sym, levels = center_order)
      )

    n_centers <- length(center_order)
    if (is.null(colours)) {
      colours <- stats::setNames(
        grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(max(min(n_centers, 8), 3), "Set2")
        )(n_centers),
        center_order
      )
    }

    auto_title <- title %||%
      sprintf("Distribution of Length of Stay by Centre (n = %s patients)",
              scales::comma(dplyr::n_distinct(clean[[patient_col]])))

    p <- ggplot2::ggplot(
      clean,
      ggplot2::aes(
        x    = LOS_days,
        y    = !!ctr_sym,
        fill = !!ctr_sym
      )
    ) +
      ggridges::geom_density_ridges(
        scale            = scale,
        alpha            = alpha,
        rel_min_height   = 0.01,
        quantile_lines   = TRUE,
        quantiles        = 2
      ) +
      ggplot2::scale_x_log10(
        breaks = x_breaks,
        labels = x_labels,
        limits = c(min_los, max_los)
      ) +
      ggplot2::scale_fill_manual(values = colours, guide = "none") +
      ggplot2::labs(
        x     = "Length of Stay (days, log scale)",
        y     = "Centre",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none") +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(
          colour = "grey85", linetype = "dashed"
        )
      )

  } else {  # mode == "single"

    clean_single <- clean %>%
      dplyr::filter(!!ctr_sym == center)

    if (nrow(clean_single) == 0)
      stop(sprintf("No data found for centre '%s'.", center))

    auto_title <- title %||%
      sprintf("Length of Stay Distribution -- %s (n = %s patients)",
              center,
              scales::comma(dplyr::n_distinct(clean_single[[patient_col]])))

    fill_col <- if (!is.null(colours) && center %in% names(colours)) {
      colours[[center]]
    } else {
      "#4DAF4A"
    }

    p <- ggplot2::ggplot(
      clean_single,
      ggplot2::aes(x = LOS_days)
    ) +
      ggplot2::geom_density(
        fill  = fill_col,
        colour = "grey30",
        alpha = alpha
      ) +
      ggplot2::scale_x_log10(
        breaks = x_breaks,
        labels = x_labels,
        limits = c(min_los, max_los)
      ) +
      ggplot2::labs(
        x     = "Length of Stay (days, log scale)",
        y     = "Density",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none") +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(
          colour = "grey85", linetype = "dashed"
        )
      )
  }

  return(p)
}



# AGE RIDGE PLOT


#' Ridge / Density Plot of Patient Age
#'
#' Visualises the distribution of patient age on the x-axis using either a
#' stacked ridge plot (\code{mode = "all"}) or a single-centre density curve
#' (\code{mode = "single"}).
#'
#' One row per unique patient x centre is retained before plotting.
#' Rows with age outside \code{[min_age, max_age]} are dropped.
#'
#' @param data A data frame.
#' @param mode Character. One of \code{"all"} (ridge plot, all centres stacked
#'   on y-axis ordered by median age ascending) or \code{"single"} (density
#'   curve for one centre specified by \code{center}).
#' @param center Character. Centre name required when \code{mode = "single"}.
#' @param patient_col Character. Column of patient identifiers.
#' @param age_col Character. Column of patient age (numeric, in years).
#'   Default \code{"age_years"}.
#' @param center_col Character. Column of centre identifiers.
#' @param min_age Numeric. Minimum age to include. Default \code{0}.
#' @param max_age Numeric. Maximum age to include. Default \code{120}.
#' @param colours Named character vector of colours for centres. If \code{NULL}
#'   (default), the \pkg{RColorBrewer} "Set2" palette is used.
#' @param scale Numeric. \code{scale} argument passed to
#'   \code{ggridges::geom_density_ridges()}. Controls overlap between ridges.
#'   Default \code{1.2}.
#' @param alpha Numeric. Fill transparency (0-1). Default \code{0.7}.
#' @param base_size Numeric. Base font size. Default \code{14}.
#' @param title Character or \code{NULL}. Plot title. Auto-generated if
#'   \code{NULL}.
#' @param syndrome_col Character or \code{NULL}. Column of syndrome labels for
#'   pre-filtering. Default \code{NULL}.
#' @param syndrome_name Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}. Default \code{NULL}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_age_ridge <- function(
    data,
    mode          = c("all", "single"),
    center        = NULL,
    patient_col   = "PatientInformation_id",
    age_col       = "age_years",
    center_col    = "center_name",
    min_age       = 0,
    max_age       = 120,
    colours       = NULL,
    scale         = 1.2,
    alpha         = 0.7,
    base_size     = 14,
    title         = NULL,
    syndrome_col  = NULL,
    syndrome_name = NULL
) {

  # -- step 0: validate arguments -----------------------------------------------
  mode <- match.arg(mode)
  if (mode == "single" && is.null(center))
    stop("'center' must be provided when mode = 'single'.")

  # -- step 1: column validation ------------------------------------------------
  required_cols <- c(patient_col, age_col, center_col)
  missing_cols  <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0)
    stop("The following required columns are missing from data: ",
         paste(missing_cols, collapse = ", "))

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter -------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- step 2: tidy-eval symbols -------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  age_sym <- rlang::sym(age_col)
  ctr_sym <- rlang::sym(center_col)

  # -- step 3: deduplicate and filter -------------------------------------------
  clean <- data %>%
    dplyr::filter(!is.na(!!age_sym)) %>%
    dplyr::mutate(!!age_col := suppressWarnings(as.numeric(!!age_sym))) %>%
    dplyr::filter(!is.na(!!age_sym)) %>%
    dplyr::distinct(!!pt_sym, !!ctr_sym, .keep_all = TRUE) %>%
    dplyr::filter(!!age_sym >= min_age, !!age_sym <= max_age)

  if (nrow(clean) == 0)
    stop("No valid age values remain after filtering.")

  # -- step 4: mode-specific prep ------------------------------------------------
  if (mode == "all") {

    # order centres by median age ascending
    center_order <- clean %>%
      dplyr::group_by(!!ctr_sym) %>%
      dplyr::summarise(med_age = stats::median(!!age_sym), .groups = "drop") %>%
      dplyr::arrange(med_age) %>%
      dplyr::pull(!!ctr_sym)

    clean <- clean %>%
      dplyr::mutate(!!center_col := factor(!!ctr_sym, levels = center_order))

    n_centers <- length(center_order)
    if (is.null(colours)) {
      colours <- stats::setNames(
        grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(max(min(n_centers, 8), 3), "Set2")
        )(n_centers),
        center_order
      )
    }

    auto_title <- title %||%
      sprintf("Age Distribution by Centre (n = %s patients)",
              scales::comma(dplyr::n_distinct(clean[[patient_col]])))

    p <- ggplot2::ggplot(
      clean,
      ggplot2::aes(
        x    = !!age_sym,
        y    = !!ctr_sym,
        fill = !!ctr_sym
      )
    ) +
      ggridges::geom_density_ridges(
        scale          = scale,
        alpha          = alpha,
        rel_min_height = 0.01,
        quantile_lines = TRUE,
        quantiles      = 2
      ) +
      ggplot2::scale_x_continuous(
        limits = c(min_age, max_age),
        breaks = scales::pretty_breaks(n = 8)
      ) +
      ggplot2::scale_fill_manual(values = colours, guide = "none") +
      ggplot2::labs(
        x     = "Age (years)",
        y     = "Centre",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none") +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(
          colour = "grey85", linetype = "dashed"
        )
      )

  } else {  # mode == "single"

    clean_single <- clean %>%
      dplyr::filter(!!ctr_sym == center)

    if (nrow(clean_single) == 0)
      stop(sprintf("No data found for centre '%s'.", center))

    auto_title <- title %||%
      sprintf("Age Distribution -- %s (n = %s patients)",
              center,
              scales::comma(dplyr::n_distinct(clean_single[[patient_col]])))

    fill_col <- if (!is.null(colours) && center %in% names(colours)) {
      colours[[center]]
    } else {
      "#4DAF4A"
    }

    p <- ggplot2::ggplot(
      clean_single,
      ggplot2::aes(x = !!age_sym)
    ) +
      ggplot2::geom_density(
        fill   = fill_col,
        colour = "grey30",
        alpha  = alpha
      ) +
      ggplot2::scale_x_continuous(
        limits = c(min_age, max_age),
        breaks = scales::pretty_breaks(n = 8)
      ) +
      ggplot2::labs(
        x     = "Age (years)",
        y     = "Density",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "none") +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(
          colour = "grey85", linetype = "dashed"
        )
      )
  }

  return(p)
}



# LOS BY AGE GROUP BOXPLOT


#' Boxplot of Length of Stay by Age Group
#'
#' Displays the distribution of LOS (y-axis, log10 scale) across patient age
#' groups (x-axis) using boxplots. Supports \code{"faceted"} (one panel per
#' centre), \code{"overall"} (all centres pooled), and \code{"single"}
#' (one centre) modes.
#'
#' LOS is computed as \code{difftime(discharge_col, admission_col, units = "days")}.
#' One row per unique patient x centre x admission x discharge is retained.
#' Rows with \code{LOS < min_los} or \code{LOS > max_los} are dropped.
#'
#' \strong{Two ways to supply age groups:}
#' \enumerate{
#'   \item \strong{Pre-existing column} (default) -- set \code{agebin_col} to an
#'     existing categorical column (e.g. \code{"Age_bin"}). Bins are sorted
#'     automatically (\code{<1}, \code{1-5}, ..., \code{85+}); override order
#'     with \code{age_levels}.
#'   \item \strong{Custom breaks from continuous age} -- set \code{age_breaks}
#'     to a numeric vector (e.g. \code{c(0, 5, 15, 30, 50, 70, Inf)}) and
#'     \code{age_col} to the continuous age column (default \code{"age_years"}).
#'     Bins are derived via \code{cut()} and labelled automatically.
#'     \code{agebin_col} is ignored when \code{age_breaks} is provided.
#' }
#'
#' @param data A data frame.
#' @param mode Character. One of \code{"faceted"}, \code{"overall"}, or
#'   \code{"single"}.
#' @param center Character. Centre name required when \code{mode = "single"}.
#' @param patient_col Character. Column of patient identifiers.
#' @param agebin_col Character. Pre-existing age group column. Used when
#'   \code{age_breaks} is \code{NULL}. Default \code{"Age_bin"}.
#' @param age_col Character. Continuous age column used when \code{age_breaks}
#'   is provided. Default \code{"age_years"}.
#' @param age_breaks Numeric vector or \code{NULL}. Breakpoints for deriving
#'   custom bins from \code{age_col} (e.g. \code{c(0, 5, 15, 30, 50, 70, Inf)}).
#'   \code{NULL} = use \code{agebin_col} directly.
#' @param age_labels Character vector or \code{NULL}. Labels for each bin
#'   created by \code{age_breaks}. Must have length
#'   \code{length(age_breaks) - 1}. \code{NULL} = auto-generate from breaks.
#' @param center_col Character. Column of centre identifiers.
#' @param admission_col Character. Admission date column. Default
#'   \code{"admission_date"}.
#' @param discharge_col Character. Discharge / final outcome date column.
#'   Default \code{"final_outcome_date"}.
#' @param min_los Numeric. Minimum LOS (days) to include. Default \code{1}.
#' @param max_los Numeric. Maximum LOS (days) to include. Default \code{365}.
#' @param age_levels Character vector or \code{NULL}. Custom bin order when
#'   using \code{agebin_col}. \code{NULL} = auto-sort.
#' @param fill_colour Character. Box fill colour. Default \code{"#2C7FB8"}.
#' @param fill_alpha Numeric. Box fill transparency (0-1). Default \code{0.7}.
#' @param outlier_alpha Numeric. Outlier point transparency (0-1). Default
#'   \code{0.25}.
#' @param box_width Numeric. Width of each box. Default \code{0.6}.
#' @param ncol Integer. Number of facet columns (\code{mode = "faceted"} only).
#'   Default \code{3}.
#' @param base_size Numeric. Base font size. Default \code{14}.
#' @param title Character or \code{NULL}. Plot title. Auto-generated if
#'   \code{NULL}.
#' @param subtitle Character or \code{NULL}. Plot subtitle. \code{NULL} = none.
#' @param syndrome_col Character or \code{NULL}. Syndrome column for
#'   pre-filtering. Default \code{NULL}.
#' @param syndrome_name Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}. Default \code{NULL}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_los_by_agebin <- function(
    data,
    mode           = c("faceted", "overall", "single"),
    center         = NULL,
    patient_col    = "PatientInformation_id",
    agebin_col     = "Age_bin",
    age_col        = "age_years",
    age_breaks     = NULL,
    age_labels     = NULL,
    center_col     = "center_name",
    admission_col  = "admission_date",
    discharge_col  = "final_outcome_date",
    min_los        = 1,
    max_los        = 365,
    age_levels     = NULL,
    fill_colour    = "#2C7FB8",
    fill_alpha     = 0.7,
    outlier_alpha  = 0.25,
    box_width      = 0.6,
    ncol           = 3,
    base_size      = 14,
    title          = NULL,
    subtitle       = NULL,
    syndrome_col   = NULL,
    syndrome_name  = NULL
) {

  # -- step 0: validate arguments -----------------------------------------------
  mode <- match.arg(mode)
  if (mode == "single" && is.null(center))
    stop("'center' must be provided when mode = 'single'.")

  # -- step 1: column validation ------------------------------------------------
  # Determine which age column to validate against
  age_source <- if (!is.null(age_breaks)) age_col else agebin_col
  required_cols <- c(patient_col, age_source, center_col, admission_col, discharge_col)
  missing_cols  <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0)
    stop("The following required columns are missing from data: ",
         paste(missing_cols, collapse = ", "))

  if (!is.null(age_breaks) && !is.null(age_labels) &&
      length(age_labels) != length(age_breaks) - 1)
    stop("'age_labels' must have length(age_breaks) - 1 elements.")

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter -------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- step 2: tidy-eval symbols -------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  ctr_sym <- rlang::sym(center_col)
  adm_sym <- rlang::sym(admission_col)
  dis_sym <- rlang::sym(discharge_col)

  # -- step 3: derive age bins if age_breaks supplied ---------------------------
  bin_col <- ".__age_bin__"   # internal working column name

  if (!is.null(age_breaks)) {
    age_sym <- rlang::sym(age_col)

    # auto-generate labels like "<5", "5-14", "15-29", ... if not supplied
    if (is.null(age_labels)) {
      n_breaks <- length(age_breaks)
      age_labels <- vapply(seq_len(n_breaks - 1), function(i) {
        lo <- age_breaks[i]
        hi <- age_breaks[i + 1]
        if (lo == 0 || lo == -Inf)
          sprintf("<%g", hi)
        else if (hi == Inf)
          sprintf("%g+", lo)
        else
          sprintf("%g\u2013%g", lo, hi - 1)   # e.g. "5-14"
      }, character(1))
    }

    data <- data %>%
      dplyr::mutate(
        !!bin_col := cut(
          suppressWarnings(as.numeric(!!age_sym)),
          breaks         = age_breaks,
          labels         = age_labels,
          right          = FALSE,
          include.lowest = TRUE
        )
      )
    lvls <- age_labels

  } else {
    # use pre-existing agebin_col, copy to working column
    ab_sym <- rlang::sym(agebin_col)
    data <- data %>% dplyr::mutate(!!bin_col := !!ab_sym)

    # auto-sort if no age_levels given
    if (!is.null(age_levels)) {
      lvls <- age_levels
    } else {
      all_bins   <- unique(data[[bin_col]])
      all_bins   <- all_bins[!is.na(all_bins)]
      lt_bins    <- all_bins[grepl("^<", all_bins)]
      plus_bins  <- all_bins[grepl("\\+$", all_bins)]
      range_bins <- setdiff(all_bins, c(lt_bins, plus_bins))
      lower_num  <- suppressWarnings(as.numeric(sub("[-\u2013].*", "", range_bins)))
      range_bins <- range_bins[order(lower_num)]
      lvls <- c(sort(lt_bins), range_bins, sort(plus_bins))
    }
  }

  bin_sym <- rlang::sym(bin_col)

  # -- step 4: compute LOS, deduplicate, filter ---------------------------------
  clean <- data %>%
    dplyr::filter(!is.na(!!adm_sym), !is.na(!!dis_sym), !is.na(!!bin_sym)) %>%
    dplyr::distinct(!!pt_sym, !!ctr_sym, !!bin_sym, !!adm_sym, !!dis_sym) %>%
    dplyr::mutate(
      LOS_days    = as.numeric(difftime(!!dis_sym, !!adm_sym, units = "days")),
      !!bin_col  := factor(!!bin_sym, levels = lvls)
    ) %>%
    dplyr::filter(LOS_days >= min_los, LOS_days <= max_los)

  if (nrow(clean) == 0)
    stop("No valid LOS values remain after filtering.")

  # -- step 5: mode filter -------------------------------------------------------
  if (mode == "single") {
    clean <- clean %>% dplyr::filter(!!ctr_sym == center)
    if (nrow(clean) == 0)
      stop(sprintf("No data found for centre '%s'.", center))
  }

  # -- step 6: y-axis breaks (include 3 to match original) ---------------------
  y_breaks <- c(1, 3, 7, 14, 30, 60, 90, 180, 365)
  y_labels <- as.character(y_breaks)

  # -- step 7: build plot -------------------------------------------------------
  n_pts <- dplyr::n_distinct(clean[[patient_col]])

  auto_title <- title %||% switch(
    mode,
    "overall" = sprintf("LOS by Age Group \u2014 All Centres (n = %s)", scales::comma(n_pts)),
    "single"  = sprintf("LOS by Age Group \u2014 %s (n = %s)", center, scales::comma(n_pts)),
    sprintf("LOS by Age Group (n = %s)", scales::comma(n_pts))
  )

  p <- ggplot2::ggplot(
    clean,
    ggplot2::aes(x = !!bin_sym, y = LOS_days)
  ) +
    ggplot2::geom_boxplot(
      fill          = fill_colour,
      alpha         = fill_alpha,
      colour        = "grey30",
      width         = box_width,
      outlier.size  = 0.8,
      outlier.alpha = outlier_alpha
    ) +
    ggplot2::scale_y_log10(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(min_los, max_los)
    ) +
    ggplot2::labs(
      x        = "Age Group",
      y        = "LOS (days, log scale)",
      title    = auto_title,
      subtitle = subtitle
    ) +
    eda_theme(base_size = base_size) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 20, hjust = 1),
      panel.grid.major.y = ggplot2::element_line(colour = "grey85",
                                                  linetype = "dashed"),
      panel.grid.major.x = ggplot2::element_blank()
    )

  if (mode == "faceted") {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", center_col)),
      ncol   = ncol,
      scales = "fixed"
    )
  }

  return(p)
}



# FINAL OUTCOME BY YEAR


#' Plot Distribution of Final Outcomes by Year
#'
#' Produces a stacked bar chart where the x-axis shows the year derived from
#' \code{date_col}, the y-axis shows the number of unique patients, and bars
#' are filled by final outcome category. Each bar segment is labelled with
#' count and percentage; a total \code{n=} label sits above each bar.
#'
#' One outcome is taken per patient per centre (the first recorded value) to
#' avoid double-counting long-format rows.
#'
#' Supports three display modes:
#' \itemize{
#'   \item \strong{"faceted"} (default) -- one panel per centre.
#'   \item \strong{"overall"} -- all centres pooled into a single chart.
#'   \item \strong{"single"} -- one specific centre; pass its name via
#'     \code{center}.
#' }
#'
#' @param data           Data frame. Long-format AMR dataset.
#' @param mode           Character. One of \code{"faceted"}, \code{"overall"},
#'   or \code{"single"}. Default \code{"faceted"}.
#' @param center         Character. Required when \code{mode = "single"}.
#'   Exact centre name.
#' @param patient_col    Character. Patient ID column.
#'   Default \code{"PatientInformation_id"}.
#' @param outcome_col    Character. Final outcome column.
#'   Default \code{"final_outcome"}.
#' @param date_col       Character. Date column used to extract the year
#'   (must be a \code{Date} or \code{POSIXct} column).
#'   Default \code{"final_outcome_date"}.
#' @param center_col     Character. Centre/facility column.
#'   Default \code{"center_name"}.
#' @param merge_referred Logical. Recode
#'   \code{"Transferred to other hospital"} to \code{"Referred"}.
#'   Default \code{TRUE}.
#' @param palette        Named character vector mapping outcome values to
#'   colours. Unmatched outcomes get a default ggplot colour.
#' @param ncol           Integer. Facet columns (\code{mode = "faceted"} only).
#'   Default 2.
#' @param base_size      Numeric. Base font size. Default 14.
#' @param title          Character. Custom plot title. Auto-generated if
#'   \code{NULL}. Default \code{NULL}.
#' @param syndrome_col   Character or \code{NULL}. Syndrome filter column.
#'   Default \code{NULL}.
#' @param syndrome_name  Character or \code{NULL}. Syndrome value to retain.
#'   Requires \code{syndrome_col}. Default \code{NULL}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
plot_outcome_by_year <- function(data,
                                 mode           = c("faceted", "overall", "single"),
                                 center         = NULL,
                                 patient_col    = "PatientInformation_id",
                                 outcome_col    = "final_outcome",
                                 date_col       = "final_outcome_date",
                                 center_col     = "center_name",
                                 merge_referred = TRUE,
                                 palette        = c(
                                   "Death"                         = "#E74C3C",
                                   "Died"                          = "#E74C3C",
                                   "Discharged"                    = "#2ECC71",
                                   "LAMA"                          = "#95A5A6",
                                   "Referred"                      = "#3498DB",
                                   "Transferred to other hospital" = "#3498DB"
                                 ),
                                 ncol           = 2,
                                 base_size      = 14,
                                 title          = NULL,
                                 syndrome_col   = NULL,
                                 syndrome_name  = NULL) {

  # -- 0. match mode -----------------------------------------------------------
  mode <- match.arg(mode)

  # -- 1. validate columns -----------------------------------------------------
  required_cols <- c(patient_col, outcome_col, date_col)
  if (mode != "overall") required_cols <- c(required_cols, center_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Column(s) not found in data: %s",
                 paste(missing_cols, collapse = ", ")))

  if (mode == "single") {
    if (is.null(center))
      stop("'center' must be provided when mode = 'single'.")
    available <- unique(data[[center_col]])
    if (!center %in% available)
      stop(sprintf("'%s' not found in column '%s'. Available values: %s",
                   center, center_col,
                   paste(sort(available), collapse = ", ")))
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("'syndrome_col' must be provided when 'syndrome_name' is set.")
  if (!is.null(syndrome_col) && !syndrome_col %in% names(data))
    stop(sprintf("syndrome_col '%s' not found in data.", syndrome_col))

  # -- syndrome pre-filter -----------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    data <- data[!is.na(data[[syndrome_col]]) &
                   data[[syndrome_col]] == syndrome_name, ]
    if (nrow(data) == 0)
      stop(sprintf("No rows found where %s == '%s'.", syndrome_col, syndrome_name))
  }

  # -- 2. tidy-eval symbols ----------------------------------------------------
  pt_sym  <- rlang::sym(patient_col)
  out_sym <- rlang::sym(outcome_col)
  dt_sym  <- rlang::sym(date_col)
  ctr_sym <- rlang::sym(center_col)

  # -- 3. clean: drop missing outcome or date ----------------------------------
  clean <- data %>%
    dplyr::filter(
      !is.na(!!out_sym), trimws(as.character(!!out_sym)) != "",
      !is.na(!!dt_sym)
    )

  # -- 4. optionally recode referred -------------------------------------------
  if (merge_referred) {
    clean <- clean %>%
      dplyr::mutate(
        !!outcome_col := dplyr::case_when(
          !!out_sym == "Transferred to other hospital" ~ "Referred",
          TRUE ~ as.character(!!out_sym)
        )
      )
    out_sym <- rlang::sym(outcome_col)
  }

  # -- 5. optional: filter to single centre ------------------------------------
  if (mode == "single")
    clean <- clean %>% dplyr::filter(!!ctr_sym == center)

  # -- 6. deduplicate: one outcome + date per patient (per centre) -------------
  if (mode == "overall") {
    deduped <- clean %>%
      dplyr::group_by(!!pt_sym) %>%
      dplyr::summarise(
        !!outcome_col := dplyr::first(!!out_sym),
        !!date_col    := dplyr::first(!!dt_sym),
        .groups = "drop"
      )
  } else {
    deduped <- clean %>%
      dplyr::group_by(!!ctr_sym, !!pt_sym) %>%
      dplyr::summarise(
        !!outcome_col := dplyr::first(!!out_sym),
        !!date_col    := dplyr::first(!!dt_sym),
        .groups = "drop"
      )
  }

  # -- 7. extract year ---------------------------------------------------------
  deduped <- deduped %>%
    dplyr::mutate(year = as.integer(format(!!dt_sym, "%Y"))) %>%
    dplyr::filter(!is.na(year))

  # -- 8. summarise: counts per year x outcome (x centre) ---------------------
  if (mode == "overall") {
    summary_yr <- deduped %>%
      dplyr::group_by(year, !!out_sym) %>%
      dplyr::summarise(unique_patients = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(year) %>%
      dplyr::mutate(
        total_year = sum(unique_patients),
        label      = paste0(unique_patients, " (",
                            round(100 * unique_patients / total_year, 1), "%)")
      ) %>%
      dplyr::ungroup()

    year_n <- summary_yr %>%
      dplyr::distinct(year, total_year)

  } else {
    summary_yr <- deduped %>%
      dplyr::group_by(!!ctr_sym, year, !!out_sym) %>%
      dplyr::summarise(unique_patients = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(!!ctr_sym, year) %>%
      dplyr::mutate(
        total_year = sum(unique_patients),
        label      = paste0(unique_patients, " (",
                            round(100 * unique_patients / total_year, 1), "%)")
      ) %>%
      dplyr::ungroup()

    year_n <- summary_yr %>%
      dplyr::distinct(!!ctr_sym, year, total_year)
  }

  # -- 9. build plot -----------------------------------------------------------
  if (mode == "overall") {
    auto_title <- title %||%
      "Distribution of Final Outcomes by Year \u2014 All Centres Pooled"

    p <- ggplot2::ggplot(
      summary_yr,
      ggplot2::aes(x = factor(year), y = unique_patients, fill = !!out_sym)
    ) +
      ggplot2::geom_col(color = "black", width = 0.75) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3.4, fontface = "bold", color = "white"
      ) +
      ggplot2::geom_text(
        data = year_n,
        ggplot2::aes(x = factor(year), y = total_year,
                     label = paste0("n=", total_year)),
        inherit.aes = FALSE,
        vjust = -0.4, size = 3.8, fontface = "bold"
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.12))
      ) +
      ggplot2::scale_fill_manual(values = palette, na.value = "grey70") +
      ggplot2::labs(
        x     = "Year",
        y     = "Number of Admissions",
        fill  = "Final Outcome",
        title = auto_title
      ) +
      eda_theme(base_size = base_size, legend_position = "top")

    return(p)
  }

  # faceted / single share the same base plot
  if (mode == "single") {
    auto_title <- title %||% sprintf(
      "Distribution of Final Outcomes by Year \u2014 %s", center
    )
  } else {
    auto_title <- title %||%
      "Distribution of Final Outcomes by Year \u2014 All Centres"
  }

  p <- ggplot2::ggplot(
    summary_yr,
    ggplot2::aes(x = factor(year), y = unique_patients, fill = !!out_sym)
  ) +
    ggplot2::geom_col(color = "black", width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3.4, fontface = "bold", color = "white"
    ) +
    ggplot2::geom_text(
      data = year_n,
      ggplot2::aes(x = factor(year), y = total_year,
                   label = paste0("n=", total_year)),
      inherit.aes = FALSE,
      vjust = -0.4, size = 3.8, fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::scale_fill_manual(values = palette, na.value = "grey70") +
    ggplot2::labs(
      x     = "Year",
      y     = "Number of Admissions",
      fill  = "Final Outcome",
      title = auto_title
    ) +
    eda_theme(base_size = base_size, legend_position = "top")

  if (mode == "faceted") {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", center_col)),
      ncol   = ncol,
      scales = "free_y"
    )
  }

  return(p)
}


# -- NULL coalescing operator (internal) --------------------------------------
# Used as: title %||% "default title"
`%||%` <- function(x, y) if (!is.null(x)) x else y
