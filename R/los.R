# los.R
# Length-of-stay distribution fitting and comparison

#' Safely Fit a Distribution
#'
#' Wrapper around \code{fitdistrplus::fitdist} that returns \code{NULL}
#' instead of throwing an error when fitting fails.
#'
#' @param x Numeric vector of positive values (e.g., LOS in days).
#' @param dist Character. Distribution name: \code{"weibull"}, \code{"lnorm"},
#'   or \code{"gamma"}.
#'
#' @return A \code{fitdist} object, or \code{NULL} on failure.
#' @export
#'
#' @examples
#' safe_fit(rlnorm(100), "lnorm")
safe_fit <- function(x, dist) {
  x <- x[is.finite(x)]
  if (length(x) < 5 || length(unique(x)) < 2) {
    return(NULL)
  }
  tryCatch(fitdistrplus::fitdist(x, dist), error = function(e) NULL)
}


#' Fit Multiple Distributions
#'
#' Fits Weibull, Lognormal, and Gamma distributions to a numeric vector,
#' returning all fit objects. Failed fits are \code{NULL}.
#'
#' @param x Numeric vector of positive values (e.g., LOS in days).
#'
#' @return Named list with elements \code{weibull}, \code{lnorm}, \code{gamma},
#'   each a \code{fitdist} object or \code{NULL}.
#' @export
#'
#' @examples
#' fits <- fit_distributions(rlnorm(200, 2, 0.5))
#' fits$lnorm$estimate
fit_distributions <- function(x) {
  list(
    weibull = safe_fit(x, "weibull"),
    lnorm   = safe_fit(x, "lnorm"),
    gamma   = safe_fit(x, "gamma")
  )
}


#' Compare Distribution Fits by AIC
#'
#' Fits Weibull, Lognormal, and Gamma distributions to a numeric vector and
#' returns their AIC values for comparison.
#'
#' @param x Numeric vector of positive values (e.g., LOS in days).
#' @param fits Optional. Pre-computed fits from \code{fit_distributions()}.
#'
#' @return Data frame with columns \code{Weibull_AIC}, \code{Lognormal_AIC},
#'   \code{Gamma_AIC}.
#' @export
#'
#' @examples
#' compare_distribution_aic(rlnorm(200, meanlog = 2, sdlog = 0.5))
compare_distribution_aic <- function(x, fits = NULL) {
  if (is.null(fits)) fits <- fit_distributions(x)

  data.frame(
    Weibull_AIC   = if (is.null(fits$weibull)) NA_real_ else fits$weibull$aic,
    Lognormal_AIC = if (is.null(fits$lnorm)) NA_real_ else fits$lnorm$aic,
    Gamma_AIC     = if (is.null(fits$gamma)) NA_real_ else fits$gamma$aic
  )
}


#' Summarise a Fitted Distribution
#'
#' Extracts mean, median, SD, and parameter values from a \code{fitdist}
#' object for Weibull, Lognormal, or Gamma distributions.
#'
#' @param fit A \code{fitdist} object from \code{fitdistrplus::fitdist()}.
#' @param dist Character. One of \code{"weibull"}, \code{"lnorm"}, or
#'   \code{"gamma"}.
#'
#' @return Data frame with columns \code{Mean_LOS}, \code{Median_LOS},
#'   \code{SD_LOS}, \code{Parameters}.
#' @export
#'
#' @examples
#' fit <- fitdistrplus::fitdist(rlnorm(200), "lnorm")
#' summarise_distribution(fit, "lnorm")
summarise_distribution <- function(fit, dist) {
  if (dist == "weibull") {
    k <- fit$estimate["shape"]
    l <- fit$estimate["scale"]
    mean_val <- l * gamma(1 + 1 / k)
    median_val <- l * (log(2))^(1 / k)
    sd_val <- sqrt(l^2 * (gamma(1 + 2 / k) - (gamma(1 + 1 / k))^2))
    params <- paste0("shape=", round(k, 3), ", scale=", round(l, 3))
  } else if (dist == "lnorm") {
    mu <- fit$estimate["meanlog"]
    sg <- fit$estimate["sdlog"]
    mean_val <- exp(mu + sg^2 / 2)
    median_val <- exp(mu)
    sd_val <- sqrt((exp(sg^2) - 1) * exp(2 * mu + sg^2))
    params <- paste0("meanlog=", round(mu, 3), ", sdlog=", round(sg, 3))
  } else if (dist == "gamma") {
    a <- fit$estimate["shape"]
    r <- fit$estimate["rate"]
    mean_val <- a / r
    median_val <- stats::qgamma(0.5, shape = a, rate = r)
    sd_val <- sqrt(a) / r
    params <- paste0("shape=", round(a, 3), ", rate=", round(r, 3))
  } else {
    stop("dist must be one of 'weibull', 'lnorm', 'gamma'")
  }

  data.frame(
    Mean_LOS   = round(mean_val, 2),
    Median_LOS = round(median_val, 2),
    SD_LOS     = round(sd_val, 2),
    Parameters = params
  )
}


#' Plot LOS Distribution with Fitted Overlays
#'
#' Creates a histogram of LOS values with Weibull, Lognormal, and Gamma
#' density curves overlaid.
#'
#' @param los_vec Numeric vector of LOS values (days).
#' @param title Character. Plot title.
#' @param bins Integer. Number of histogram bins. Default 35.
#' @param fits Optional. Pre-computed fits from \code{fit_distributions()}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' plot_los_distributions(rlnorm(200, 2, 0.5), "Example LOS Distribution")
plot_los_distributions <- function(los_vec, title, bins = 35, fits = NULL) {
  if (is.null(fits)) fits <- fit_distributions(los_vec)

  overlay_cfg <- list(
    weibull = list(
      dfun = stats::dweibull, color = "red",
      linetype = "solid", label = "Weibull"
    ),
    lnorm = list(
      dfun = stats::dlnorm, color = "blue",
      linetype = "dashed", label = "Lognormal"
    ),
    gamma = list(
      dfun = stats::dgamma, color = "darkgreen",
      linetype = "dotdash", label = "Gamma"
    )
  )

  p <- ggplot2::ggplot(data.frame(los = los_vec), ggplot2::aes(x = .data$los)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins, fill = "#d9d9d9", color = "black"
    )

  fitted_labels <- character(0)
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    if (is.null(fit)) next
    cfg <- overlay_cfg[[nm]]
    p <- p + ggplot2::stat_function(
      fun = cfg$dfun,
      args = as.list(fit$estimate),
      color = cfg$color, linetype = cfg$linetype, linewidth = 1.2
    )
    fitted_labels <- c(
      fitted_labels,
      paste0(cfg$label, " (", cfg$color, ")")
    )
  }

  subtitle <- if (length(fitted_labels) > 0) {
    paste(fitted_labels, collapse = " | ")
  } else {
    "No distributions could be fitted"
  }

  p + ggplot2::labs(
    title = title, subtitle = subtitle,
    x = "LOS (days)", y = "Density"
  ) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 18),
      plot.subtitle = ggplot2::element_text(size = 14, color = "grey40")
    )
}


#' Prepare LOS Dataset
#'
#' Calculates length of stay from admission and outcome dates, filters to
#' discharged patients with valid LOS within a plausible range.
#'
#' @param data Data frame with patient-level records.
#' @param admission_col Character. Column name for admission date.
#'   Default \code{"date_of_admission"}.
#' @param outcome_date_col Character. Column name for outcome/discharge date.
#'   Default \code{"final_outcome_date"}.
#' @param outcome_col Character. Column name for outcome status.
#'   Default \code{"final_outcome"}.
#' @param patient_id_col Character. Column name for patient identifier.
#'   Default \code{"PatientInformation_id"}.
#' @param max_los Numeric. Maximum plausible LOS in days. Default 200.
#'
#' @return Data frame with one row per patient-admission, including
#'   \code{LOS_days}.
#' @export
prepare_los_data <- function(data,
                             admission_col = "date_of_admission",
                             outcome_date_col = "final_outcome_date",
                             outcome_col = "final_outcome",
                             patient_id_col = "PatientInformation_id",
                             max_los = 200) {
  data %>%
    dplyr::mutate(
      .adm_date = as.Date(.data[[admission_col]]),
      .out_date = as.Date(.data[[outcome_date_col]]),
      LOS_days  = as.numeric(.out_date - .adm_date)
    ) %>%
    dplyr::filter(
      .data[[outcome_col]] == "Discharged",
      !is.na(.data$LOS_days),
      .data$LOS_days > 0,
      .data$LOS_days <= max_los
    ) %>%
    dplyr::distinct(
      .data[[patient_id_col]],
      .data$.adm_date,
      .data$.out_date,
      .data$LOS_days
    ) %>%
    dplyr::select(-".adm_date", -".out_date")
}


#' Extract LOS Vectors by Resistance Status
#'
#' For a given organism (and optionally center), computes per-patient resistance
#' status and returns separate LOS vectors for resistant and susceptible
#' exposures.
#'
#' @param abx_data Data frame of antibiotic susceptibility results with columns
#'   for patient ID, organism, antibiotic name, and susceptibility value
#'   (already cleaned to "R"/"S").
#' @param los_data Data frame from \code{prepare_los_data()}.
#' @param organism Character. Organism name (lowercase, trimmed) to filter on.
#' @param center Optional character. Center name to filter on. Default
#'   \code{NULL} (all centers).
#' @param patient_id_col Character. Patient ID column. Default
#'   \code{"PatientInformation_id"}.
#' @param organism_col Character. Organism column. Default
#'   \code{"organism_clean"}.
#' @param center_col Character. Center column. Default \code{"center_name"}.
#' @param antibiotic_col Character. Antibiotic name column. Default
#'   \code{"antibiotic_name"}.
#' @param sus_col Character. Susceptibility column (values "R"/"S"). Default
#'   \code{"sus"}.
#'
#' @return List with elements \code{R} (LOS vector for resistant exposures) and
#'   \code{S} (LOS vector for susceptible exposures).
#' @export
get_los_by_resistance <- function(abx_data,
                                  los_data,
                                  organism,
                                  center = NULL,
                                  patient_id_col = "PatientInformation_id",
                                  organism_col = "organism_clean",
                                  center_col = "center_name",
                                  antibiotic_col = "antibiotic_name",
                                  sus_col = "sus") {
  org_abx <- abx_data %>%
    dplyr::filter(.data[[organism_col]] == organism)

  if (!is.null(center)) {
    org_abx <- dplyr::filter(org_abx, .data[[center_col]] == center)
    los_data <- dplyr::filter(los_data, .data[[center_col]] == center)
  }

  patient_resistance <- org_abx %>%
    dplyr::group_by(.data[[patient_id_col]], .data[[antibiotic_col]]) %>%
    dplyr::summarise(
      final_sus = ifelse(any(.data[[sus_col]] == "R"), "R", "S"),
      .groups = "drop"
    ) %>%
    dplyr::count(.data[[patient_id_col]], .data$final_sus) %>%
    tidyr::pivot_wider(
      names_from = "final_sus", values_from = "n",
      values_fill = 0
    )

  # Handle case where R or S column may not exist
  if (!"R" %in% names(patient_resistance)) {
    patient_resistance$R <- 0L
  }
  if (!"S" %in% names(patient_resistance)) {
    patient_resistance$S <- 0L
  }

  patient_resistance <- dplyr::rename(patient_resistance,
    R_count = "R", S_count = "S"
  )

  org_los <- dplyr::inner_join(los_data, patient_resistance,
    by = patient_id_col
  )

  los_R <- org_los %>%
    dplyr::filter(.data$R_count > 0) %>%
    tidyr::uncount(.data$R_count) %>%
    dplyr::pull("LOS_days")

  los_S <- org_los %>%
    dplyr::filter(.data$S_count > 0) %>%
    tidyr::uncount(.data$S_count) %>%
    dplyr::pull("LOS_days")

  list(R = los_R, S = los_S)
}
