# report.R
# Data provenance tracking and reporting for AMR preprocessing

#' Generate Preprocessing Report
#'
#' Creates a comprehensive report documenting all transformations applied
#' during preprocessing, including column mappings, data quality metrics,
#' and transformation summaries.
#'
#' @param raw_data Original input data frame (before preprocessing)
#' @param processed_data Final output data frame (after preprocessing)
#' @param preprocessing_log List containing logs from each preprocessing step
#' @param config AMR configuration object used
#' @param output_file Optional. Path to save report. If NULL, returns list.
#'   Supports .html, .pdf, .txt, .rds formats.
#' @param include_plots Logical. Include data quality plots. Default TRUE.
#'
#' @return If output_file is NULL, returns a list with report components.
#'   Otherwise, saves report to file and returns file path.
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate report
#' report <- generate_preprocessing_report(
#'   raw_data = my_raw_data,
#'   processed_data = result$data,
#'   preprocessing_log = result$log,
#'   config = result$config,
#'   output_file = "preprocessing_report.html"
#' )
#' }
generate_preprocessing_report <- function(raw_data,
                                          processed_data,
                                          preprocessing_log,
                                          config,
                                          output_file = NULL,
                                          include_plots = TRUE) {
  report <- list(
    metadata = list(
      report_generated = Sys.time(),
      package_version = utils::packageVersion("anumaan"),
      r_version = R.version.string
    )
  )

  # Section 1: Raw Data Summary
  report$raw_data_summary <- summarize_raw_data(raw_data)

  # Section 2: Column Mapping Details
  report$column_mapping <- summarize_column_mapping(
    raw_data,
    processed_data,
    preprocessing_log$column_mapping
  )

  # Section 3: Data Transformation Summary
  report$transformations <- summarize_transformations(
    raw_data,
    processed_data,
    preprocessing_log
  )

  # Section 4: Data Quality Metrics
  report$data_quality <- summarize_data_quality(
    raw_data,
    processed_data,
    include_plots = include_plots
  )

  # Section 5: Configuration Used
  report$config_used <- config

  # Section 6: Processing Statistics
  report$processing_stats <- compute_processing_statistics(
    raw_data,
    processed_data,
    preprocessing_log
  )

  # Section 7: Warnings and Issues
  report$warnings <- preprocessing_log$warnings %||% list()
  report$errors <- preprocessing_log$errors %||% list()

  # Add class
  class(report) <- c("amr_preprocessing_report", "list")

  # Save to file if requested
  if (!is.null(output_file)) {
    export_report(report, output_file)
    message(sprintf("Report saved to: %s", output_file))
    return(invisible(output_file))
  }

  return(report)
}


#' Summarize Raw Data
#'
#' Creates summary statistics for the original input dataset.
#'
#' @param data Raw input data frame
#'
#' @return List with raw data summaries
#' @export
summarize_raw_data <- function(data) {
  summary <- list(
    n_rows = nrow(data),
    n_cols = ncol(data),
    column_names = names(data),
    column_types = sapply(data, class),
    memory_size = format(object.size(data), units = "auto")
  )

  # Missing data summary
  summary$missing_data <- data.frame(
    column = names(data),
    n_missing = sapply(data, function(x) sum(is.na(x))),
    pct_missing = sapply(data, function(x) 100 * mean(is.na(x))),
    stringsAsFactors = FALSE
  )
  rownames(summary$missing_data) <- NULL

  # Identify likely key columns
  summary$key_columns <- detect_key_columns(data)

  # Sample data (first 5 rows)
  summary$sample_data <- utils::head(data, 5)

  return(summary)
}


#' Summarize Column Mapping
#'
#' Documents how original column names were mapped to standard names.
#'
#' @param raw_data Original data
#' @param processed_data Processed data
#' @param mapping_log Column mapping log from preprocessing
#'
#' @return Data frame with mapping details
#' @export
summarize_column_mapping <- function(raw_data, processed_data, mapping_log) {
  raw_cols <- names(raw_data)
  processed_cols <- names(processed_data)

  # Build mapping table
  mapping_df <- data.frame(
    original_name = character(),
    standardized_name = character(),
    mapping_method = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )

  # Add mapped columns
  if (!is.null(mapping_log)) {
    for (std_name in names(mapping_log)) {
      entry <- mapping_log[[std_name]]
      mapping_df <- rbind(mapping_df, data.frame(
        original_name = entry$original,
        standardized_name = std_name,
        mapping_method = entry$method,
        notes = if (!is.null(entry$distance)) {
          sprintf("fuzzy match (distance: %.2f)", entry$distance)
        } else {
          "exact match"
        },
        stringsAsFactors = FALSE
      ))
    }
  }

  # Add unmapped columns (stayed the same)
  unmapped_original <- setdiff(raw_cols, mapping_df$original_name)
  unmapped_in_processed <- intersect(unmapped_original, processed_cols)

  if (length(unmapped_in_processed) > 0) {
    unmapped_df <- data.frame(
      original_name = unmapped_in_processed,
      standardized_name = unmapped_in_processed,
      mapping_method = "no_mapping",
      notes = "kept original name",
      stringsAsFactors = FALSE
    )
    mapping_df <- rbind(mapping_df, unmapped_df)
  }

  # Add new columns (created during processing)
  new_cols <- setdiff(processed_cols, raw_cols)
  if (length(new_cols) > 0) {
    new_df <- data.frame(
      original_name = NA_character_,
      standardized_name = new_cols,
      mapping_method = "derived",
      notes = "created during preprocessing",
      stringsAsFactors = FALSE
    )
    mapping_df <- rbind(mapping_df, new_df)
  }

  # Add dropped columns (in raw but not in processed)
  dropped_cols <- setdiff(raw_cols, processed_cols)
  if (length(dropped_cols) > 0) {
    dropped_df <- data.frame(
      original_name = dropped_cols,
      standardized_name = NA_character_,
      mapping_method = "dropped",
      notes = "removed during preprocessing",
      stringsAsFactors = FALSE
    )
    mapping_df <- rbind(mapping_df, dropped_df)
  }

  return(mapping_df)
}


#' Summarize Data Transformations
#'
#' Documents all transformations applied to data values.
#'
#' @param raw_data Original data
#' @param processed_data Processed data
#' @param preprocessing_log Full preprocessing log
#'
#' @return List with transformation summaries
#' @export
summarize_transformations <- function(raw_data, processed_data, preprocessing_log) {
  transformations <- list()

  # Gender standardization
  if ("gender" %in% names(processed_data)) {
    transformations$gender <- summarize_value_mapping(
      raw_data,
      processed_data,
      find_raw_column(raw_data, c("gender", "Gender", "sex", "Sex")),
      "gender"
    )
  }

  # Outcome standardization
  if ("final_outcome" %in% names(processed_data)) {
    transformations$outcome <- summarize_value_mapping(
      raw_data,
      processed_data,
      find_raw_column(raw_data, c("final_outcome", "Final.outcome", "outcome")),
      "final_outcome"
    )
  }

  # Susceptibility standardization
  if ("antibiotic_value_std" %in% names(processed_data)) {
    transformations$susceptibility <- summarize_value_mapping(
      raw_data,
      processed_data,
      find_raw_column(raw_data, c("antibiotic_value", "remarks", "Result")),
      "antibiotic_value_std"
    )
  }

  # Date parsing
  date_cols <- c("date_of_admission", "date_of_culture", "date_of_final_outcome", "DOB")
  transformations$dates <- list()
  for (col in date_cols) {
    if (col %in% names(processed_data)) {
      raw_col <- find_raw_column(raw_data, col)
      if (!is.null(raw_col)) {
        transformations$dates[[col]] <- list(
          raw_class = class(raw_data[[raw_col]])[1],
          processed_class = class(processed_data[[col]])[1],
          n_parsed = sum(!is.na(processed_data[[col]])),
          n_failed = sum(is.na(processed_data[[col]]))
        )
      }
    }
  }

  # Organism normalization
  if ("organism_normalized" %in% names(processed_data)) {
    raw_org_col <- find_raw_column(
      raw_data,
      c("organism_name", "Organism", "organism", "org_name")
    )
    if (!is.null(raw_org_col)) {
      n_raw_org <- length(unique(na.omit(raw_data[[raw_org_col]])))
      n_norm_org <- length(unique(na.omit(processed_data$organism_normalized)))
      transformations$organism <- list(
        n_unique_raw        = n_raw_org,
        n_unique_normalized = n_norm_org,
        reduction_pct       = 100 * (1 - n_norm_org / n_raw_org)
      )
    }
  }

  # Antibiotic normalization
  if ("antibiotic_normalized" %in% names(processed_data)) {
    raw_abx_col <- find_raw_column(
      raw_data,
      c("antibiotic_name", "Antibiotic", "antibiotic", "drug_name", "Drug")
    )
    if (!is.null(raw_abx_col)) {
      n_raw_abx <- length(unique(na.omit(raw_data[[raw_abx_col]])))
      n_norm_abx <- length(unique(na.omit(processed_data$antibiotic_normalized)))
      transformations$antibiotic <- list(
        n_unique_raw        = n_raw_abx,
        n_unique_normalized = n_norm_abx,
        reduction_pct       = 100 * (1 - n_norm_abx / n_raw_abx)
      )
    }
  }

  return(transformations)
}


#' Summarize Data Quality
#'
#' Computes data quality metrics before and after preprocessing.
#'
#' @param raw_data Original data
#' @param processed_data Processed data
#' @param include_plots Logical. Generate plots. Default TRUE.
#'
#' @return List with quality metrics
#' @export
summarize_data_quality <- function(raw_data, processed_data, include_plots = TRUE) {
  quality <- list()

  # Completeness metrics
  quality$completeness <- data.frame(
    stage = c("raw", "processed"),
    n_rows = c(nrow(raw_data), nrow(processed_data)),
    n_cols = c(ncol(raw_data), ncol(processed_data)),
    total_values = c(
      nrow(raw_data) * ncol(raw_data),
      nrow(processed_data) * ncol(processed_data)
    ),
    missing_values = c(
      sum(is.na(raw_data)),
      sum(is.na(processed_data))
    ),
    pct_complete = c(
      100 * (1 - sum(is.na(raw_data)) / (nrow(raw_data) * ncol(raw_data))),
      100 * (1 - sum(is.na(processed_data)) / (nrow(processed_data) * ncol(processed_data)))
    ),
    stringsAsFactors = FALSE
  )

  # Data type consistency
  quality$data_types <- list(
    raw = table(sapply(raw_data, class)),
    processed = table(sapply(processed_data, class))
  )

  # Duplicate detection
  quality$duplicates <- list(
    raw_duplicate_rows = sum(duplicated(raw_data)),
    processed_duplicate_rows = sum(duplicated(processed_data))
  )

  # Outlier detection (numeric columns only)
  quality$outliers <- detect_outliers(processed_data)

  return(quality)
}


#' Compute Processing Statistics
#'
#' Calculates performance and processing statistics.
#'
#' @param raw_data Original data
#' @param processed_data Processed data
#' @param preprocessing_log Full log
#'
#' @return List with statistics
#' @export
compute_processing_statistics <- function(raw_data, processed_data, preprocessing_log) {
  stats <- list(
    runtime = preprocessing_log$runtime %||% NA,
    rows_input = nrow(raw_data),
    rows_output = nrow(processed_data),
    rows_removed = nrow(raw_data) - nrow(processed_data),
    pct_retained = 100 * nrow(processed_data) / nrow(raw_data),
    cols_input = ncol(raw_data),
    cols_output = ncol(processed_data),
    cols_added = ncol(processed_data) - ncol(raw_data),
    variables_derived = sum(grepl(
      "_method$|_confidence$|_derived$",
      names(processed_data)
    ))
  )

  return(stats)
}


#' Print Preprocessing Report
#'
#' @param x An amr_preprocessing_report object
#' @param ... Additional arguments (not used)
#'
#' @export
print.amr_preprocessing_report <- function(x, ...) {
  cat("===================================================\n")
  cat("   AMR Data Preprocessing Report\n")
  cat("===================================================\n\n")

  cat("Report Generated:", format(x$metadata$report_generated), "\n")
  cat("Package Version:", as.character(x$metadata$package_version), "\n")
  cat("R Version:", x$metadata$r_version, "\n\n")

  cat("--- Raw Data Summary ---\n")
  cat(sprintf("  Rows: %s\n", format(x$raw_data_summary$n_rows, big.mark = ",")))
  cat(sprintf("  Columns: %d\n", x$raw_data_summary$n_cols))
  cat(sprintf("  Memory: %s\n", x$raw_data_summary$memory_size))
  cat(sprintf(
    "  Missing Data: %.1f%%\n\n",
    mean(x$raw_data_summary$missing_data$pct_missing)
  ))

  cat("--- Column Mapping ---\n")
  mapping_summary <- table(x$column_mapping$mapping_method)
  for (method in names(mapping_summary)) {
    cat(sprintf("  %s: %d columns\n", method, mapping_summary[method]))
  }
  cat("\n")

  cat("--- Processing Statistics ---\n")
  cat(sprintf("  Runtime: %.2f seconds\n", x$processing_stats$runtime))
  cat(sprintf(
    "  Rows: %s -> %s (%.1f%% retained)\n",
    format(x$processing_stats$rows_input, big.mark = ","),
    format(x$processing_stats$rows_output, big.mark = ","),
    x$processing_stats$pct_retained
  ))
  cat(sprintf(
    "  Columns: %d -> %d (+%d derived)\n",
    x$processing_stats$cols_input,
    x$processing_stats$cols_output,
    x$processing_stats$cols_added
  ))
  cat(sprintf(
    "  Variables Derived: %d\n\n",
    x$processing_stats$variables_derived
  ))

  cat("--- Data Quality ---\n")
  cat(sprintf(
    "  Completeness: %.1f%% -> %.1f%%\n",
    x$data_quality$completeness$pct_complete[1],
    x$data_quality$completeness$pct_complete[2]
  ))
  cat(sprintf(
    "  Duplicate Rows: %d -> %d\n",
    x$data_quality$duplicates$raw_duplicate_rows,
    x$data_quality$duplicates$processed_duplicate_rows
  ))

  if (length(x$warnings) > 0) {
    cat(sprintf("\n[!] Warnings: %d (see $warnings for details)\n", length(x$warnings)))
  }

  if (length(x$errors) > 0) {
    cat(sprintf("[x] Errors: %d (see $errors for details)\n", length(x$errors)))
  }

  cat("\n===================================================\n")
  cat("Use summary(report) for detailed sections\n")
  cat("Use export_report(report, 'file.html') to save\n")
  cat("===================================================\n")

  invisible(x)
}


#' Export Report to File
#'
#' Saves preprocessing report to HTML, PDF, TXT, or RDS format.
#'
#' @param report An amr_preprocessing_report object
#' @param file Path to output file
#'
#' @return Invisibly returns file path
#' @export
export_report <- function(report, file) {
  ext <- tools::file_ext(file)

  switch(tolower(ext),
    html = export_report_html(report, file),
    pdf = export_report_pdf(report, file),
    txt = export_report_txt(report, file),
    rds = saveRDS(report, file),
    stop("Unsupported file format. Use .html, .pdf, .txt, or .rds")
  )

  invisible(file)
}


#' Export Report as HTML
#'
#' @param report Report object
#' @param file Output file path
#' @keywords internal
export_report_html <- function(report, file) {
  html <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "<meta charset='UTF-8'>",
    "<title>AMR Preprocessing Report</title>",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 40px; }",
    "h1 { color: #2c3e50; border-bottom: 2px solid #3498db; }",
    "h2 { color: #34495e; margin-top: 30px; }",
    "table { border-collapse: collapse; width: 100%; margin: 20px 0; }",
    "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
    "th { background-color: #3498db; color: white; }",
    "tr:nth-child(even) { background-color: #f2f2f2; }",
    ".metric { background-color: #ecf0f1; padding: 10px; margin: 10px 0; border-radius: 5px; }",
    ".warning { color: #e74c3c; }",
    ".success { color: #27ae60; }",
    "</style>",
    "</head>",
    "<body>",
    sprintf("<h1>AMR Data Preprocessing Report</h1>"),
    sprintf("<p><strong>Generated:</strong> %s</p>", report$metadata$report_generated),
    sprintf("<p><strong>Package Version:</strong> %s</p>", report$metadata$package_version),
    "<h2>1. Raw Data Summary</h2>",
    sprintf("<div class='metric'>Rows: %s</div>", format(report$raw_data_summary$n_rows, big.mark = ",")),
    sprintf("<div class='metric'>Columns: %d</div>", report$raw_data_summary$n_cols),
    "<h2>2. Column Mapping</h2>",
    convert_df_to_html_table(report$column_mapping),
    "<h2>3. Processing Statistics</h2>",
    sprintf("<div class='metric'>Runtime: %.2f seconds</div>", report$processing_stats$runtime),
    sprintf("<div class='metric'>Data Retention: %.1f%%</div>", report$processing_stats$pct_retained),
    "<h2>4. Data Quality</h2>",
    convert_df_to_html_table(report$data_quality$completeness),
    "</body>",
    "</html>"
  )

  writeLines(html, file)
}


#' Export Report as PDF
#'
#' @param report Report object
#' @param file Output file path
#' @keywords internal
export_report_pdf <- function(report, file) {
  stop(
    "PDF export requires rmarkdown and a LaTeX installation. ",
    "Use .html, .txt, or .rds format instead."
  )
}


#' Export Report as Text
#'
#' @param report Report object
#' @param file Output file path
#' @keywords internal
export_report_txt <- function(report, file) {
  sink(file)
  on.exit(sink())
  print(report)
  cat("\n\n=== DETAILED COLUMN MAPPING ===\n")
  print(report$column_mapping)
  cat("\n\n=== DATA QUALITY METRICS ===\n")
  print(report$data_quality$completeness)
  sink()
}


# ===== Helper Functions =====

#' Detect Key Columns
#' @keywords internal
detect_key_columns <- function(data) {
  key_patterns <- c("id", "patient", "event", "subject", "mrn", "uhid")
  potential_keys <- character()

  for (col in names(data)) {
    if (any(sapply(key_patterns, function(p) grepl(p, col, ignore.case = TRUE)))) {
      potential_keys <- c(potential_keys, col)
    }
  }

  return(potential_keys)
}


#' Find Raw Column by Aliases
#' @keywords internal
find_raw_column <- function(data, aliases) {
  for (alias in aliases) {
    if (alias %in% names(data)) {
      return(alias)
    }
  }
  return(NULL)
}


#' Summarize Value Mapping
#' @keywords internal
summarize_value_mapping <- function(raw_data, processed_data, raw_col, processed_col) {
  if (is.null(raw_col) || !processed_col %in% names(processed_data)) {
    return(NULL)
  }

  list(
    raw_unique_values = unique(raw_data[[raw_col]]),
    processed_unique_values = unique(processed_data[[processed_col]]),
    value_mapping = table(
      raw = raw_data[[raw_col]],
      processed = processed_data[[processed_col]]
    )
  )
}


#' Detect Outliers
#' @keywords internal
detect_outliers <- function(data) {
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  outlier_summary <- list()

  for (col in numeric_cols) {
    q1 <- quantile(data[[col]], 0.25, na.rm = TRUE)
    q3 <- quantile(data[[col]], 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower <- q1 - 1.5 * iqr
    upper <- q3 + 1.5 * iqr

    outliers <- data[[col]] < lower | data[[col]] > upper
    n_outliers <- sum(outliers, na.rm = TRUE)

    if (n_outliers > 0) {
      outlier_summary[[col]] <- list(
        n_outliers = n_outliers,
        pct_outliers = 100 * n_outliers / length(data[[col]]),
        range = range(data[[col]], na.rm = TRUE),
        outlier_range = range(data[[col]][outliers], na.rm = TRUE)
      )
    }
  }

  return(outlier_summary)
}


#' Convert Data Frame to HTML Table
#' @keywords internal
convert_df_to_html_table <- function(df) {
  if (nrow(df) == 0) {
    return("<p>No data</p>")
  }

  html <- c("<table>", "<tr>")

  # Header
  for (col in names(df)) {
    html <- c(html, sprintf("<th>%s</th>", col))
  }
  html <- c(html, "</tr>")

  # Rows (limit to first 100)
  for (i in seq_len(min(nrow(df), 100))) {
    html <- c(html, "<tr>")
    for (col in names(df)) {
      html <- c(html, sprintf("<td>%s</td>", df[i, col]))
    }
    html <- c(html, "</tr>")
  }

  html <- c(html, "</table>")
  return(paste(html, collapse = "\n"))
}
