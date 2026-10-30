# profile_render_report.R
# Fixed-order evidence report + release gate. See profile_fit.R's header
# for how this fits into the pipeline; this is the LAST stage, assembling
# outputs already produced elsewhere in this family -- it computes
# nothing new itself except the release gate's own pass/fail logic and
# the artifact manifest.
#
# The report order (Section M of the implementation spec) is FIXED, not
# a caller preference:
#   1. Scientific question, panel and cohort
#   2. Data flow, exclusions and hospital support
#   3. Assumption register and prior predictive assessment
#   4. Model equation, shared/hospital-specific terms, coding and priors
#   5. Sampling configuration, progress/completion and computational diagnostics
#   6. Hospital and pooled PPC, with unavailable families visible
#   7. Held-out validation and its scope/uncertainty
#   8. Profile estimates with numerical CrIs and plots, exploratory or approved
#   9. Sensitivity results and unresolved limitations
#   10. Named-use decision and burden handoff
#
# A report whose evidence is incomplete, whose sampler status is "fail",
# or whose named-use review is absent is never silently presented as
# approved: release_status is "blocked" (missing/failed evidence),
# "exploratory" (evidence present, no named-use decision yet), or
# "approved" (evidence present AND an explicit, recorded named-use
# decision) -- exploratory output stays available, just visibly
# qualified, never hidden.


#' (Internal) Escape Text for Inclusion in HTML
#' @keywords internal
.profile_html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}


#' (Internal) Render a Data Frame as a Plain HTML Table
#' @keywords internal
.profile_html_table <- function(df, max_rows = 200L) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return("<p><em>Not available.</em></p>")
  }
  shown <- utils::head(df, max_rows)
  header <- paste0("<tr>", paste0("<th>", .profile_html_escape(names(shown)), "</th>", collapse = ""), "</tr>")
  body_rows <- apply(shown, 1L, function(r) {
    paste0("<tr>", paste0("<td>", .profile_html_escape(r), "</td>", collapse = ""), "</tr>")
  })
  note <- if (nrow(df) > max_rows) sprintf("<p><em>Showing first %d of %d rows.</em></p>", max_rows, nrow(df)) else ""
  paste0("<table border='1' cellpadding='4' cellspacing='0'>", header, paste(body_rows, collapse = ""), "</table>", note)
}


#' (Internal) Release-Gate Decision
#' @keywords internal
.profile_release_gate <- function(assumptions, sampler_check, named_use_review) {
  blocking_reasons <- character(0)
  if (is.null(sampler_check)) {
    blocking_reasons <- c(blocking_reasons, "sampler_check not supplied -- computational status unknown")
  } else if (identical(sampler_check$status, "fail")) {
    blocking_reasons <- c(blocking_reasons, "sampler_check status is \"fail\"")
  }
  if (!is.null(assumptions)) {
    violated <- assumptions$assumption_id[assumptions$status == "violated"]
    if (length(violated) > 0L) {
      blocking_reasons <- c(blocking_reasons, sprintf("assumption(s) violated: %s", paste(violated, collapse = ", ")))
    }
  }

  if (length(blocking_reasons) > 0L) {
    return(list(release_status = "blocked", reasons = blocking_reasons))
  }
  if (is.null(named_use_review) || !isTRUE(named_use_review$approved)) {
    return(list(release_status = "exploratory", reasons = "No recorded named-use approval -- evidence may be complete, but release requires an explicit, named decision (see named_use_review)."))
  }
  list(release_status = "approved", reasons = character(0))
}


#' Fixed-Order Evidence Report, With a Release Gate
#'
#' Assembles whatever this family's other functions have ALREADY
#' produced into the fixed 10-section report order (see this file's
#' header), writes every tabular piece supplied to \code{output_dir} as
#' CSV/RDS plus an artifact manifest (relative paths and md5 content
#' hashes), and renders one navigable HTML file. Computes nothing new
#' except the release gate (see \code{release_status} below) and the
#' manifest -- every number in the report comes from an argument the
#' caller already computed elsewhere in this family.
#'
#' @param cohort_summary Output of
#'   \code{\link{profile_summarize_cohort_eligibility}}, or \code{NULL}.
#'   Feeds sections 1-2.
#' @param assumptions Output of \code{\link{profile_check_assumptions}},
#'   or \code{NULL}. Feeds section 3, and the release gate (any
#'   \code{status == "violated"} row blocks release).
#' @param prior_predictive Output of
#'   \code{\link{profile_prior_predictive}}, or \code{NULL}. Feeds
#'   section 3.
#' @param design Output of \code{\link{profile_build_design}}, or
#'   \code{NULL}. Feeds section 4 (\code{$manifest}, the term dictionary).
#' @param sampler_check Output of \code{\link{profile_check_sampler}}, or
#'   \code{NULL}. Feeds section 5, and the release gate (a
#'   \code{status == "fail"}, or this argument being entirely absent,
#'   blocks release -- computational status unknown is never treated as
#'   computational status fine).
#' @param ppc Output of \code{\link{profile_ppc}}, or \code{NULL}. Feeds
#'   section 6. Any family not supplied (e.g. \code{count_groups} never
#'   computed) is reported as "Not available", never silently absent
#'   from the report's structure.
#' @param heldout Output of \code{\link{profile_validate_masked}},
#'   \code{\link{profile_validate_heldout}}, \code{\link{profile_validate_newpatient}},
#'   or \code{\link{profile_validate_temporal}}, or \code{NULL}. Feeds
#'   section 7.
#' @param estimates Output of \code{\link{profile_estimate}}, or
#'   \code{NULL}. Feeds section 8.
#' @param sensitivity Output of \code{\link{profile_sensitivity}}, or
#'   \code{NULL}. Feeds section 9.
#' @param release_table Output of \code{\link{profile_release}}, or
#'   \code{NULL}. Feeds sections 9-10 (the per-group decision table).
#' @param named_use_review List or \code{NULL}. The human reviewer's
#'   recorded decision: \code{list(approved = TRUE/FALSE, named_use =
#'   "<text>", reviewer = "<text>", date = "<text>", notes = "<text>")}.
#'   \code{NULL} (default) means no review has been recorded yet --
#'   release_status is then at best \code{"exploratory"}, never
#'   \code{"approved"}, no matter how clean the rest of the evidence is.
#' @param output_dir Character. Directory for the HTML report, every
#'   supplied table's CSV/RDS, and \code{manifest.csv}.
#'
#' @return A list: \code{release_status} (\code{"blocked"},
#'   \code{"exploratory"}, or \code{"approved"}), \code{release_reasons}
#'   (character vector, empty only when \code{"approved"}),
#'   \code{report_path} (the written HTML file), \code{manifest} (data
#'   frame: \code{file}, \code{relative_path}, \code{md5}).
#' @export
profile_render_report <- function(cohort_summary = NULL,
                                  assumptions = NULL,
                                  prior_predictive = NULL,
                                  design = NULL,
                                  sampler_check = NULL,
                                  ppc = NULL,
                                  heldout = NULL,
                                  estimates = NULL,
                                  sensitivity = NULL,
                                  release_table = NULL,
                                  named_use_review = NULL,
                                  output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  gate <- .profile_release_gate(assumptions, sampler_check, named_use_review)

  manifest_rows <- list()
  .write_artifact <- function(obj, name) {
    if (is.null(obj) || !is.data.frame(obj)) return(invisible(NULL))
    csv_path <- file.path(output_dir, paste0(name, ".csv"))
    rds_path <- file.path(output_dir, paste0(name, ".rds"))
    utils::write.csv(obj, csv_path, row.names = FALSE)
    saveRDS(obj, rds_path)
    for (p in c(csv_path, rds_path)) {
      manifest_rows[[length(manifest_rows) + 1]] <<- data.frame(
        file = name, relative_path = basename(p),
        md5 = unname(tools::md5sum(p)), stringsAsFactors = FALSE
      )
    }
  }

  .write_artifact(cohort_summary$summary, "cohort_summary")
  .write_artifact(assumptions, "assumptions")
  .write_artifact(prior_predictive$class_margins, "prior_predictive_class_margins")
  .write_artifact(if (!is.null(design)) design$manifest else NULL, "term_dictionary")
  .write_artifact(if (!is.null(sampler_check)) sampler_check$checks else NULL, "sampler_checks")
  .write_artifact(if (!is.null(ppc)) ppc$marginal else NULL, "ppc_marginal")
  .write_artifact(if (!is.null(ppc)) ppc$pairwise else NULL, "ppc_pairwise")
  .write_artifact(if (!is.null(heldout)) heldout$summary else NULL, "heldout_summary")
  .write_artifact(if (!is.null(estimates)) as.data.frame(estimates$profile_probability) else NULL, "profile_estimates")
  .write_artifact(sensitivity, "sensitivity")
  .write_artifact(release_table, "release_table")

  manifest <- if (length(manifest_rows) > 0L) do.call(rbind, manifest_rows) else {
    data.frame(file = character(0), relative_path = character(0), md5 = character(0), stringsAsFactors = FALSE)
  }
  rownames(manifest) <- NULL
  utils::write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names = FALSE)

  gate_banner <- switch(gate$release_status,
    blocked = "<div style='background:#fee;border:2px solid #c00;padding:10px;'><strong>BLOCKED</strong> -- not released. Reasons: ",
    exploratory = "<div style='background:#ffe;border:2px solid #cc0;padding:10px;'><strong>EXPLORATORY</strong> -- evidence present, no named-use approval recorded yet. ",
    approved = "<div style='background:#efe;border:2px solid #0a0;padding:10px;'><strong>APPROVED</strong> for the recorded named use. "
  )
  gate_html <- paste0(
    gate_banner,
    .profile_html_escape(paste(gate$reasons, collapse = "; ")),
    "</div>"
  )

  sections <- c(
    "<h2>1. Scientific question, panel and cohort</h2>", .profile_html_table(cohort_summary$summary),
    "<h2>2. Data flow, exclusions and hospital support</h2>",
    .profile_html_table(cohort_summary$class), .profile_html_table(cohort_summary$pair),
    "<h2>3. Assumption register and prior predictive assessment</h2>",
    .profile_html_table(assumptions), .profile_html_table(prior_predictive$class_margins),
    .profile_html_table(prior_predictive$concentration),
    "<h2>4. Model equation, terms, coding and priors</h2>",
    if (!is.null(design)) .profile_html_table(design$manifest) else "<p><em>Not available.</em></p>",
    "<h2>5. Sampling configuration and computational diagnostics</h2>",
    if (!is.null(sampler_check)) .profile_html_table(sampler_check$checks) else "<p><em>Not available.</em></p>",
    if (!is.null(sampler_check)) .profile_html_table(sampler_check$mcse) else "",
    "<h2>6. Hospital and pooled posterior predictive checks</h2>",
    if (!is.null(ppc)) .profile_html_table(ppc$marginal) else "<p><em>Not available.</em></p>",
    if (!is.null(ppc)) .profile_html_table(ppc$pairwise) else "",
    if (!is.null(ppc)) .profile_html_table(ppc$count_groups) else "<p><em>count_groups not available.</em></p>",
    if (!is.null(ppc)) .profile_html_table(ppc$partial_patterns) else "<p><em>partial_patterns not available.</em></p>",
    "<h2>7. Held-out validation</h2>",
    if (!is.null(heldout)) .profile_html_table(heldout$summary) else "<p><em>Not available.</em></p>",
    if (!is.null(heldout)) .profile_html_table(heldout$calibration) else "",
    "<h2>8. Profile estimates</h2>",
    if (!is.null(estimates)) sprintf("<p>Estimand: %s</p>", .profile_html_escape(estimates$estimand)) else "<p><em>Not available.</em></p>",
    "<h2>9. Sensitivity results and unresolved limitations</h2>",
    .profile_html_table(sensitivity),
    "<h2>10. Named-use decision and burden handoff</h2>",
    .profile_html_table(release_table),
    if (!is.null(named_use_review)) {
      sprintf("<p>Named use: %s | Reviewer: %s | Date: %s | Notes: %s</p>",
              .profile_html_escape(named_use_review$named_use %||% "NA"),
              .profile_html_escape(named_use_review$reviewer %||% "NA"),
              .profile_html_escape(named_use_review$date %||% "NA"),
              .profile_html_escape(named_use_review$notes %||% "NA"))
    } else {
      "<p><em>No named-use review recorded.</em></p>"
    }
  )

  html <- paste0(
    "<!DOCTYPE html><html><head><meta charset='utf-8'><title>Resistance Profile Evidence Report</title>",
    "<style>body{font-family:sans-serif;margin:30px;} table{border-collapse:collapse;margin-bottom:20px;} th,td{padding:4px 8px;}</style>",
    "</head><body>",
    "<h1>Resistance Profile Evidence Report</h1>",
    gate_html,
    paste(sections, collapse = "\n"),
    "</body></html>"
  )
  report_path <- file.path(output_dir, "report.html")
  writeLines(html, report_path)

  list(release_status = gate$release_status, release_reasons = gate$reasons, report_path = report_path, manifest = manifest)
}
