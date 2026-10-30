# profile_release.R
# Final per-group release table, combining pre-fit support with post-fit
# model evidence. See profile_fit.R's header for how this fits into the
# pipeline, and the module header below for why these two kinds of
# evidence are never collapsed into one field.
#
# "A hospital with strong initial support may still have a poorly fitting
# model. A hospital with weak support may have narrow model-based
# intervals that depend strongly on shared structure." Keep both:
#   Initial data-support tier       (profile_define_cohort()'s manifest)
#   Final estimate qualification    (profile_check_sampler() + profile_ppc()
#                                     + profile_validate_masked()/heldout()
#                                     + profile_sensitivity())
# Manual judgement is recorded as its own field with a reason, never by
# overwriting any of the above.


#' (Internal) Roll Up PPC Status for One Group
#' @keywords internal
.profile_rollup_ppc_group <- function(ppc, g) {
  flagged <- character(0)
  any_evaluated <- FALSE
  for (lvl in c("marginal", "pairwise", "complete_profile")) {
    tbl <- ppc[[lvl]]
    if (is.null(tbl)) next
    sub <- tbl[tbl$group == g & tbl$status == "evaluated", , drop = FALSE]
    if (nrow(sub) == 0L) next
    any_evaluated <- TRUE
    bad <- sub[sub$interval_contains_observed == FALSE, , drop = FALSE]
    if (nrow(bad) > 0L) {
      labels <- if (lvl == "marginal") {
        bad$class
      } else if (lvl == "pairwise") {
        paste(bad$class_a, bad$class_b, sep = "-")
      } else {
        bad$state
      }
      flagged <- c(flagged, paste0(lvl, ":", labels))
    }
  }
  status <- if (!any_evaluated) {
    "not_computable"
  } else if (length(flagged) == 0L) {
    "adequate"
  } else {
    "miscalibrated"
  }
  list(status = status, flagged = paste(flagged, collapse = "; "))
}


#' (Internal) Roll Up Held-Out/Masked Validation for One Group
#' @keywords internal
.profile_rollup_heldout_group <- function(heldout, g) {
  if (is.null(heldout) || is.null(heldout$predictions)) {
    return(list(n_masked = NA_integer_, accuracy = NA_real_, brier = NA_real_, auroc = NA_real_))
  }
  sub <- heldout$predictions[heldout$predictions$group == g, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(list(n_masked = 0L, accuracy = NA_real_, brier = NA_real_, auroc = NA_real_))
  }
  list(
    n_masked = nrow(sub),
    accuracy = mean(sub$predicted_class == sub$true_value),
    brier = mean(sub$brier_score),
    auroc = .auroc(sub$true_value, sub$predicted_resistance_probability)
  )
}


#' Final Per-Group Release Table
#'
#' Combines \code{\link{profile_define_cohort}}'s pre-fit support tier
#' with every post-fit check this family produces --
#' \code{\link{profile_check_sampler}}, \code{\link{profile_ppc}},
#' \code{\link{profile_validate_masked}}/\code{\link{profile_validate_heldout}},
#' \code{\link{profile_sensitivity}} -- into one per-group table, WITHOUT
#' collapsing any of them into a single trust label. Every post-fit
#' argument is optional and independently \code{NULL}-able: a check that
#' was not run is reported as \code{NA}/\code{"not_evaluated"}, never
#' silently omitted or inferred.
#'
#' \code{manual_judgment}, if supplied, is recorded as its own column with
#' its own reason text -- it never overwrites \code{initial_tier_reason},
#' \code{sampler_status}, \code{ppc_status}, or any other field computed
#' from actual evidence above it.
#'
#' @param cohort_manifest Data frame. The \code{$manifest} element of
#'   \code{\link{profile_define_cohort}}'s return value (needs
#'   \code{group}, \code{included}, \code{reason}).
#' @param sampler_check List or \code{NULL}. Return value of
#'   \code{\link{profile_check_sampler}} (one verdict for the whole fit,
#'   repeated on every row -- a sampler either worked or it didn't,
#'   regardless of group).
#' @param ppc List or \code{NULL}. Return value of \code{\link{profile_ppc}}
#'   (\code{$marginal}/\code{$pairwise}/\code{$complete_profile}).
#' @param heldout List or \code{NULL}. Return value of
#'   \code{\link{profile_validate_masked}} or
#'   \code{\link{profile_validate_heldout}} (reads \code{$predictions}).
#' @param sensitivity Data frame or \code{NULL}. Return value of
#'   \code{\link{profile_sensitivity}}.
#' @param manual_judgment Named list (\code{group -> reason}), data frame
#'   (\code{group}, \code{reason} columns), or \code{NULL}. Analyst
#'   judgment, recorded verbatim with its stated reason -- never generated
#'   by this function.
#' @param output_dir Character or \code{NULL}. If supplied, writes the
#'   returned table as \code{release_table.csv} and \code{.rds}.
#'
#' @return A data frame, one row per group appearing in ANY of the
#'   supplied inputs (a group the cohort manifest excluded but that still
#'   shows up in, say, \code{ppc} is not dropped): \code{group},
#'   \code{initial_included}, \code{initial_tier_reason},
#'   \code{sampler_status}, \code{sampler_reasons}, \code{ppc_status},
#'   \code{ppc_flagged_checks}, \code{heldout_n_masked},
#'   \code{heldout_accuracy}, \code{heldout_brier}, \code{heldout_auroc},
#'   \code{max_sensitivity_contrast}, \code{manual_judgment_reason}.
#' @export
profile_release <- function(cohort_manifest,
                            sampler_check = NULL,
                            ppc = NULL,
                            heldout = NULL,
                            sensitivity = NULL,
                            manual_judgment = NULL,
                            output_dir = NULL) {
  if (!is.data.frame(cohort_manifest) || !all(c("group", "included") %in% names(cohort_manifest))) {
    stop(
      "cohort_manifest must be the manifest data frame from profile_define_cohort() ",
      "(needs group, included, reason)."
    )
  }

  all_groups <- unique(c(
    cohort_manifest$group,
    if (!is.null(ppc$marginal)) ppc$marginal$group,
    if (!is.null(ppc$pairwise)) ppc$pairwise$group,
    if (!is.null(ppc$complete_profile)) ppc$complete_profile$group,
    if (!is.null(heldout$predictions)) heldout$predictions$group,
    if (!is.null(sensitivity)) sensitivity$group
  ))

  rows <- lapply(all_groups, function(g) {
    cm <- cohort_manifest[cohort_manifest$group == g, , drop = FALSE]
    initial_included <- if (nrow(cm) > 0L) cm$included[1] else NA
    initial_tier_reason <- if (nrow(cm) > 0L) cm$reason[1] else "not_in_cohort_manifest"

    sampler_status <- if (!is.null(sampler_check)) sampler_check$status else NA_character_
    sampler_reasons <- if (!is.null(sampler_check)) paste(sampler_check$reasons, collapse = "; ") else NA_character_

    ppc_roll <- if (!is.null(ppc)) {
      .profile_rollup_ppc_group(ppc, g)
    } else {
      list(status = "not_evaluated", flagged = NA_character_)
    }

    held_roll <- .profile_rollup_heldout_group(heldout, g)

    sens_val <- if (!is.null(sensitivity)) {
      sub <- sensitivity[sensitivity$group == g, , drop = FALSE]
      if (nrow(sub) > 0L) max(sub$sensitivity_contrast, na.rm = TRUE) else NA_real_
    } else {
      NA_real_
    }

    mj <- NA_character_
    if (!is.null(manual_judgment)) {
      if (is.data.frame(manual_judgment)) {
        m <- manual_judgment[manual_judgment$group == g, , drop = FALSE]
        if (nrow(m) > 0L) mj <- m$reason[1]
      } else if (is.list(manual_judgment) && g %in% names(manual_judgment)) {
        mj <- manual_judgment[[g]]
      }
    }

    data.frame(
      group = g,
      initial_included = initial_included,
      initial_tier_reason = initial_tier_reason,
      sampler_status = sampler_status,
      sampler_reasons = sampler_reasons,
      ppc_status = ppc_roll$status,
      ppc_flagged_checks = ppc_roll$flagged,
      heldout_n_masked = held_roll$n_masked,
      heldout_accuracy = held_roll$accuracy,
      heldout_brier = held_roll$brier,
      heldout_auroc = held_roll$auroc,
      max_sensitivity_contrast = sens_val,
      manual_judgment_reason = mj,
      stringsAsFactors = FALSE
    )
  })

  release_tbl <- do.call(rbind, rows)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(release_tbl, file.path(output_dir, "release_table.csv"), row.names = FALSE)
    saveRDS(release_tbl, file.path(output_dir, "release_table.rds"))
  }

  release_tbl
}
