# profile_eligibility.R
# Hospital/stratum eligibility for joint resistance-profile models, derived
# from actual data support rather than a hand-typed list.
#
# This is the piece the J2-L model-development history never actually closed
# the loop on: anumaan-analysis/Scripts/stewardship/resistance_profiles/
# methodology/j2_support_precision/ computes rigorous Wilson/Goodman
# precision for candidate assurance scenarios, but the 20-hospital cohort
# itself (joint_bernoulli_20_hospital's `j2l20_hospitals`) was a hardcoded
# literal list, built independently of that audit. These three functions
# make the threshold computation the thing that PRODUCES a roster, instead
# of commentary alongside one.
#
# Fully generic: `group_col` is a real argument, never assumed to be
# "hospital". Any stratification variable (district, lab, time period, ...)
# works identically.
#
# Functions, in pipeline order:
#   - profile_select_cohort  (stage 0: restrict to the modeling cohort --
#                             e.g. a discharged/non-fatal outcome filter --
#                             BEFORE any support/eligibility calculation, so
#                             eligibility is never computed on a bigger
#                             cohort than the one actually fit)
#   - profile_audit_support  (raw tested/resistant/susceptible counts, by
#                             group x class, by group x class-pair, and by
#                             group x complete-panel)
#   - profile_plan_precision (Wilson/Goodman interval width + pass/fail at
#                             a chosen set of precision tiers, per support row)
#   - profile_define_cohort  (the actual group roster + a traceable
#                             inclusion/exclusion manifest)
#   - profile_precision_bound (distribution-free minimum complete-event
#                             count, n_min, to simultaneously bound every
#                             profile probability of a D-class panel within
#                             a stated error -- Hoeffding + union bound,
#                             optionally Bonferroni-corrected across strata.
#                             A mathematical consequence of a chosen
#                             (epsilon, alpha), not a fitted or tuned value)
#   - profile_check_panel_feasibility (pathogen x candidate-panel x group
#                             feasibility table: n_min computed once per
#                             candidate panel, checked against each group's
#                             own complete-event count. No coverage-fraction
#                             aggregation and no auto-selected panel size --
#                             replaces the retired profile_define_panel(),
#                             which picked "the" panel size via an arbitrary
#                             0.60 coverage cutoff; see anumaan-analysis's
#                             R/stewardship/panel_feasibility.R::
#                             score_candidate_panels() for the even older,
#                             raw-count version of that same mistake)
#   - profile_audit_interaction_support (which interaction terms of a
#                             candidate panel have ANY direct co-tested
#                             support, per group, vs. relying almost
#                             entirely on the model's shared cross-class
#                             prior -- the cheap, pre-fit screen; a
#                             structural-support audit, explicitly NOT an
#                             identifiability claim. A formal
#                             identifiability analysis against a fitted
#                             model's information matrix is a separate,
#                             later concern)
#   - profile_prepare_model_data (stage 3: builds the Y/ids/covariates
#                             wide model dataset for one chosen panel,
#                             ready for profile_build_design())
#   - profile_summarize_cohort_eligibility (reporting layer: one per-group
#                             table -- event/patient counts, the
#                             complete/partial/zero-observed breakdown, and
#                             BOTH independent eligibility verdicts
#                             (direct_support_pass, included_in_fit) --
#                             plus the class/pair/interaction support
#                             tables, assembled from the functions above
#                             rather than computed fresh)


#' Select and Record the Modeling Cohort From Raw Event Data
#'
#' Stage 0 of the resistance-profiling pipeline. Every eligibility,
#' precision, and panel-feasibility calculation downstream of this function
#' must be computed on the SAME cohort that will eventually be fit -- never
#' on the full dataset while fitting only a restricted subset (e.g.
#' discharged-only events). This function makes that restriction an
#' explicit, traceable step instead of an assumption buried in a script.
#'
#' Filters are applied in the order supplied, each an exact-match inclusion
#' rule against one column. The returned manifest shows how many rows
#' survived EACH filter in sequence, so a caller can see exactly how much a
#' given restriction removed on its own, rather than only one before/after
#' total that hides which criterion did the work.
#'
#' This filters rows; it does not attach a clinical label to the result.
#' Restricting to \code{final_outcome == "Discharged"} describes exactly
#' that cohort, not automatically "all non-fatal infections" -- transfers,
#' ongoing admissions, or other surviving-but-not-discharged outcomes are
#' excluded by that specific restriction, and the caller should label the
#' resulting cohort accordingly rather than assuming the broader label.
#'
#' @param data Data frame of events/isolates.
#' @param filters Named list. Each name is a column in \code{data}; each
#'   value is a character/numeric vector of values that column must match
#'   (via \code{\%in\%}) for a row to survive. Applied in list order (e.g.
#'   \code{list(final_outcome = "Discharged", center_name = my_hospitals)}).
#' @param event_id_col Character or \code{NULL}. Reserved for future
#'   event-level manifest detail; currently unused beyond validating it is
#'   present if supplied.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{data}}{the filtered data frame}
#'     \item{\code{manifest}}{one row per filter applied, in order:
#'       \code{column}, \code{values} (comma-joined), \code{n_before},
#'       \code{n_after}, \code{n_excluded}}
#'   }
#' @export
profile_select_cohort <- function(data, filters, event_id_col = NULL) {
  if (length(filters) == 0L || is.null(names(filters)) || any(names(filters) == "")) {
    stop("filters must be a non-empty named list (column name -> allowed values)")
  }
  if (!is.null(event_id_col) && !event_id_col %in% names(data)) {
    stop(sprintf("event_id_col '%s' not found in data", event_id_col))
  }

  manifest_rows <- list()
  cur <- data
  for (col in names(filters)) {
    if (!col %in% names(cur)) {
      stop(sprintf("filter column '%s' not found in data", col))
    }
    vals <- filters[[col]]
    n_before <- nrow(cur)
    cur <- cur[cur[[col]] %in% vals, , drop = FALSE]
    n_after <- nrow(cur)
    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      column = col,
      values = paste(vals, collapse = ", "),
      n_before = n_before,
      n_after = n_after,
      n_excluded = n_before - n_after,
      stringsAsFactors = FALSE
    )
  }

  list(data = cur, manifest = do.call(rbind, manifest_rows))
}


#' (Internal) Classify and Pivot Long/Wide AST Data for the profile_* Family
#'
#' Shared by every \code{profile_*} function that needs a group x event x
#' class S/R/I matrix: \code{\link{profile_audit_support}},
#' \code{\link{profile_audit_interaction_support}}, and
#' \code{\link{profile_prepare_model_data}} all call this once instead of
#' each re-implementing the same long-to-wide classify-and-pivot. Not
#' exported -- this is specifically the classification step (mapping raw
#' values to canonical R/S/I via \code{resistant_values}/etc.), which is
#' resistance-profile-domain logic, not a generic preprocessing pivot; it is
#' deliberately NOT the same job as \code{anumaan}'s
#' \code{prep_create_wide_ast_matrix()}, which pivots raw values verbatim
#' with no classification and is for preprocessing-stage use.
#'
#' One AST policy, applied identically everywhere this helper is used. Raw
#' value parsing tries the caller's declared \code{resistant_values}/
#' \code{susceptible_values}/\code{intermediate_values} vocabulary first
#' (exact match, case-insensitive), then falls back to
#' \code{\link{prep_clean_ast_values}}'s messy-string extraction (handles
#' \code{"S (HIGH LEVEL)"}, keyword scans, etc.) for anything that didn't
#' match -- the same parsing logic used everywhere else in \code{anumaan},
#' not a separate, cruder exact-match-only implementation kept just here.
#' blank/\code{NA}/declared \code{untested_values} tokens are UNTESTED, not
#' an outcome; a non-blank value that matches none of the above (declared
#' vocabulary OR \code{prep_clean_ast_values}'s extraction) is UNKNOWN,
#' handled per \code{unknown_action} (\code{"error"}, the default -- a typo
#' or an unanticipated lab code should be loud, not silently swept into
#' "untested"; \code{"treat_as_untested"} to opt out when that is a
#' deliberate, declared choice). Duplicate rows for the SAME group x event x
#' class are resolved by delegating to \code{\link{prep_deduplicate_ast}}
#' (via \code{conflict_action}, passed through as its \code{strategy}) --
#' not a separate, parallel implementation: identical duplicates always
#' collapse to one first (an \code{n_duplicates_collapsed} count is
#' returned, never silently inflating support), and duplicates with
#' DIFFERING R/S/I results are a genuine data conflict, never resolved by
#' taking the first row unless that policy is explicitly chosen --
#' \code{"error"} (the default) stops with every conflicting cell named;
#' \code{"drop_conflicts"} removes every row for a conflicting cell (never
#' guesses which one is right); \code{"resistant_wins"}/
#' \code{"susceptible_wins"}/\code{"first"} resolve by the same clinical
#' priority or row-order policy \code{prep_deduplicate_ast} already offers
#' elsewhere in the package. Grouping is always by \code{(.group, .event)}
#' jointly, so an \code{event_id_col} value that repeats ACROSS groups (e.g.
#' a hospital-local patient id reused at a different hospital) is never
#' merged across hospitals -- group-scoped identifiers are the default, not
#' an opt-in.
#'
#' When \code{intermediate_as_resistant} is supplied (not \code{NULL}), an
#' additional \code{.result_folded} column (and its own wide pivot,
#' \code{wide_folded}/\code{result_col_folded}) is built by delegating to
#' \code{\link{prep_recode_intermediate_ast}} -- the SAME Colistin-aware
#' clinical convention already used elsewhere in \code{anumaan} (Colistin
#' \code{I -> S}, every other class \code{I -> R}), rather than a separate,
#' cruder symmetric \code{I -> R} rule. \code{FALSE} leaves
#' \code{.result_folded} identical to \code{.result} (\code{I} stays its own
#' category, as before). \code{NULL} (the default) skips this step entirely
#' -- only the callers that fold at all (\code{\link{profile_audit_support}},
#' \code{\link{profile_prepare_model_data}}) pass a non-\code{NULL} value.
#'
#' @keywords internal
.profile_tested_long_wide <- function(data, group_col, class_cols, class_col, value_col,
                                      resistant_values, susceptible_values, intermediate_values,
                                      event_id_col,
                                      untested_values = c("NT", "Not tested", "not tested", "NA", "N/A"),
                                      unknown_action = c("error", "treat_as_untested"),
                                      conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins", "first"),
                                      intermediate_as_resistant = NULL) {
  unknown_action <- match.arg(unknown_action)
  conflict_action <- match.arg(conflict_action)

  if (!group_col %in% names(data)) {
    stop(sprintf("group_col '%s' not found in data", group_col))
  }

  is_wide <- all(class_cols %in% names(data))
  is_long <- !is_wide && !is.null(class_col) && class_col %in% names(data) && value_col %in% names(data)
  if (!is_wide && !is_long) {
    stop(
      "Could not resolve class_cols: either all of class_cols must be ",
      "columns in data (wide format), or class_col + value_col must be ",
      "supplied and present (long format)."
    )
  }

  .classify <- function(x) {
    x_chr <- trimws(as.character(x))
    is_blank <- is.na(x_chr) | x_chr == "" | toupper(x_chr) %in% toupper(untested_values)
    result <- rep(NA_character_, length(x_chr))
    result[!is_blank & toupper(x_chr) %in% toupper(resistant_values)] <- "R"
    result[!is_blank & toupper(x_chr) %in% toupper(susceptible_values)] <- "S"
    result[!is_blank & toupper(x_chr) %in% toupper(intermediate_values)] <- "I"

    # Anything non-blank that didn't exact-match the caller's declared
    # vocabulary falls back to prep_clean_ast_values()'s messy-string
    # parsing (handles "S (HIGH LEVEL)", "R   Escherichia coli", keyword
    # scans for RESISTANT/SUSCEPTIBLE/SENSITIVE/INTERMEDIATE, etc.) -- the
    # SAME extraction logic used everywhere else in anumaan, instead of a
    # separate, cruder exact-match-only implementation kept just here.
    # Exact-match is tried FIRST so a caller's custom vocabulary (e.g.
    # resistant_values including a site-specific code prep_clean_ast_values
    # would never recognize) always wins; this is purely a fallback for
    # values neither side already resolved.
    needs_fallback <- !is_blank & is.na(result)
    if (any(needs_fallback)) {
      cleaned <- suppressMessages(prep_clean_ast_values(
        data.frame(.v = x_chr[needs_fallback], stringsAsFactors = FALSE),
        value_col = ".v", strict = FALSE
      ))$.v
      result[needs_fallback] <- cleaned
    }

    is_invalid <- !is_blank & is.na(result)
    if (any(is_invalid)) {
      if (unknown_action == "error") {
        stop(sprintf(
          "[.profile_tested_long_wide] Unrecognized AST value(s) -- not in resistant_values/susceptible_values/intermediate_values/untested_values, and not parseable by prep_clean_ast_values(): %s. Pass unknown_action = \"treat_as_untested\" only if that is a deliberate, declared policy.",
          paste(unique(x_chr[is_invalid]), collapse = ", ")
        ), call. = FALSE)
      }
    }
    list(result = result, is_invalid = is_invalid)
  }

  if (is_wide) {
    long <- do.call(rbind, lapply(class_cols, function(cl) {
      cls <- .classify(data[[cl]])
      data.frame(
        .group = data[[group_col]],
        .class = cl,
        .result = cls$result,
        .invalid = cls$is_invalid,
        .event = if (!is.null(event_id_col) && event_id_col %in% names(data)) {
          data[[event_id_col]]
        } else {
          seq_len(nrow(data))
        },
        stringsAsFactors = FALSE
      )
    }))
  } else {
    keep <- data[[class_col]] %in% class_cols
    cls <- .classify(data[[value_col]][keep])
    long <- data.frame(
      .group = data[[group_col]][keep],
      .class = as.character(data[[class_col]][keep]),
      .result = cls$result,
      .invalid = cls$is_invalid,
      .event = if (!is.null(event_id_col) && event_id_col %in% names(data)) {
        data[[event_id_col]][keep]
      } else {
        seq_len(sum(keep))
      },
      stringsAsFactors = FALSE
    )
  }

  tested <- long[!is.na(long$.result), , drop = FALSE] # R, I, or S all count as tested
  all_groups <- unique(long$.group)

  # ---- duplicate / conflict resolution, BEFORE the wide pivot ----
  # reshape() requires (.group, .event, .class) to be unique; silently
  # taking the first duplicate row would let a conflicting R/S pair vanish
  # without a trace (the AST-policy review's point), so this delegates to
  # the shared prep_deduplicate_ast() policy instead of a separate,
  # parallel implementation: group_col = ".group" is the key addition that
  # makes it safe to reuse here -- without it, a patient/event id reused
  # across two different hospitals would read as one conflicting record.
  dedup <- tryCatch(
    suppressMessages(prep_deduplicate_ast(
      tested,
      mode = "remove", strategy = conflict_action,
      group_col = ".group", patient_col = ".event", organism_col = NULL,
      antibiotic_col = ".class", date_col = NULL, ast_col = ".result",
      return_diagnostics = TRUE
    )),
    error = function(e) {
      # prep_deduplicate_ast()'s error names its OWN (internal, dot-prefixed)
      # column names -- translate those back to the caller-facing
      # group_col/event_id_col/class_col names before re-raising, so this
      # delegation stays invisible from the outside.
      msg <- conditionMessage(e)
      msg <- gsub(".group", group_col, msg, fixed = TRUE)
      msg <- gsub(".event", event_id_col %||% "event", msg, fixed = TRUE)
      msg <- gsub(".class", "class", msg, fixed = TRUE)
      msg <- sub("[prep_deduplicate_ast]", "[.profile_tested_long_wide]", msg, fixed = TRUE)
      msg <- gsub("strategy = ", "conflict_action = ", msg, fixed = TRUE)
      stop(msg, call. = FALSE)
    }
  )
  tested <- dedup$data
  n_duplicates_collapsed <- nrow(dedup$exact_duplicates_removed)
  dup_log <- rbind(
    if (nrow(dedup$exact_duplicates_removed) > 0L) {
      as.data.frame(do.call(rbind, lapply(split(dedup$exact_duplicates_removed, list(dedup$exact_duplicates_removed$.group, dedup$exact_duplicates_removed$.class), drop = TRUE), function(z) {
        data.frame(group = z$.group[1], class = z$.class[1], type = "collapsed", n = nrow(z))
      })))
    } else {
      data.frame(group = character(0), class = character(0), type = character(0), n = integer(0))
    },
    if (nrow(dedup$conflict_groups) > 0L) {
      data.frame(group = dedup$conflict_groups$.group, class = dedup$conflict_groups$.class, type = "conflict", n = dedup$conflict_groups$n_rows)
    } else {
      data.frame(group = character(0), class = character(0), type = character(0), n = integer(0))
    }
  )

  # One group x event x class -> wide pivot, reused by every caller for
  # pairwise/higher-order/complete-panel support. `.result` is the value
  # column being spread across per-class columns named ".result.<class>";
  # omitting it here (as an earlier version of profile_audit_support did)
  # leaves nothing to pivot and every downstream count silently comes back
  # zero.
  wide_tested <- reshape(
    tested[c(".group", ".event", ".class", ".result")],
    idvar = c(".group", ".event"), timevar = ".class", direction = "wide"
  )

  wide_folded <- NULL
  if (!is.null(intermediate_as_resistant)) {
    tested$.result_folded <- tested$.result
    if (isTRUE(intermediate_as_resistant)) {
      # Delegate the I-fold to the same Colistin-aware clinical convention
      # used elsewhere in anumaan, instead of a separate, cruder symmetric
      # I -> R rule maintained only here.
      tested <- suppressMessages(prep_recode_intermediate_ast(
        tested, antibiotic_col = ".class", value_col = ".result_folded",
        colistin_to_s = TRUE, others_to_r = TRUE
      ))
    }
    wide_folded <- reshape(
      tested[c(".group", ".event", ".class", ".result_folded")],
      idvar = c(".group", ".event"), timevar = ".class", direction = "wide"
    )
  }

  list(
    long = long, tested = tested, wide_tested = wide_tested, all_groups = all_groups,
    result_col = function(cl) paste0(".result.", cl),
    wide_folded = wide_folded, result_col_folded = function(cl) paste0(".result_folded.", cl),
    n_duplicates_collapsed = n_duplicates_collapsed, conflict_cells = dedup$conflict_groups, dup_log = dup_log
  )
}


#' Tabulate Class-Level Testing Support, by Group
#'
#' Pure counting -- no statistics. For every group (hospital, district,
#' whatever \code{group_col} holds) and every class in \code{class_cols},
#' counts how many events were tested, how many were resistant, and how many
#' were susceptible. Also produces the same counts for every class pair
#' (co-tested) and for the full class set at once (complete panel), since a
#' hospital can have plenty of events yet too few co-tested pairs or complete
#' panels for a joint model.
#'
#' @param data Data frame, one row per event x class (long) or one row per
#'   event with one column per class (wide) -- auto-detected from whether
#'   \code{class_cols} resolve to one shared value column or to several
#'   distinct columns.
#' @param group_col Character. Column identifying the stratum (hospital,
#'   district, lab, ...). Not assumed to mean "hospital".
#' @param class_cols Character vector. Antibiotic classes in the panel.
#'   For long data, these are the values \code{class_col} takes. For wide
#'   data, these are column names, each holding an S/I/R-type result.
#' @param class_col Character. For long-format data, the column naming
#'   which class a row's result belongs to. Ignored for wide data (when
#'   \code{class_cols} resolve to real column names in \code{data}).
#' @param value_col Character. Column holding the susceptibility result.
#'   Resistant is matched case-insensitively against \code{resistant_values};
#'   susceptible against \code{susceptible_values}; anything else (including
#'   \code{NA}) counts as untested for that class.
#' @param resistant_values,susceptible_values,intermediate_values Character
#'   vectors of values in \code{value_col} counted as resistant /
#'   susceptible / intermediate. Defaults cover the common `"R"`/`"S"`/`"I"`
#'   and `"Resistant"`/`"Susceptible"`/`"Intermediate"` spellings.
#'   Intermediate results still count toward \code{tested} regardless of
#'   \code{intermediate_as_resistant} -- an ambiguous result is still a
#'   result, not a missing test.
#' @param intermediate_as_resistant Logical. When \code{TRUE} (default,
#'   matching the clinical convention already used by
#'   \code{prep_recode_intermediate_ast()} elsewhere in \code{anumaan}),
#'   intermediate results are folded into \code{resistant} for the
#'   class-level counts and into the "R" side of the pairwise RR/RS/SR/SS
#'   breakdown. When \code{FALSE}, intermediate results are reported only
#'   in their own \code{intermediate} column/count and excluded from the
#'   binary resistant/susceptible tallies (class-level) and from the
#'   pairwise cell breakdown (they still count toward \code{co_tested},
#'   since co-testing is about coverage, not the result).
#' @param event_id_col Character or \code{NULL}. If supplied, used to dedupe
#'   rows before counting complete-panel support (one row per event x class
#'   is assumed otherwise).
#'
#' @param unknown_action,conflict_action Passed through to the shared AST
#'   policy in \code{\link{.profile_tested_long_wide}}: \code{unknown_action}
#'   (default \code{"error"}) controls what happens to a non-blank
#'   \code{value_col} value that matches none of \code{resistant_values}/
#'   \code{susceptible_values}/\code{intermediate_values}; \code{conflict_action}
#'   (default \code{"error"}) controls what happens when the same group x
#'   event x class cell has duplicate rows with DIFFERING results.
#'
#' @return A list with three data frames:
#'   \describe{
#'     \item{\code{class}}{one row per group x class: \code{tested},
#'       \code{resistant}, \code{susceptible}, \code{intermediate} (the raw
#'       intermediate count, always reported regardless of how it was
#'       folded into \code{resistant}) -- kept for backward compatibility --
#'       plus the separated AST-policy counts \code{n_raw_tested} (every
#'       non-blank raw value seen, including unrecognized ones, before
#'       dedup), \code{n_binary_observed} (= \code{tested}, the resolved
#'       R/S/I count after dedup/conflict handling), \code{n_R}/\code{n_S}/
#'       \code{n_I} (raw per-category counts, before
#'       \code{intermediate_as_resistant} folding), \code{n_untested}
#'       (genuinely blank/\code{NA}/declared untested-token rows),
#'       \code{n_invalid} (non-blank, unrecognized values -- only nonzero
#'       when \code{unknown_action = "treat_as_untested"}; otherwise such a
#'       value stops the function), and \code{n_duplicates_collapsed}
#'       (identical-result duplicate rows collapsed to one)}
#'     \item{\code{pair}}{one row per group x unordered class pair:
#'       \code{co_tested} (both classes have any result, R/I/S), and the
#'       \code{n_RR}/\code{n_RS}/\code{n_SR}/\code{n_SS} 2x2 breakdown plus
#'       \code{near_separation} (\code{TRUE} when any one of those four
#'       cells is below 5, the same historical near-separation rule used
#'       for the J2-L hospital tiers)}
#'     \item{\code{profile}}{one row per group: \code{complete_panel}, the
#'       count of events with every class in \code{class_cols} tested}
#'   }
#' @export
profile_audit_support <- function(data,
                             group_col,
                             class_cols,
                             class_col = NULL,
                             value_col = "antibiotic_value",
                             resistant_values = c("R", "Resistant", "resistant"),
                             susceptible_values = c("S", "Susceptible", "susceptible"),
                             intermediate_values = c("I", "Intermediate", "intermediate"),
                             intermediate_as_resistant = TRUE,
                             event_id_col = NULL,
                             unknown_action = c("error", "treat_as_untested"),
                             conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins", "first")) {
  unknown_action <- match.arg(unknown_action)
  conflict_action <- match.arg(conflict_action)
  piv <- .profile_tested_long_wide(
    data, group_col, class_cols, class_col, value_col,
    resistant_values, susceptible_values, intermediate_values, event_id_col,
    unknown_action = unknown_action, conflict_action = conflict_action,
    intermediate_as_resistant = intermediate_as_resistant
  )
  tested <- piv$tested
  wide_tested <- piv$wide_tested
  wide_folded <- piv$wide_folded
  all_groups <- piv$all_groups
  result_col <- piv$result_col
  result_col_folded <- piv$result_col_folded

  # ---- class-level support ----
  # resistant/susceptible are read off `.result_folded`, which -- when
  # intermediate_as_resistant = TRUE -- was already resolved by
  # prep_recode_intermediate_ast()'s Colistin-aware rule (Colistin I -> S,
  # every other class I -> R), not a separate symmetric I -> R rule kept
  # only here. `intermediate` stays the RAW I count regardless, so it's
  # always reported however it was folded (unchanged contract).
  class_support <- do.call(rbind, lapply(split(tested, list(tested$.group, tested$.class), drop = TRUE), function(z) {
    n_R <- sum(z$.result == "R")
    n_I <- sum(z$.result == "I")
    n_S <- sum(z$.result == "S")
    data.frame(
      group = z$.group[1], class = z$.class[1],
      tested = nrow(z),
      resistant = sum(z$.result_folded == "R"),
      susceptible = sum(z$.result_folded == "S"),
      intermediate = n_I,
      n_binary_observed = nrow(z),
      n_R = n_R, n_S = n_S, n_I = n_I,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(class_support)) {
    # Can legitimately happen (e.g. every row for every group x class was
    # dropped as a conflicting duplicate): still return the full zero-filled
    # grid below, not a merge() crash on a NULL right-hand side.
    class_support <- data.frame(
      group = character(0), class = character(0), tested = integer(0),
      resistant = integer(0), susceptible = integer(0), intermediate = integer(0),
      n_binary_observed = integer(0), n_R = integer(0), n_S = integer(0), n_I = integer(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(class_support) <- NULL

  # ---- raw-tested / untested / invalid counts, from the UNDEDUPED long
  # table (every raw row, before the dedup/conflict resolution above) ----
  raw_counts <- do.call(rbind, lapply(split(piv$long, list(piv$long$.group, piv$long$.class), drop = TRUE), function(z) {
    data.frame(
      group = z$.group[1], class = z$.class[1],
      n_raw_tested = sum(!is.na(z$.result) | z$.invalid),
      n_untested = sum(is.na(z$.result) & !z$.invalid),
      n_invalid = sum(z$.invalid),
      stringsAsFactors = FALSE
    )
  }))
  rownames(raw_counts) <- NULL

  dup_counts <- if (nrow(piv$dup_log) > 0L) {
    do.call(rbind, lapply(split(piv$dup_log, list(piv$dup_log$group, piv$dup_log$class), drop = TRUE), function(z) {
      data.frame(
        group = z$group[1], class = z$class[1],
        n_duplicates_collapsed = sum(z$n[z$type == "collapsed"]),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(group = character(0), class = character(0), n_duplicates_collapsed = integer(0))
  }

  # ensure every group x class combination appears, zero-filled, so a class
  # with NO testing anywhere in a group isn't silently absent from the table
  full_grid <- expand.grid(group = all_groups, class = class_cols, stringsAsFactors = FALSE)
  class_support <- merge(full_grid, class_support, by = c("group", "class"), all.x = TRUE)
  class_support <- merge(class_support, raw_counts, by = c("group", "class"), all.x = TRUE)
  class_support <- merge(class_support, dup_counts, by = c("group", "class"), all.x = TRUE)
  count_cols <- c("tested", "resistant", "susceptible", "intermediate",
                   "n_binary_observed", "n_R", "n_S", "n_I",
                   "n_raw_tested", "n_untested", "n_invalid", "n_duplicates_collapsed")
  class_support[count_cols][is.na(class_support[count_cols])] <- 0L

  # ---- pairwise co-testing support ----
  # Binarize a FOLDED result column to 1 (R)/0 (S)/NA (untested, or I when
  # not folded into resistant) -- kept separate from "co_tested" above,
  # since co-testing is about coverage (any result at all) while the
  # RR/RS/SR/SS breakdown needs a clean binary outcome on both sides.
  # `.result_folded` already applied the Colistin-aware fold (when
  # intermediate_as_resistant = TRUE), so no separate fold logic belongs
  # here.
  .binarize <- function(x) {
    out <- rep(NA_integer_, length(x))
    out[x == "R"] <- 1L
    out[x == "S"] <- 0L
    out
  }

  if (length(class_cols) >= 2L) {
    pairs <- utils::combn(class_cols, 2, simplify = FALSE)
    pair_rows <- list()
    for (p in pairs) {
      c1 <- result_col_folded(p[1]); c2 <- result_col_folded(p[2])
      have_cols <- c1 %in% names(wide_folded) && c2 %in% names(wide_folded)
      both <- if (have_cols) {
        !is.na(wide_folded[[c1]]) & !is.na(wide_folded[[c2]])
      } else {
        rep(FALSE, nrow(wide_folded))
      }
      b1 <- if (have_cols) .binarize(wide_folded[[c1]]) else rep(NA_integer_, nrow(wide_folded))
      b2 <- if (have_cols) .binarize(wide_folded[[c2]]) else rep(NA_integer_, nrow(wide_folded))
      resolved <- !is.na(b1) & !is.na(b2)  # excludes I-when-not-folded from the 2x2 breakdown
      for (g in all_groups) {
        gg <- wide_folded$.group == g
        n_RR <- sum(resolved & gg & b1 == 1L & b2 == 1L)
        n_RS <- sum(resolved & gg & b1 == 1L & b2 == 0L)
        n_SR <- sum(resolved & gg & b1 == 0L & b2 == 1L)
        n_SS <- sum(resolved & gg & b1 == 0L & b2 == 0L)
        n_resolved <- n_RR + n_RS + n_SR + n_SS
        pair_rows[[length(pair_rows) + 1]] <- data.frame(
          group = g, class_a = p[1], class_b = p[2],
          co_tested = sum(both & gg),
          n_RR = n_RR, n_RS = n_RS, n_SR = n_SR, n_SS = n_SS,
          near_separation = n_resolved > 0L && min(n_RR, n_RS, n_SR, n_SS) < 5L,
          stringsAsFactors = FALSE
        )
      }
    }
    pair_support <- do.call(rbind, pair_rows)
    rownames(pair_support) <- NULL
  } else {
    pair_support <- data.frame(group = character(0), class_a = character(0),
                               class_b = character(0), co_tested = integer(0),
                               n_RR = integer(0), n_RS = integer(0), n_SR = integer(0), n_SS = integer(0),
                               near_separation = logical(0))
  }

  # ---- complete-panel support ----
  class_cols_present <- intersect(vapply(class_cols, result_col, character(1)), names(wide_tested))
  complete <- if (length(class_cols_present) == length(class_cols)) {
    rowSums(!is.na(wide_tested[class_cols_present])) == length(class_cols)
  } else {
    rep(FALSE, nrow(wide_tested))
  }
  profile_support <- do.call(rbind, lapply(all_groups, function(g) {
    data.frame(group = g, complete_panel = sum(complete & wide_tested$.group == g), stringsAsFactors = FALSE)
  }))
  rownames(profile_support) <- NULL

  list(class = class_support, pair = pair_support, profile = profile_support)
}


#' Evaluate Precision Against a Set of Named Tiers
#'
#' Takes the raw counts from \code{\link{profile_audit_support}} and computes, for
#' each row, how wide a confidence interval a reported proportion from that
#' count would actually have -- then marks whether that width clears each
#' requested precision tier. Class-level rows use the Wilson score
#' interval (single proportion); pairwise and profile rows use the
#' conservative simultaneous Goodman interval (several proportions reported
#' from one group at once), matching the same two interval families used in
#' \code{methodology/j2_support_precision/} in anumaan-analysis.
#'
#' @param support One element of the list returned by
#'   \code{\link{profile_audit_support}} (i.e. \code{$class}, \code{$pair}, or
#'   \code{$profile}), or the whole list (all three are evaluated).
#' @param tiers Named numeric vector, half-width in percentage points for
#'   each named tier, ordered from most to least demanding. Default
#'   \code{c(tier1 = 5, tier2 = 7.5, tier3 = 10)} -- the same three
#'   percentage-point bars used historically for J2-L, just named plainly
#'   by rank (\code{tier1} = narrowest interval required = hardest to
#'   clear) instead of labels like "strict"/"moderate"/"loose" that imply
#'   an interpretation this function doesn't make. Any number of tiers is
#'   fine -- add a \code{tier4}, rename them, whatever a given use needs.
#' @param conf_level Numeric. Confidence level for the interval. Default
#'   0.95.
#' @param n_categories Integer. Number of simultaneous categories for the
#'   Goodman interval on pair/profile rows (e.g. 4 for the RR/RS/SR/SS cells
#'   of one pair, or \code{2^D} for a full D-class profile). Required when
#'   evaluating \code{$pair}/\code{$profile} directly; ignored for
#'   \code{$class}.
#'
#' @return The same data frame(s) as the input, with an added
#'   \code{interval_half_width} column (percentage points, at the
#'   worst-case/conservative proportion for Goodman rows) and one logical
#'   \code{passes_<tier>} column per entry in \code{tiers} (e.g.
#'   \code{passes_tier1}, \code{passes_tier2}, ...).
#' @export
profile_plan_precision <- function(support,
                                    tiers = c(tier1 = 5, tier2 = 7.5, tier3 = 10),
                                    conf_level = 0.95,
                                    n_categories = NULL) {
  if (is.list(support) && all(c("class", "pair", "profile") %in% names(support))) {
    # `$pair` is ALWAYS exactly 4 categories (the RR/RS/SR/SS 2x2 breakdown),
    # by construction, regardless of panel size D -- hardcoded, never
    # derived from the caller's `n_categories` (which is for `$profile`,
    # a `2^D`-category table whose category count DOES depend on D). A
    # previous version used `n_categories %||% 4L` here, which meant
    # passing n_categories for `$profile` silently also overrode `$pair`'s
    # category count to the same (wrong) value.
    return(list(
      class = profile_plan_precision(support$class, tiers, conf_level),
      pair = profile_plan_precision(support$pair, tiers, conf_level, n_categories = 4L),
      profile = profile_plan_precision(support$profile, tiers, conf_level,
                                        n_categories = n_categories)
    ))
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  .wilson_half_width <- function(n, p = 0.5) {
    if (n <= 0) return(NA_real_)
    denom <- 1 + z^2 / n
    center_adj <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
    100 * center_adj
  }

  # Conservative simultaneous Goodman (1965) half-width for k categories out
  # of n trials: a = qchisq(1 - alpha/k, df=1) (Bonferroni-corrected across
  # the k categories), half-width(p) = sqrt(a^2 + 4*a*n*p*(1-p)) / (2*(n+a)).
  # Evaluated at p = 1/2, the TRUE least-favorable (width-maximizing)
  # proportion -- p*(1-p) is a parabola maximized at p=1/2 for any n/k/a, not
  # at p=1/k (an earlier version of this function used p=1/k AND a different,
  # non-standard formula shape; both were wrong -- verified against an
  # independent reference calculation, see profile_plan_precision()'s tests).
  .goodman_half_width <- function(n, k) {
    if (n <= 0 || is.na(k) || k < 2) return(NA_real_)
    alpha <- 1 - conf_level
    a <- stats::qchisq(1 - alpha / k, df = 1)
    p <- 0.5
    100 * sqrt(a^2 + 4 * a * n * p * (1 - p)) / (2 * (n + a))
  }

  if ("tested" %in% names(support)) {
    support$interval_half_width <- vapply(support$tested, .wilson_half_width, numeric(1))
  } else if ("co_tested" %in% names(support)) {
    k <- n_categories %||% 4L
    support$interval_half_width <- vapply(support$co_tested, .goodman_half_width, numeric(1), k = k)
  } else if ("complete_panel" %in% names(support)) {
    if (is.null(n_categories)) {
      stop("n_categories is required to evaluate profile-level (complete_panel) rows -- ",
           "pass 2^D for a D-class panel.")
    }
    support$interval_half_width <- vapply(support$complete_panel, .goodman_half_width, numeric(1), k = n_categories)
  } else {
    stop("support must be a data frame from profile_audit_support() ($class, $pair, or $profile)")
  }

  for (nm in names(tiers)) {
    support[[paste0("passes_", nm)]] <- !is.na(support$interval_half_width) &
      support$interval_half_width <= tiers[[nm]]
  }
  support
}


#' Build an Eligible-Group Roster From Precision-Evaluated Support
#'
#' The function that was missing from the J2-L history: takes the pass/fail
#' table from \code{\link{profile_plan_precision}} and actually decides
#' which groups (hospitals, ...) are in, producing both the roster vector
#' and a traceable manifest -- replacing a hand-typed list with a real,
#' reproducible, data-derived decision.
#'
#' @param thresholds Output of \code{\link{profile_plan_precision}} (the
#'   list form, with \code{$class}/\code{$pair}/\code{$profile}).
#' @param min_tier Character. Which \code{passes_<tier>} column must be
#'   \code{TRUE} -- one of the names used in \code{tiers} when
#'   \code{\link{profile_plan_precision}} was run (e.g. \code{"tier1"},
#'   \code{"tier2"}, ...). Default \code{"tier2"} (the middle of the
#'   default three-tier set, \code{tier1=5/tier2=7.5/tier3=10} points).
#' @param estimands Character vector, any of \code{"class"}, \code{"pair"},
#'   \code{"profile"}. A group must pass \code{min_tier} on every requested
#'   estimand, for every class/pair/the profile, to be included. Default
#'   all three (the strictest, most defensible choice -- narrow this
#'   explicitly if e.g. only class margins matter for a given use).
#' @param mode Character, \code{"all"} or \code{"complete_case"}. Recorded
#'   in the manifest only -- this function does not itself subset event-level
#'   data; it documents which cohort definition the roster is meant for.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{roster}}{character vector of included group names}
#'     \item{\code{manifest}}{data frame, one row per group, with a logical
#'       \code{included} column and a \code{reason} string naming the first
#'       failing estimand/class when excluded}
#'   }
#' @export
profile_define_cohort <- function(thresholds,
                            min_tier = "tier2",
                            estimands = c("class", "pair", "profile"),
                            mode = c("all", "complete_case")) {
  mode <- match.arg(mode)
  pass_col <- paste0("passes_", min_tier)

  all_groups <- unique(unlist(lapply(estimands, function(e) thresholds[[e]]$group)))

  manifest <- do.call(rbind, lapply(all_groups, function(g) {
    failures <- character(0)
    for (e in estimands) {
      tbl <- thresholds[[e]]
      if (!pass_col %in% names(tbl)) {
        stop(sprintf("'%s' not found -- was profile_plan_precision() run with this tier?", pass_col))
      }
      rows <- tbl[tbl$group == g, , drop = FALSE]
      bad <- rows[!isTRUE(all(rows[[pass_col]])) & !rows[[pass_col]], , drop = FALSE]
      if (nrow(bad) > 0L) {
        label <- if ("class" %in% names(bad)) bad$class[1]
        else if ("class_a" %in% names(bad)) paste(bad$class_a[1], bad$class_b[1], sep = ":")
        else "profile"
        failures <- c(failures, sprintf(
          "%s estimand '%s': interval half-width %.1f exceeds the %s bar",
          e, label, bad$interval_half_width[1], min_tier
        ))
      }
    }
    data.frame(
      group = g,
      included = length(failures) == 0L,
      reason = if (length(failures)) failures[1] else sprintf("clears %s on all requested estimands", min_tier),
      mode = mode,
      stringsAsFactors = FALSE
    )
  }))
  rownames(manifest) <- NULL

  list(roster = manifest$group[manifest$included], manifest = manifest)
}


#' Minimum Complete-Event Count for a Simultaneous Profile-Precision Bound
#'
#' A distribution-free sample-size benchmark for estimating ALL \code{K = 2^D}
#' joint-profile probabilities of a \code{D}-class panel at once, each within
#' absolute error \code{epsilon}, from independent complete events. Derived
#' from Hoeffding's inequality applied to each of the \code{K} profile
#' indicators plus a union bound across profiles (optionally also across
#' \code{n_strata} strata, e.g. hospitals, via a further Bonferroni
#' correction):
#' \deqn{n_{\min} = \lceil \log(2 K S / \alpha) / (2 \epsilon^2) \rceil}
#' where \eqn{S} = \code{n_strata} (or 1 if \code{NULL}).
#'
#' This is a different quantity from the Wilson/Goodman interval widths in
#' \code{\link{profile_plan_precision}}: Wilson/Goodman answer "how wide is
#' the interval around one estimated proportion", while this answers "how
#' many complete events guarantee every one of the \code{K} profile estimates
#' is simultaneously within \code{epsilon} of truth, with probability
#' \code{1 - alpha}". Neither subsumes the other; both are reported.
#'
#' The bound is deliberately conservative (distribution-free, no credit for
#' partial events, no correction for clustering or non-representative
#' complete-case selection) and assumes independent, representative complete
#' events -- stated limitations, not footnotes to ignore. It assesses direct
#' data support only; a stratum below the bound may still be modelable via
#' partial events and shared cross-class structure, but that requires a
#' separate model-based assessment (\code{profile_sensitivity}), not this
#' function.
#'
#' @param D Integer. Number of antibiotic classes in the candidate panel.
#' @param epsilon Numeric in \code{(0,1)}. Maximum acceptable absolute error
#'   per profile probability -- a scientific/clinical choice, not something
#'   this function derives.
#' @param alpha Numeric in \code{(0,1)}. Significance level (so
#'   \code{1 - alpha} is the per-stratum confidence). Default \code{0.05}.
#' @param n_strata Integer or \code{NULL}. When supplied, applies a
#'   Bonferroni correction across this many strata (e.g. hospitals) so the
#'   \code{1 - alpha} guarantee holds jointly across all of them, not just
#'   for each one individually. \code{NULL} (default) gives a per-stratum
#'   guarantee only -- the usual case when each hospital is checked on its
#'   own against a shared benchmark.
#'
#' @return A list: \code{D}, \code{K} (\code{= 2^D}), \code{epsilon},
#'   \code{alpha}, \code{n_strata}, \code{n_min}.
#' @export
profile_precision_bound <- function(D, epsilon, alpha = 0.05, n_strata = NULL) {
  if (D < 1L) stop("D must be >= 1")
  if (epsilon <= 0 || epsilon >= 1) stop("epsilon must be in (0, 1)")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")

  K <- 2^D
  S <- if (is.null(n_strata)) 1 else n_strata
  n_min <- ceiling(log(2 * K * S / alpha) / (2 * epsilon^2))

  list(D = D, K = K, epsilon = epsilon, alpha = alpha, n_strata = n_strata, n_min = n_min)
}


#' Pathogen x Candidate-Panel x Hospital Feasibility Table
#'
#' Builds the feasibility table the way \code{\link{profile_precision_bound}}
#' says it has to be built: the \code{n_min} benchmark is computed ONCE per
#' candidate panel (shared across strata), and checked against each
#' stratum's OWN complete-event count for that exact panel. There is no
#' aggregation step, no coverage-fraction cutoff, and no auto-selected "best"
#' panel size -- this reports \code{meets_benchmark} per
#' pathogen-panel-stratum row and leaves the choice of how many classes to
#' use, per pathogen, to be made by inspecting that table (together with
#' \code{\link{profile_plan_precision}}'s Wilson/Goodman widths and clinical
#' judgment about which classes matter), not by this function.
#'
#' \code{candidate_panels} is supplied explicitly, never auto-generated: a
#' panel of a given size is not interchangeable with any other same-size
#' panel (A-B-C-D and A-B-D-E can have very different joint support even at
#' the same \code{D}), so which combinations to test is a decision for the
#' caller -- a clinically pre-screened shortlist, or an exhaustive
#' enumeration the caller builds and passes in.
#'
#' @param data,group_col,class_col,value_col,resistant_values,susceptible_values,intermediate_values,event_id_col
#'   Passed through to
#'   \code{\link{profile_audit_support}} once per candidate panel (its
#'   \code{class_cols} argument is the candidate panel itself).
#' @param candidate_panels A named list of character vectors, each a
#'   candidate set of classes to evaluate (e.g.
#'   \code{list(core4 = c("Carbapenems","Fluoroquinolones","Penicillins",
#'   "Aminoglycosides"), alt4 = c(...))}). Unnamed entries are labelled
#'   \code{panel_1}, \code{panel_2}, ...
#' @param epsilon,alpha Passed to \code{\link{profile_precision_bound}}.
#' @param simultaneous_across_groups Logical. When \code{TRUE}, applies the
#'   Bonferroni correction across the number of groups found for each panel
#'   (a joint guarantee across all of them). Default \code{FALSE} (each
#'   group checked against its own \code{1 - alpha} guarantee).
#' @param unknown_action,conflict_action Same meaning as in
#'   \code{\link{profile_audit_support}}, passed through once per candidate
#'   panel -- the same AST policy is applied to every panel being screened.
#'
#' @return A data frame, one row per candidate panel x group:
#'   \code{panel}, \code{classes}, \code{D}, \code{K}, \code{group},
#'   \code{n_complete}, \code{n_min}, \code{meets_benchmark}.
#' @export
profile_check_panel_feasibility <- function(data,
                                            group_col,
                                            candidate_panels,
                                            class_col = NULL,
                                            value_col = "antibiotic_value",
                                            resistant_values = c("R", "Resistant", "resistant"),
                                            susceptible_values = c("S", "Susceptible", "susceptible"),
                                            intermediate_values = c("I", "Intermediate", "intermediate"),
                                            event_id_col = NULL,
                                            epsilon,
                                            alpha = 0.05,
                                            simultaneous_across_groups = FALSE,
                                            unknown_action = c("error", "treat_as_untested"),
                                            conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins", "first")) {
  unknown_action <- match.arg(unknown_action)
  conflict_action <- match.arg(conflict_action)
  if (is.null(names(candidate_panels)) || any(names(candidate_panels) == "")) {
    nm <- names(candidate_panels)
    if (is.null(nm)) nm <- rep("", length(candidate_panels))
    nm[nm == ""] <- paste0("panel_", which(nm == ""))
    names(candidate_panels) <- nm
  }

  rows <- lapply(names(candidate_panels), function(panel_name) {
    panel <- candidate_panels[[panel_name]]
    D <- length(panel)

    support <- profile_audit_support(
      data = data, group_col = group_col, class_cols = panel,
      class_col = class_col, value_col = value_col,
      resistant_values = resistant_values, susceptible_values = susceptible_values,
      intermediate_values = intermediate_values, event_id_col = event_id_col,
      unknown_action = unknown_action, conflict_action = conflict_action
    )
    prof <- support$profile
    n_groups <- nrow(prof)

    bound <- profile_precision_bound(
      D = D, epsilon = epsilon, alpha = alpha,
      n_strata = if (isTRUE(simultaneous_across_groups)) n_groups else NULL
    )

    data.frame(
      panel = panel_name,
      classes = paste(panel, collapse = "; "),
      D = D,
      K = bound$K,
      group = prof$group,
      n_complete = prof$complete_panel,
      n_min = bound$n_min,
      meets_benchmark = prof$complete_panel >= bound$n_min,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Audit Direct Co-Tested Support for Each Interaction Order of a Candidate Panel
#'
#' A candidate panel's model includes one interaction term per class subset
#' of order 2 up to \code{max_order} (pairs, triples, ...). This counts,
#' for each such subset and each group, how many events have EVERY class in
#' that subset tested together -- i.e. whether that interaction term has
#' any direct co-tested support at all, as opposed to being estimated
#' almost entirely from the model's shared cross-class prior (\code{eta_sd}
#' in \code{\link{profile_build_design}}).
#'
#' \strong{This is a structural-support audit, not an identifiability
#' check, and does not claim to be one.} \code{n_co_tested >= min_co_tested}
#' establishes neither statistical identification nor adequate information
#' for a subset's interaction term -- it is an operational, configurable
#' warning threshold on raw co-testing counts, nothing more. A formal
#' identifiability assessment would examine a FITTED model's observed-data
#' information matrix (or the rank of the observed-probability map) for
#' unidentified or very weak directions, and whether posterior conclusions
#' are prior-dominated rather than data-driven -- that requires an actual
#' fit, not just the raw co-testing pattern, and is a separate, later
#' concern (closer to \code{profile_sensitivity} once a fit exists, not
#' yet built into this family). This function is the cheap, pre-fit screen
#' only: subsets with zero or sparse direct co-tested support are flagged
#' before any time is spent fitting at all, nothing stronger is claimed.
#'
#' @param data,group_col,class_col,value_col,resistant_values,susceptible_values,intermediate_values,event_id_col
#'   Same as
#'   \code{\link{profile_audit_support}}.
#' @param panel Character vector. The candidate panel's classes.
#' @param max_order Integer or \code{NULL}. Highest interaction order to
#'   check. Default: \code{length(panel)} (every order up to the complete
#'   panel).
#' @param min_co_tested Integer. An operational warning threshold, not an
#'   identifiability criterion: an interaction term is flagged
#'   \code{has_direct_support = FALSE} for a group when fewer than this
#'   many events in that group have every class in the subset tested
#'   together. Default \code{5}, the same historical near-separation
#'   threshold used elsewhere in this family
#'   (\code{\link{profile_audit_support}}'s \code{near_separation} flag).
#'
#' @param unknown_action,conflict_action Same meaning as in
#'   \code{\link{profile_audit_support}} -- the same shared AST policy in
#'   \code{\link{.profile_tested_long_wide}}, applied identically here.
#'
#' @return A data frame, one row per (order, class subset, group):
#'   \code{order}, \code{classes} (the subset, "; "-joined), \code{group},
#'   \code{n_co_tested}, \code{has_direct_support}.
#' @export
profile_audit_interaction_support <- function(data,
                                          group_col,
                                          panel,
                                          class_col = NULL,
                                          value_col = "antibiotic_value",
                                          resistant_values = c("R", "Resistant", "resistant"),
                                          susceptible_values = c("S", "Susceptible", "susceptible"),
                                          intermediate_values = c("I", "Intermediate", "intermediate"),
                                          event_id_col = NULL,
                                          max_order = NULL,
                                          min_co_tested = 5L,
                                          unknown_action = c("error", "treat_as_untested"),
                                          conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins", "first")) {
  unknown_action <- match.arg(unknown_action)
  conflict_action <- match.arg(conflict_action)
  if (is.null(max_order)) max_order <- length(panel)
  max_order <- min(max_order, length(panel))
  if (max_order < 2L) {
    stop("max_order must be >= 2 (a structural-support audit only applies to interaction terms)")
  }

  piv <- .profile_tested_long_wide(
    data, group_col, panel, class_col, value_col,
    resistant_values, susceptible_values, intermediate_values, event_id_col,
    unknown_action = unknown_action, conflict_action = conflict_action
  )
  wide_tested <- piv$wide_tested
  all_groups <- piv$all_groups
  result_col <- piv$result_col

  rows <- list()
  for (k in 2:max_order) {
    subsets <- utils::combn(panel, k, simplify = FALSE)
    for (sub in subsets) {
      cols <- vapply(sub, result_col, character(1))
      present <- cols %in% names(wide_tested)
      co_tested_mask <- if (all(present)) {
        rowSums(!is.na(wide_tested[cols])) == k
      } else {
        rep(FALSE, nrow(wide_tested))
      }
      for (g in all_groups) {
        n <- sum(co_tested_mask & wide_tested$.group == g)
        rows[[length(rows) + 1]] <- data.frame(
          order = k,
          classes = paste(sub, collapse = "; "),
          group = g,
          n_co_tested = n,
          has_direct_support = n >= min_co_tested,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, rows)
}


#' Build the Wide Model Dataset for a Chosen Panel
#'
#' Stage 3 of the resistance-profiling pipeline: turns long- or wide-format
#' AST data, restricted to one chosen panel of antibiotic classes, into the
#' \code{Y} matrix \code{\link{profile_build_design}} expects (events x
#' classes, entries in \code{c(0, 1, NA)}), plus everything needed to trace
#' and model it: group/event identifiers, a covariate table, and an
#' exclusion manifest for events with nothing in the panel tested at all.
#'
#' Covariates are returned as a plain table, not a model-matrix design --
#' identifiers (patient id, hospital, admission id) are carried for
#' grouping and provenance only and do not automatically become predictors
#' or random effects. Build the actual \code{X} passed to
#' \code{\link{profile_build_design}} from this table with
#' \code{model.matrix()} and whatever formula the model specification
#' calls for; that choice belongs to the caller (the analysis repository's
#' model config), not to this function.
#'
#' Every event counted as model-eligible here is the set that
#' \code{\link{profile_audit_support}}/\code{\link{profile_check_panel_feasibility}}
#' should be RE-RUN against once this function has dropped events with
#' nothing tested in the panel at all (or, if a covariate's completeness is
#' later enforced downstream, re-run again on whatever survives that too) --
#' support computed only on the pre-cohort-selection or pre-covariate-
#' reconciliation data describes a bigger cohort than the one that actually
#' reaches the model.
#'
#' @param data,group_col,class_col,value_col,resistant_values,susceptible_values,intermediate_values,event_id_col
#'   Same as
#'   \code{\link{profile_audit_support}}. \code{event_id_col} is required
#'   (not optional) for long-format \code{data} -- a real per-event key is
#'   needed to merge covariates back onto the wide outcome matrix.
#' @param panel Character vector. The chosen panel's classes -- becomes
#'   \code{colnames(Y)}, in this order.
#' @param intermediate_as_resistant Logical. Same meaning as in
#'   \code{\link{profile_audit_support}}; controls whether \code{Y} folds
#'   \code{"I"} into \code{1} (resistant) or leaves it untested
#'   (\code{NA}).
#' @param covariates Character vector of column names in \code{data} to
#'   carry into the returned covariate table, one value per event. Rows
#'   sharing an \code{event_id_col} value are checked for agreement on
#'   every requested covariate first; see \code{covariate_conflict_action}.
#'   Default \code{NULL} (no covariates carried -- just \code{Y} and
#'   identifiers).
#' @param unknown_action,conflict_action Same meaning as in
#'   \code{\link{profile_audit_support}} -- the shared AST policy for
#'   susceptibility RESULTS, applied identically here.
#' @param covariate_conflict_action Character, \code{"error"} (default) or
#'   \code{"first"}. A DIFFERENT conflict from \code{conflict_action} above:
#'   this one is about COVARIATE values (e.g. age, ward) disagreeing across
#'   the several raw rows that share one \code{event_id_col} value, not
#'   about AST results. \code{"error"} stops, naming every event id and
#'   covariate where rows disagree -- silently keeping whichever row
#'   happened to come first is exactly the failure mode this guards
#'   against. \code{"first"} restores the old behavior (take the first
#'   row's value) for callers who have already verified, outside this
#'   function, that their covariates are event-invariant.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{Y}}{matrix, events x \code{length(panel)}, entries in
#'       \code{c(0, 1, NA)}, column names = \code{panel}}
#'     \item{\code{ids}}{data frame, one row per event matching \code{Y}'s
#'       rows: \code{event_id}, \code{group}}
#'     \item{\code{covariates}}{data frame, one row per event matching
#'       \code{Y}'s rows, the requested \code{covariates} columns (zero
#'       columns if none requested)}
#'     \item{\code{status}}{character vector, one per event:
#'       \code{"complete"} (every panel class tested) or \code{"partial"}}
#'     \item{\code{manifest}}{one-row data frame: \code{n_events_input},
#'       \code{n_events_excluded_untested} (no panel class tested at all),
#'       \code{n_events_output} (\code{= nrow(Y)}), \code{n_complete},
#'       \code{n_partial}}
#'   }
#' @export
profile_prepare_model_data <- function(data,
                                       group_col,
                                       panel,
                                       class_col = NULL,
                                       value_col = "antibiotic_value",
                                       resistant_values = c("R", "Resistant", "resistant"),
                                       susceptible_values = c("S", "Susceptible", "susceptible"),
                                       intermediate_values = c("I", "Intermediate", "intermediate"),
                                       intermediate_as_resistant = TRUE,
                                       event_id_col = NULL,
                                       covariates = NULL,
                                       unknown_action = c("error", "treat_as_untested"),
                                       conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins", "first"),
                                       covariate_conflict_action = c("error", "first")) {
  unknown_action <- match.arg(unknown_action)
  conflict_action <- match.arg(conflict_action)
  covariate_conflict_action <- match.arg(covariate_conflict_action)
  is_wide <- all(panel %in% names(data))
  if (!is_wide && is.null(event_id_col)) {
    stop("event_id_col is required for long-format data (needed to merge covariates onto one row per event).")
  }
  if (!is.null(covariates) && length(covariates) > 0L && is.null(event_id_col)) {
    stop("event_id_col is required when covariates is supplied (needed to look up one row per event).")
  }

  piv <- .profile_tested_long_wide(
    data, group_col, panel, class_col, value_col,
    resistant_values, susceptible_values, intermediate_values, event_id_col,
    unknown_action = unknown_action, conflict_action = conflict_action,
    intermediate_as_resistant = intermediate_as_resistant
  )
  wide_tested <- piv$wide_folded
  result_col <- piv$result_col_folded

  for (cl in panel) {
    if (!result_col(cl) %in% names(wide_tested)) wide_tested[[result_col(cl)]] <- NA_character_
  }

  # `.result_folded` already applied the Colistin-aware I-fold (when
  # intermediate_as_resistant = TRUE, delegated to
  # prep_recode_intermediate_ast() inside .profile_tested_long_wide()), so Y
  # only needs the plain R -> 1 / S -> 0 mapping here, not a separate fold.
  Y <- matrix(NA_real_, nrow(wide_tested), length(panel), dimnames = list(NULL, panel))
  for (cl in panel) {
    v <- wide_tested[[result_col(cl)]]
    y <- rep(NA_real_, length(v))
    y[v == "R"] <- 1
    y[v == "S"] <- 0
    Y[, cl] <- y
  }

  n_total <- nrow(Y)
  keep <- rowSums(!is.na(Y)) > 0L
  n_excluded_untested <- sum(!keep)
  Y <- Y[keep, , drop = FALSE]
  wide_tested <- wide_tested[keep, , drop = FALSE]

  ids <- data.frame(event_id = wide_tested$.event, group = wide_tested$.group, stringsAsFactors = FALSE)

  cov_tbl <- if (!is.null(covariates) && length(covariates) > 0L) {
    missing_cov <- setdiff(covariates, names(data))
    if (length(missing_cov) > 0L) {
      stop(sprintf("covariate column(s) not found in data: %s", paste(missing_cov, collapse = ", ")))
    }

    # A real agreement check, not a silent !duplicated() first-row pick: for
    # every event id with MORE than one raw row, every requested covariate
    # must take the same value (NA treated as equal to NA, never silently
    # equal to a non-NA value) across those rows, or this is a genuine data
    # conflict -- the same failure mode as the AST-result conflicts above,
    # just for covariates instead of susceptibility results.
    event_ids_all <- data[[event_id_col]]
    dup_ids <- unique(event_ids_all[duplicated(event_ids_all)])
    if (length(dup_ids) > 0L) {
      conflicts <- character(0)
      for (eid in dup_ids) {
        rows <- data[event_ids_all == eid, covariates, drop = FALSE]
        bad_cov <- covariates[vapply(covariates, function(cv) {
          v <- rows[[cv]]
          !all(is.na(v) == is.na(v[1])) || !all(v[!is.na(v)] == v[!is.na(v)][1])
        }, logical(1))]
        if (length(bad_cov) > 0L) {
          conflicts <- c(conflicts, sprintf("event_id=%s, covariate(s): %s", eid, paste(bad_cov, collapse = ", ")))
        }
      }
      if (length(conflicts) > 0L) {
        if (covariate_conflict_action == "error") {
          stop(sprintf(
            "[profile_prepare_model_data] %d event(s) have covariate values that DISAGREE across their raw rows (never resolved by taking the first row). First conflict: %s. Pass covariate_conflict_action = \"first\" only once you have verified, outside this function, that these covariates really are event-invariant.",
            length(conflicts), conflicts[1]
          ), call. = FALSE)
        }
      }
    }

    first_by_event <- data[!duplicated(data[[event_id_col]]), c(event_id_col, covariates), drop = FALSE]
    names(first_by_event)[1] <- "event_id"
    lookup <- data.frame(event_id = ids$event_id, .order = seq_len(nrow(ids)), stringsAsFactors = FALSE)
    merged <- merge(lookup, first_by_event, by = "event_id", all.x = TRUE, sort = FALSE)
    merged <- merged[order(merged$.order), , drop = FALSE]
    out <- merged[covariates]
    rownames(out) <- NULL
    out
  } else {
    as.data.frame(matrix(nrow = nrow(Y), ncol = 0))
  }

  status <- ifelse(rowSums(!is.na(Y)) == length(panel), "complete", "partial")

  manifest <- data.frame(
    n_events_input = n_total,
    n_events_excluded_untested = n_excluded_untested,
    n_events_output = nrow(Y),
    n_complete = sum(status == "complete"),
    n_partial = sum(status == "partial"),
    stringsAsFactors = FALSE
  )

  list(Y = Y, ids = ids, covariates = cov_tbl, status = status, manifest = manifest)
}


#' Per-Group Cohort Eligibility Summary (One Report, Every Eligibility Fact)
#'
#' Assembles, in one call, the single per-group eligibility report a
#' pipeline review actually wants to read: event/patient counts, the
#' complete/partial/zero-observed breakdown, BOTH independent eligibility
#' verdicts this family computes, and the full class/pair/higher-order
#' support tables -- instead of a caller re-deriving this by hand from five
#' separate calls. This is purely a REPORTING layer: it calls
#' \code{\link{profile_audit_support}},
#' \code{\link{profile_audit_interaction_support}},
#' \code{\link{profile_plan_precision}}, \code{\link{profile_define_cohort}},
#' and \code{\link{profile_precision_bound}} once each, on the SAME
#' already-cohort-filtered \code{data} -- it computes nothing new beyond
#' assembling their outputs into one table.
#'
#' The two eligibility verdicts are genuinely independent, not duplicates of
#' each other -- see \code{\link{profile_precision_bound}}'s docs for why
#' neither subsumes the other:
#' \describe{
#'   \item{\code{direct_support_pass}}{this group's \code{n_complete} meets
#'     \code{\link{profile_precision_bound}}'s distribution-free,
#'     Hoeffding/union-bound \code{n_min} for this exact panel.}
#'   \item{\code{included_in_fit}}{this group clears \code{min_tier} on
#'     every requested estimand, per the Wilson/Goodman precision-tier
#'     roster decision in \code{\link{profile_define_cohort}}.}
#' }
#' A group can pass one and fail the other.
#'
#' \code{n_zero_observed} needs care: the shared pivot this family uses
#' internally only produces a row for events with AT LEAST ONE tested panel
#' class, so a naive count over that pivot would silently miss zero-observed
#' events entirely rather than counting them. This function reconstructs the
#' full event universe from \code{data} itself (via \code{group_col} +
#' \code{event_id_col}) and left-joins the tested-class counts onto it, so
#' an event with nothing tested on this panel is correctly tallied as
#' \code{n_zero_observed}, not silently dropped.
#'
#' @param data Data frame, already restricted to the modeling cohort (e.g.
#'   \code{profile_select_cohort(...)$data}) -- every count here describes
#'   THIS cohort, never the full unfiltered dataset.
#' @param group_col,class_col,value_col,resistant_values,susceptible_values,intermediate_values,intermediate_as_resistant,unknown_action,conflict_action
#'   Same meaning as in \code{\link{profile_audit_support}}.
#' @param panel Character vector. The candidate panel's classes (what
#'   \code{\link{profile_check_panel_feasibility}} calls one entry of
#'   \code{candidate_panels}, and \code{\link{profile_prepare_model_data}}
#'   calls \code{panel}).
#' @param event_id_col Character. Required (no default, unlike elsewhere in
#'   this family) -- event counting is this function's whole job, and a
#'   fallback row-index id would silently miscount events for long-format
#'   data (one row per class per event, not one row per event).
#' @param patient_id_col Character or \code{NULL}. Column identifying the
#'   patient. When supplied, \code{n_patients} is the count of DISTINCT
#'   patients per group; when \code{NULL} (default), \code{n_patients} is
#'   \code{NA_integer_} -- no other function in this family tracks patient
#'   identity, so it is never inferred from something else.
#' @param epsilon Numeric in \code{(0,1)}. Passed to
#'   \code{\link{profile_precision_bound}} for \code{direct_support_pass}.
#' @param alpha Numeric in \code{(0,1)}. Passed to
#'   \code{\link{profile_precision_bound}}. Default \code{0.05}.
#' @param tiers,min_tier,estimands,conf_level Passed to
#'   \code{\link{profile_plan_precision}}/\code{\link{profile_define_cohort}}
#'   for \code{included_in_fit}. Defaults match those functions' own
#'   defaults.
#' @param max_interaction_order Integer or \code{NULL}. Passed to
#'   \code{\link{profile_audit_interaction_support}} as \code{max_order}.
#'   Ignored (that table is returned empty) when \code{length(panel) < 2}.
#' @param min_co_tested Integer. Passed to
#'   \code{\link{profile_audit_interaction_support}}. Default \code{5L}.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{summary}}{one row per group: \code{group},
#'       \code{n_pathogen_events}, \code{n_patients}, \code{n_complete},
#'       \code{n_partial}, \code{n_zero_observed}, \code{direct_support_pass},
#'       \code{included_in_fit}}
#'     \item{\code{class}}{\code{\link{profile_audit_support}}'s
#'       \code{$class} table, on this cohort/panel}
#'     \item{\code{pair}}{\code{\link{profile_audit_support}}'s \code{$pair}
#'       table, on this cohort/panel}
#'     \item{\code{interaction}}{\code{\link{profile_audit_interaction_support}}'s
#'       table (higher-order support); empty data frame if
#'       \code{length(panel) < 2}}
#'   }
#' @export
profile_summarize_cohort_eligibility <- function(data,
                                                  group_col,
                                                  panel,
                                                  class_col = NULL,
                                                  value_col = "antibiotic_value",
                                                  resistant_values = c("R", "Resistant", "resistant"),
                                                  susceptible_values = c("S", "Susceptible", "susceptible"),
                                                  intermediate_values = c("I", "Intermediate", "intermediate"),
                                                  intermediate_as_resistant = TRUE,
                                                  event_id_col,
                                                  patient_id_col = NULL,
                                                  epsilon,
                                                  alpha = 0.05,
                                                  tiers = c(tier1 = 5, tier2 = 7.5, tier3 = 10),
                                                  min_tier = "tier2",
                                                  estimands = c("class", "pair", "profile"),
                                                  conf_level = 0.95,
                                                  max_interaction_order = NULL,
                                                  min_co_tested = 5L,
                                                  unknown_action = c("error", "treat_as_untested"),
                                                  conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins", "first")) {
  unknown_action <- match.arg(unknown_action)
  conflict_action <- match.arg(conflict_action)
  if (!event_id_col %in% names(data)) {
    stop(sprintf("event_id_col '%s' not found in data", event_id_col))
  }
  if (!is.null(patient_id_col) && !patient_id_col %in% names(data)) {
    stop(sprintf("patient_id_col '%s' not found in data", patient_id_col))
  }

  piv <- .profile_tested_long_wide(
    data, group_col, panel, class_col, value_col,
    resistant_values, susceptible_values, intermediate_values, event_id_col,
    unknown_action = unknown_action, conflict_action = conflict_action
  )
  wide_tested <- piv$wide_tested
  result_col <- piv$result_col

  for (cl in panel) {
    if (!result_col(cl) %in% names(wide_tested)) wide_tested[[result_col(cl)]] <- NA_character_
  }
  class_cols_present <- vapply(panel, result_col, character(1))

  # Full event universe from `data` itself -- NOT from wide_tested, which
  # structurally only contains events with >= 1 tested panel class. A
  # zero-observed event never gets a row in wide_tested at all, so it would
  # silently vanish from a count taken directly off that pivot.
  all_event_keys <- unique(data.frame(
    .group = data[[group_col]], .event = data[[event_id_col]], stringsAsFactors = FALSE
  ))
  status_lookup <- data.frame(
    .group = wide_tested$.group, .event = wide_tested$.event,
    n_tested_cols = rowSums(!is.na(wide_tested[class_cols_present])),
    stringsAsFactors = FALSE
  )
  event_status <- merge(all_event_keys, status_lookup, by = c(".group", ".event"), all.x = TRUE)
  event_status$n_tested_cols[is.na(event_status$n_tested_cols)] <- 0L
  event_status$status <- ifelse(
    event_status$n_tested_cols == length(panel), "complete",
    ifelse(event_status$n_tested_cols == 0L, "zero_observed", "partial")
  )

  all_groups <- unique(all_event_keys$.group)
  event_counts <- do.call(rbind, lapply(all_groups, function(g) {
    gg <- event_status$.group == g
    data.frame(
      group = g,
      n_pathogen_events = sum(gg),
      n_complete = sum(gg & event_status$status == "complete"),
      n_partial = sum(gg & event_status$status == "partial"),
      n_zero_observed = sum(gg & event_status$status == "zero_observed"),
      stringsAsFactors = FALSE
    )
  }))

  n_patients_tbl <- if (!is.null(patient_id_col)) {
    do.call(rbind, lapply(all_groups, function(g) {
      data.frame(
        group = g,
        n_patients = length(unique(data[[patient_id_col]][data[[group_col]] == g])),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(group = all_groups, n_patients = NA_integer_, stringsAsFactors = FALSE)
  }

  # ---- class / pair support (reused, not recomputed) ----
  support <- profile_audit_support(
    data, group_col = group_col, class_cols = panel, class_col = class_col, value_col = value_col,
    resistant_values = resistant_values, susceptible_values = susceptible_values,
    intermediate_values = intermediate_values, intermediate_as_resistant = intermediate_as_resistant,
    event_id_col = event_id_col, unknown_action = unknown_action, conflict_action = conflict_action
  )

  # ---- interaction / higher-order support (optional) ----
  interaction_support <- if (length(panel) >= 2L) {
    profile_audit_interaction_support(
      data, group_col = group_col, panel = panel, class_col = class_col, value_col = value_col,
      resistant_values = resistant_values, susceptible_values = susceptible_values,
      intermediate_values = intermediate_values, event_id_col = event_id_col,
      max_order = max_interaction_order, min_co_tested = min_co_tested,
      unknown_action = unknown_action, conflict_action = conflict_action
    )
  } else {
    data.frame(order = integer(0), classes = character(0), group = character(0),
               n_co_tested = integer(0), has_direct_support = logical(0))
  }

  # ---- direct_support_pass: Hoeffding/union-bound n_min vs this panel's n_complete ----
  bound <- profile_precision_bound(D = length(panel), epsilon = epsilon, alpha = alpha)
  direct_support <- data.frame(
    group = event_counts$group,
    direct_support_pass = event_counts$n_complete >= bound$n_min,
    stringsAsFactors = FALSE
  )

  # ---- included_in_fit: Wilson/Goodman precision-tier roster decision ----
  thresholds <- profile_plan_precision(support, tiers = tiers, conf_level = conf_level, n_categories = 2^length(panel))
  cohort <- profile_define_cohort(thresholds, min_tier = min_tier, estimands = estimands)
  included <- data.frame(
    group = cohort$manifest$group,
    included_in_fit = cohort$manifest$included,
    stringsAsFactors = FALSE
  )

  summary_tbl <- merge(event_counts, n_patients_tbl, by = "group", all.x = TRUE)
  summary_tbl <- merge(summary_tbl, direct_support, by = "group", all.x = TRUE)
  summary_tbl <- merge(summary_tbl, included, by = "group", all.x = TRUE)
  summary_tbl <- summary_tbl[c(
    "group", "n_pathogen_events", "n_patients", "n_complete",
    "n_partial", "n_zero_observed", "direct_support_pass", "included_in_fit"
  )]
  rownames(summary_tbl) <- NULL

  list(summary = summary_tbl, class = support$class, pair = support$pair, interaction = interaction_support)
}


`%||%` <- function(a, b) if (is.null(a)) b else a
