# profile_check_assumptions.R
# Structured pre-fit assumption register. See profile_fit.R's header for
# how this fits into the pipeline.
#
# A deliberately NARROW function: it assembles a fixed 12-row catalog
# (one row per assumption question this family's design review requires),
# auto-fills whatever rows are genuinely computable from diagnostics
# ALREADY BUILT elsewhere in this family (duplicate/conflict detection,
# design rank, direct-support counts, prior-predictive concentration),
# and leaves every other row "not_evaluated" -- never silently marked
# "proved" for a question that is, in fact, statistically untestable from
# the observed data alone (the core discipline Section C of the
# implementation spec requires). A human reviewer closes the remaining
# rows explicitly via `reviewer_notes`.


#' (Internal) The Fixed 12-Row Pre-Fit Assumption Catalog
#' @keywords internal
.profile_assumption_catalog <- function() {
  data.frame(
    assumption_id = c(
      "correct_unit", "outcome_well_defined", "observations_comparable",
      "repeated_events_independent", "design_estimable", "panel_clinically_meaningful",
      "direct_support_present", "testing_selective", "common_interactions_plausible",
      "priors_appropriate", "targets_within_support", "computation_fits_resources"
    ),
    claim = c(
      "Hospital/patient/admission/event keys are correct, and duplicate/conflicting AST results are resolved explicitly, not silently averaged or first-row-wins.",
      "The binary resistant/susceptible outcome is well defined: raw AST tokens, the intermediate-result policy, and drug-to-class aggregation are explicit and documented.",
      "Events being pooled are comparable across hospital, time, specimen, species and testing method.",
      "Repeated events for the same patient/admission are independent enough for the likelihood's independence assumption.",
      "The design matrix is estimable: full rank, no accidental aliasing, no constant/collinear predictors.",
      "The chosen panel is clinically meaningful, not selected solely by complete-count or interval width.",
      "Each class, pair, requested higher-order subset, and complete state/mask has direct (not just prior-borrowed) support, locally and pooled.",
      "Testing (which classes get tested) is not informatively selective on the unobserved resistance value.",
      "Interactions shared across hospitals/strata are a plausible assumption for this panel and roster.",
      "The declared priors do not induce implausible profile probabilities at the model's actual covariates.",
      "Prediction targets (new events/hospitals) are within the training covariate/hospital support.",
      "The configured computation (N, D, 2^D, parameter count, draw/output size) fits available resources."
    ),
    stringsAsFactors = FALSE
  )
}


#' Structured Pre-Fit Assumption Register
#'
#' Assembles the fixed 12-question assumption register this family's
#' design review requires (Section C of the implementation spec),
#' auto-filling whatever rows are genuinely computable from diagnostics
#' ALREADY produced elsewhere in this package -- never re-deriving them --
#' and leaving every other row \code{status = "not_evaluated"}.
#' \strong{A statistically untestable assumption (e.g. whether testing is
#' informatively selective on the unobserved resistance value, or whether
#' shared interactions are scientifically plausible) is never marked
#' \code{"verified_by_data_contract"} by this function}, no matter what
#' inputs are supplied -- only \code{correct_unit}, \code{design_estimable},
#' and \code{direct_support_present}/\code{priors_appropriate} have any
#' automatic path to that status, and only when the relevant hard check
#' actually passed.
#'
#' Rows this function CAN auto-fill, and from what:
#' \describe{
#'   \item{\code{correct_unit}}{from \code{support} --
#'     \code{\link{profile_audit_support}}'s \code{$class$n_duplicates_collapsed}
#'     (identical duplicates collapsed, never inflating a count). Reported
#'     as \code{"diagnosed"}, not \code{"verified_by_data_contract"} --
#'     \code{profile_audit_support()}'s default \code{conflict_action =
#'     "error"} means a CONFLICTING duplicate would already have stopped
#'     whatever call produced \code{support} (so its mere existence is
#'     good evidence), but this function cannot see which
#'     \code{conflict_action} that call actually used, so it stops short
#'     of "verified".}
#'   \item{\code{design_estimable}}{from \code{fixed_effect_design} --
#'     \code{\link{profile_build_fixed_effect_design}}'s return value;
#'     checks \code{rank == ncol(X)}.}
#'   \item{\code{direct_support_present}}{from \code{support} --
#'     \code{\link{profile_audit_support}}'s return value; flags any
#'     zero-count group x class cell. Always \code{status = "diagnosed"},
#'     never \code{"verified"} -- a nonzero count is not identifiability
#'     or adequate precision (see \code{\link{profile_plan_precision}}).}
#'   \item{\code{priors_appropriate}}{from \code{prior_predictive} --
#'     \code{\link{profile_prior_predictive}}'s return value; flags
#'     extreme (near-0 or near-max) entropy ratios. Always
#'     \code{status = "diagnosed"} -- flags only extreme concentration, a
#'     reviewer still has to judge clinical plausibility of the actual
#'     implied margins/pair probabilities.}
#' }
#' Every other row requires an explicit human judgment, supplied via
#' \code{reviewer_notes}.
#'
#' @param support Output of \code{\link{profile_audit_support}}, or
#'   \code{NULL} to leave \code{direct_support_present} unfilled.
#' @param fixed_effect_design Output of
#'   \code{\link{profile_build_fixed_effect_design}}, or \code{NULL} to
#'   leave \code{design_estimable} unfilled.
#' @param prior_predictive Output of \code{\link{profile_prior_predictive}},
#'   or \code{NULL} to leave \code{priors_appropriate} unfilled.
#' @param reviewer_notes Named list, keyed by \code{assumption_id} (an
#'   error names the offending key if any is unrecognized -- this catches
#'   a typo rather than silently dropping the note). Each element is a
#'   list with any of \code{note} (free text), \code{status} (one of
#'   \code{"verified_by_data_contract"}, \code{"diagnosed"},
#'   \code{"assumed"}, \code{"not_evaluated"}, \code{"violated"}),
#'   \code{sensitivity_plan}. Applied AFTER the automatic fills, so a
#'   reviewer can override an automatic verdict with their own judgment
#'   (e.g. downgrading \code{"verified_by_data_contract"} to
#'   \code{"assumed"} if they have reservations the automatic check
#'   cannot see).
#'
#' @return A data frame, one row per \code{assumption_id}: \code{claim},
#'   \code{evidence}, \code{diagnostic}, \code{limitation},
#'   \code{sensitivity_plan}, \code{severity} (\code{"not_assessed"},
#'   \code{"low"}, \code{"moderate"}, \code{"blocking"}), \code{status},
#'   \code{reviewer_note}.
#' @export
profile_check_assumptions <- function(support = NULL,
                                      fixed_effect_design = NULL,
                                      prior_predictive = NULL,
                                      reviewer_notes = list()) {
  reg <- .profile_assumption_catalog()
  reg$evidence <- NA_character_
  reg$diagnostic <- NA_character_
  reg$limitation <- NA_character_
  reg$sensitivity_plan <- NA_character_
  reg$severity <- "not_assessed"
  reg$status <- "not_evaluated"
  reg$reviewer_note <- NA_character_

  .set <- function(id, status = NA, evidence = NA, diagnostic = NA, limitation = NA, severity = NA) {
    i <- which(reg$assumption_id == id)
    if (!is.na(status)) reg$status[i] <<- status
    if (!is.na(evidence)) reg$evidence[i] <<- evidence
    if (!is.na(diagnostic)) reg$diagnostic[i] <<- diagnostic
    if (!is.na(limitation)) reg$limitation[i] <<- limitation
    if (!is.na(severity)) reg$severity[i] <<- severity
  }

  if (!is.null(support)) {
    n_dup <- sum(support$class$n_duplicates_collapsed, na.rm = TRUE)
    .set("correct_unit",
      status = "diagnosed",
      evidence = sprintf(
        "%d identical-duplicate row(s) collapsed across %d group x class cell(s). profile_audit_support()'s default conflict_action = \"error\" means a CONFLICTING duplicate would already have stopped the call that produced this table.",
        n_dup, nrow(support$class)
      ),
      diagnostic = "profile_audit_support()$class$n_duplicates_collapsed.",
      limitation = "Detects duplicate/conflicting rows at the (group, event, class) key only, and cannot see whether that call actually used conflict_action = \"error\" -- does not itself verify that the keys are semantically correct.",
      severity = if (n_dup > 0L) "moderate" else "low"
    )
  }

  if (!is.null(fixed_effect_design)) {
    ok <- fixed_effect_design$rank == ncol(fixed_effect_design$X)
    .set("design_estimable",
      status = if (ok) "verified_by_data_contract" else "violated",
      evidence = sprintf(
        "Design has %d column(s), rank %d (%s).",
        ncol(fixed_effect_design$X), fixed_effect_design$rank, if (ok) "full rank" else "RANK-DEFICIENT"
      ),
      severity = if (ok) "low" else "blocking"
    )
  }

  if (!is.null(support)) {
    cls <- support$class
    n_zero <- sum(cls$tested == 0L)
    .set("direct_support_present",
      status = "diagnosed",
      evidence = sprintf("%d of %d group x class cell(s) have zero direct support.", n_zero, nrow(cls)),
      diagnostic = "profile_audit_support()$class/$pair/$profile counts.",
      limitation = "A nonzero count establishes SOME direct support, not identifiability or adequate precision -- see profile_plan_precision()/profile_precision_bound() for that.",
      severity = if (n_zero > 0L) "moderate" else "low"
    )
  }

  if (!is.null(prior_predictive)) {
    conc <- prior_predictive$concentration
    flag <- any(conc$entropy_ratio < 0.05) || any(conc$entropy_ratio > 0.98)
    .set("priors_appropriate",
      status = "diagnosed",
      evidence = sprintf("Entropy ratio by group: %s.", paste(sprintf("%s=%.3f", conc$group, conc$entropy_ratio), collapse = "; ")),
      diagnostic = "profile_prior_predictive()$concentration.",
      limitation = "Flags only extreme concentration -- does not itself judge clinical plausibility of the implied margins/pair probabilities; a reviewer must still inspect profile_prior_predictive()'s full output.",
      severity = if (flag) "moderate" else "low"
    )
  }

  for (id in names(reviewer_notes)) {
    i <- which(reg$assumption_id == id)
    if (length(i) == 0L) {
      stop(sprintf(
        "reviewer_notes name '%s' is not a known assumption_id. Valid ids: %s",
        id, paste(reg$assumption_id, collapse = ", ")
      ), call. = FALSE)
    }
    entry <- reviewer_notes[[id]]
    if (!is.null(entry$note)) reg$reviewer_note[i] <- entry$note
    if (!is.null(entry$status)) reg$status[i] <- entry$status
    if (!is.null(entry$sensitivity_plan)) reg$sensitivity_plan[i] <- entry$sensitivity_plan
  }

  reg
}
