# profile_validate_holdout_tasks.R
# New-patient and temporal held-out validation. See profile_fit.R's
# header for how this fits into the pipeline; see profile_validate_masked.R
# for the masked-CELL holdout task this file does NOT duplicate.
#
# Genuinely different holdout mechanics from profile_validate_heldout()
# (which holds out individual CELLS within otherwise-training events):
# these hold out WHOLE EVENTS, predict them with NOTHING known about that
# event at all (the "unconditional" estimand -- information actually
# available at deployment, before any AST result exists for a brand-new
# patient or a not-yet-occurred period), then score whichever classes
# were ACTUALLY tested in the real data for those events.
#
# Leave-one-hospital-out is deliberately NOT implemented here: this
# family's hospital fixed effects are unobserved for an unseen hospital,
# so predicting for one needs an explicitly different model/rule, not a
# default validation task.


#' (Internal) Refit on a Training Subset and Score Whole Held-Out Events
#'
#' Shared machinery for \code{\link{profile_validate_newpatient}} and
#' \code{\link{profile_validate_temporal}}: refits on \code{Y[train_idx,]}/
#' \code{X[train_idx,]} via \code{\link{profile_fit_joint_bernoulli}}, then
#' predicts EVERY held-out event's per-class probabilities using ONLY its
#' covariates (the \code{"unconditional"} estimand via
#' \code{\link{profile_estimate}} -- no AST information from that event is
#' used, matching what is actually available at deployment), and scores
#' against whichever classes were ACTUALLY tested in the real (unmasked)
#' \code{Y} for those events via the SAME scoring machinery
#' \code{\link{profile_validate_masked}} uses. A naturally untested cell
#' contributes no score.
#'
#' @keywords internal
.profile_validate_event_holdout <- function(Y, X, ids, output_dir, holdout_idx,
                                             sampler_config, n_draws, seed, validation_mode) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  N <- nrow(Y)
  train_idx <- setdiff(seq_len(N), holdout_idx)
  if (length(train_idx) == 0L) {
    stop(sprintf("[%s] holdout_idx covers every event -- nothing left to train on.", validation_mode), call. = FALSE)
  }
  if (length(holdout_idx) == 0L) {
    stop(sprintf("[%s] holdout_idx is empty -- nothing to validate.", validation_mode), call. = FALSE)
  }

  design_train <- profile_build_design(Y[train_idx, , drop = FALSE], X[train_idx, , drop = FALSE])
  refit_result <- profile_fit_joint_bernoulli(design = design_train, output_dir = output_dir, sampler_config = sampler_config)

  # B/V are structural (function of D only, not N), so reusing design_train's
  # un-changed; only X (and therefore N) differs for the held-out prediction.
  pred_design <- design_train
  pred_design$stan_data$X <- X[holdout_idx, , drop = FALSE]

  est <- profile_estimate(
    refit_result, pred_design, Y = Y[holdout_idx, , drop = FALSE], ids = ids[holdout_idx, , drop = FALSE],
    estimand = "unconditional", n_draws = n_draws, seed = seed
  )
  B <- design_train$stan_data$B
  D <- design_train$stan_data$D
  class_prob <- vapply(
    seq_len(D), function(d) rowSums(est$profile_probability[, B[, d] == 1, drop = FALSE]),
    numeric(length(holdout_idx))
  )

  classes <- colnames(Y)
  if (is.null(classes)) classes <- paste0("class", seq_len(D))
  rows <- list()
  for (d in seq_len(D)) {
    tested <- which(!is.na(Y[holdout_idx, d]))
    if (length(tested) == 0L) next
    rows[[length(rows) + 1]] <- data.frame(
      group = ids$group[holdout_idx][tested], row = holdout_idx[tested], class = classes[d],
      value = Y[holdout_idx, d][tested], predicted = class_prob[tested, d],
      stringsAsFactors = FALSE
    )
  }
  scored_df <- do.call(rbind, rows)
  if (is.null(scored_df) || nrow(scored_df) == 0L) {
    stop(sprintf(
      "[%s] No held-out event had any genuinely tested class to score against -- naturally untested outcomes have no truth and cannot be used to score.",
      validation_mode
    ), call. = FALSE)
  }

  result <- .profile_score_mask_predictions(
    data.frame(group = scored_df$group, row = scored_df$row, class = scored_df$class, value = scored_df$value,
               stringsAsFactors = FALSE),
    scored_df$predicted, validation_mode
  )
  result$refit_result <- refit_result
  result$n_holdout_events <- length(holdout_idx)
  result$n_train_events <- length(train_idx)
  result$joint_profile_scores <- .profile_joint_scores(est, Y, holdout_idx, ids)
  result
}


#' (Internal) Joint-Profile Log Loss and Multiclass Brier, for Held-Out
#' Events With a Fully Known True Profile
#'
#' Cell-level scores (AUROC, Brier, average precision, calibration -- see
#' \code{\link{.profile_score_mask_predictions}}) can agree perfectly
#' between two models that disagree completely about the joint profile:
#' two models both giving P(classD=R)=P(classC=R)=0.5 can place all of
#' their JOINT mass on entirely different complete states. This scores
#' the WHOLE profile as one unit instead, for exactly the held-out events
#' where the true complete profile is actually known (every class
#' genuinely tested) -- partially-tested held-out events contribute no
#' row here, since there is no single true state to score against.
#'
#' \code{est$profile_probability} (from \code{\link{profile_estimate}},
#' \code{estimand = "unconditional"}) is already the posterior-MEAN
#' probability per event x state -- i.e. already \eqn{\bar q_{e,s}},
#' averaged over draws; no further draw-level averaging happens here.
#'
#' @keywords internal
.profile_joint_scores <- function(est, Y, holdout_idx, ids) {
  D <- ncol(Y)
  state_labels <- est$state_labels
  S <- length(state_labels)
  Y_holdout <- Y[holdout_idx, , drop = FALSE]
  is_complete <- rowSums(!is.na(Y_holdout)) == D
  if (!any(is_complete)) {
    return(data.frame(
      group = character(0), row = integer(0), true_state = character(0),
      joint_log_loss = numeric(0), multiclass_brier = numeric(0), stringsAsFactors = FALSE
    ))
  }
  complete_rows <- which(is_complete)
  true_labels <- apply(Y_holdout[complete_rows, , drop = FALSE], 1L,
                       function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))
  s_true_idx <- match(true_labels, state_labels)
  if (anyNA(s_true_idx)) {
    stop("[.profile_joint_scores] A held-out event's true R/S pattern did not match any known state_labels -- class order mismatch.", call. = FALSE)
  }

  qbar <- est$profile_probability[complete_rows, , drop = FALSE] # n_complete x S
  rows <- lapply(seq_along(complete_rows), function(i) {
    q_true <- qbar[i, s_true_idx[i]]
    onehot <- as.integer(seq_len(S) == s_true_idx[i])
    data.frame(
      group = ids$group[holdout_idx][complete_rows[i]],
      row = holdout_idx[complete_rows[i]],
      true_state = true_labels[i],
      joint_log_loss = -log(max(q_true, 1e-12)),
      multiclass_brier = sum((qbar[i, ] - onehot)^2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' New-Patient Held-Out Validation (Within Represented Hospitals)
#'
#' Holds out EVERY event for a selected set of patients entirely from
#' training (refits on the remaining events via
#' \code{\link{profile_fit_joint_bernoulli}}), then predicts each held-out
#' event using ONLY its covariates -- the \code{"unconditional"} estimand,
#' i.e. nothing known about that patient's AST results, matching what is
#' actually available for a genuinely new patient at deployment -- and
#' scores against whichever classes were actually tested for those
#' events. This is validation WITHIN the hospitals already represented in
#' training; it is NOT leave-one-hospital-out (this family does not
#' implement that as a default task -- a hospital's fixed effect is
#' unobserved for an unseen hospital, needing an explicitly different
#' prediction rule/model).
#'
#' @param Y,X,ids Same as \code{\link{profile_validate_heldout}} -- the
#'   FULL (not yet masked/split) data.
#' @param patient_id Character/integer vector, length \code{nrow(Y)}, in
#'   the SAME row order as \code{Y}/\code{ids} -- the patient identifier
#'   for every event.
#' @param strata Character/integer vector, length \code{nrow(Y)} (same row
#'   order as \code{patient_id}), or \code{NULL}. When supplied and
#'   \code{holdout_patient_ids} is \code{NULL}, \code{holdout_fraction} of
#'   DISTINCT patients is sampled WITHIN EACH STRATUM separately (e.g. one
#'   value per event giving that event's hospital) rather than pooled
#'   across the whole cohort -- avoids a hospital-blind draw leaving a
#'   small hospital with zero (or disproportionately many) held-out
#'   patients. A patient must have the SAME stratum value across all of
#'   their events (checked explicitly; errors otherwise, naming the
#'   patient). Default \code{NULL}: pooled sampling, the original
#'   behaviour, unchanged.
#' @param holdout_patient_ids Character/integer vector or \code{NULL}.
#'   Explicit patient ids to hold out. Default \code{NULL}: sampled via
#'   \code{holdout_fraction} instead.
#' @param holdout_fraction Numeric in (0, 1). Used only when
#'   \code{holdout_patient_ids} is \code{NULL}: this fraction of DISTINCT
#'   patients (not events) is sampled for holdout. Default \code{0.1}.
#' @param output_dir Character. Passed to
#'   \code{\link{profile_fit_joint_bernoulli}} for the training-only refit.
#' @param sampler_config Named list. Passed to
#'   \code{\link{profile_fit_joint_bernoulli}} for the training-only refit.
#' @param n_draws,seed Passed to \code{\link{profile_estimate}} for the
#'   held-out prediction (also the patient-sampling seed, when
#'   \code{holdout_patient_ids} is \code{NULL}).
#'
#' @return Named list: \code{predictions}, \code{summary}, \code{calibration}
#'   (same shape as \code{\link{profile_validate_masked}}),
#'   \code{joint_profile_scores} (data frame, one row per held-out event
#'   whose true complete profile is actually known -- i.e. every class
#'   was genuinely tested for that event: \code{joint_log_loss},
#'   \code{multiclass_brier} scoring the WHOLE profile as one unit, which
#'   cell-level scores in \code{summary} cannot do -- empty with 0 rows
#'   if no held-out event had a complete true profile),
#'   \code{refit_result} (the training-only fit -- check its sampler
#'   diagnostics via \code{\link{profile_check_sampler}} before trusting
#'   these scores), \code{n_holdout_events}, \code{n_train_events},
#'   \code{holdout_patient_ids} (echoed, for provenance).
#' @export
profile_validate_newpatient <- function(Y, X, ids, patient_id, output_dir,
                                        strata = NULL,
                                        holdout_patient_ids = NULL,
                                        holdout_fraction = 0.1,
                                        sampler_config = list(),
                                        n_draws = 500L, seed = 123L) {
  N <- nrow(as.matrix(Y))
  if (length(patient_id) != N) stop("patient_id must have length nrow(Y).")
  if (!is.null(strata) && length(strata) != N) {
    stop("strata must have length nrow(Y) (one value per event, matching patient_id).", call. = FALSE)
  }

  if (is.null(holdout_patient_ids)) {
    .local_seed(as.integer(seed))
    distinct_patients <- unique(patient_id)
    if (is.null(strata)) {
      n_hold <- max(1L, round(holdout_fraction * length(distinct_patients)))
      # sample(x, n) on a length-1 NUMERIC x samples from 1:x instead of
      # returning x itself (R's documented single-number special case) --
      # indexing via sample.int() is safe regardless of patient_id's type
      # or length.
      holdout_patient_ids <- sort(distinct_patients[sample.int(length(distinct_patients), n_hold)])
    } else {
      patient_stratum <- vapply(distinct_patients, function(p) {
        s <- unique(strata[patient_id == p])
        if (length(s) > 1L) {
          stop(sprintf(
            "Patient '%s' has more than one distinct stratum value across their events (%s) -- strata must be constant within a patient.",
            p, paste(s, collapse = ", ")
          ), call. = FALSE)
        }
        as.character(s)
      }, character(1))
      by_stratum <- split(distinct_patients, patient_stratum)
      holdout_patient_ids <- sort(unlist(lapply(names(by_stratum), function(st) {
        ps <- by_stratum[[st]]
        if (length(ps) < 2L) {
          message(sprintf(
            "[profile_validate_newpatient] Stratum '%s' has only %d patient(s) -- holding any out would leave this stratum's fixed effect completely unsupported in training; keeping all of it in training instead (0 held out from this stratum).",
            st, length(ps)
          ))
          return(character(0))
        }
        n_hold <- max(1L, round(holdout_fraction * length(ps)))
        # Same sample.int() indexing fix as above -- a stratum of exactly
        # one numeric patient ID (e.g. 42) would otherwise have sample()
        # draw from 1:42 instead of returning patient 42.
        chosen <- ps[sample.int(length(ps), n_hold)]
        message(sprintf(
          "[profile_validate_newpatient] Stratum '%s': holding out %d of %d patient(s).",
          st, n_hold, length(ps)
        ))
        chosen
      })))
    }
  }
  holdout_idx <- which(patient_id %in% holdout_patient_ids)
  if (length(holdout_idx) == 0L) {
    stop("No events match holdout_patient_ids -- check patient_id/holdout_patient_ids alignment.", call. = FALSE)
  }
  message(sprintf(
    "[profile_validate_newpatient] Holding out %d patient(s), %d event(s); refitting on the remainder...",
    length(unique(holdout_patient_ids)), length(holdout_idx)
  ))

  result <- .profile_validate_event_holdout(
    Y, X, ids, output_dir, holdout_idx, sampler_config, n_draws, seed, "new_patient_holdout"
  )
  result$holdout_patient_ids <- holdout_patient_ids
  result
}


#' Temporal Held-Out Validation
#'
#' Trains on events before \code{cutoff_date} and predicts events on or
#' after it, with NOTHING known about the held-out event's AST results
#' (the \code{"unconditional"} estimand) -- the SAME prediction and
#' scoring mechanics as \code{\link{profile_validate_newpatient}}, just a
#' date-based split instead of a patient-based one. Reports the training/
#' holdout patient-overlap count explicitly when \code{patient_id} is
#' supplied, rather than assuming temporal and patient independence.
#'
#' @param Y,X,ids Same as \code{\link{profile_validate_heldout}}.
#' @param event_date Date vector (or coercible via \code{as.Date}), length
#'   \code{nrow(Y)}, in the SAME row order as \code{Y}/\code{ids}.
#' @param cutoff_date Date (or coercible). Events with \code{event_date <
#'   cutoff_date} train; \code{>= cutoff_date} are held out.
#' @param patient_id Character/integer vector or \code{NULL}. When
#'   supplied, the returned \code{patient_overlap} documents how many
#'   holdout-period patients also appear in the training period --
#'   temporal validation is not automatically patient-independent.
#'   Default \code{NULL} (no overlap report).
#' @param output_dir Character. Passed to
#'   \code{\link{profile_fit_joint_bernoulli}} for the training-only refit.
#' @param sampler_config Named list. Passed to
#'   \code{\link{profile_fit_joint_bernoulli}} for the training-only refit.
#' @param n_draws,seed Passed to \code{\link{profile_estimate}}.
#'
#' @return Named list: \code{predictions}, \code{summary}, \code{calibration},
#'   \code{joint_profile_scores} (see \code{\link{profile_validate_newpatient}}
#'   for its exact shape -- identical here), \code{refit_result},
#'   \code{n_holdout_events}, \code{n_train_events},
#'   \code{cutoff_date} (echoed), and \code{patient_overlap} (list:
#'   \code{n_overlapping_patients}, \code{fraction_of_holdout_patients} --
#'   both \code{NA} when \code{patient_id} is not supplied).
#' @export
profile_validate_temporal <- function(Y, X, ids, event_date, cutoff_date,
                                      patient_id = NULL,
                                      output_dir, sampler_config = list(),
                                      n_draws = 500L, seed = 123L) {
  N <- nrow(as.matrix(Y))
  event_date <- as.Date(event_date)
  cutoff_date <- as.Date(cutoff_date)
  if (length(event_date) != N) stop("event_date must have length nrow(Y).")
  if (!is.null(patient_id) && length(patient_id) != N) stop("patient_id must have length nrow(Y).")

  holdout_idx <- which(event_date >= cutoff_date)
  train_idx <- which(event_date < cutoff_date)
  if (length(holdout_idx) == 0L) stop("No events on/after cutoff_date -- nothing to hold out.", call. = FALSE)
  if (length(train_idx) == 0L) stop("No events before cutoff_date -- nothing to train on.", call. = FALSE)
  message(sprintf(
    "[profile_validate_temporal] Training on %d event(s) before %s, holding out %d event(s) on/after it; refitting...",
    length(train_idx), cutoff_date, length(holdout_idx)
  ))

  overlap <- list(n_overlapping_patients = NA_integer_, fraction_of_holdout_patients = NA_real_)
  if (!is.null(patient_id)) {
    train_patients <- unique(patient_id[train_idx])
    holdout_patients <- unique(patient_id[holdout_idx])
    n_overlap <- length(intersect(train_patients, holdout_patients))
    overlap <- list(
      n_overlapping_patients = n_overlap,
      fraction_of_holdout_patients = n_overlap / length(holdout_patients)
    )
    message(sprintf(
      "[profile_validate_temporal] %d of %d holdout-period patient(s) (%.1f%%) also appear in the training period.",
      n_overlap, length(holdout_patients), 100 * overlap$fraction_of_holdout_patients
    ))
  }

  result <- .profile_validate_event_holdout(
    Y, X, ids, output_dir, holdout_idx, sampler_config, n_draws, seed, "temporal_holdout"
  )
  result$cutoff_date <- cutoff_date
  result$patient_overlap <- overlap
  result
}
