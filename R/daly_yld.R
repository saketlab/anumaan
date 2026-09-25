# daly_yld.R
# YLD (Years Lived with Disability) calculation functions for AMR burden estimation

# -- PAF_LOS : population attributable fraction for length of stay ---------------

#' Compute PAF for length of stay per resistance profile
#'
#' Computes the population attributable fraction (PAF) for length of stay
#' from a named list of profile data frames, each containing resistance-profile
#' probabilities and profile-level relative risk of LOS.
#'
#' For each pathogen the formula is:
#'   PAF = sum_d R'_kd * (RR_kd - 1) / (1 + sum_d R'_kd * (RR_kd - 1))
#'
#' where R'_kd is the probability of resistance profile d and RR_kd is the
#' corresponding relative LOS multiplier.
#'
#' @param profiles_with_rr Named list returned by \code{assign_rr_to_profiles()}
#'   or \code{filter_profiles_to_rr_classes()}. Each entry is a data frame of
#'   resistance profiles for one pathogen.
#' @param probability_col Character. Column name for profile probability
#'   (must sum to 1 within each pathogen). Default \code{"probability"}.
#' @param rr_profile_col Character. Column name for profile-level LOS relative
#'   risk. Default \code{"RR_LOS_profile"}.
#' @param profile_col Character. Column name for the profile identifier.
#'   Default \code{"profile"}.
#'
#' @return Named list (one entry per pathogen) with columns \code{profile},
#'   \code{probability}, \code{rr_profile_col}, \code{numerator}, \code{PAF_LOS},
#'   and \code{denominator} added to the input profile data frame.
#' @export

daly_calc_paf_los <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col, profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_kd  (sums to 1)
    rr <- df[[rr_profile_col]] # RR_LOS_kd
    numerator_vec <- p * (rr - 1.0) # R'_kd * (RR_kd - 1)
    denom <- 1.0 + sum(numerator_vec, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) {
      warning(sprintf(
        "'%s': PAF denominator = %.6g (must be > 0) -- all relative LOS may be < 1 or NA; skipping.",
        path, denom
      ))
      next
    }
    paf_vec <- numerator_vec / denom

    df$numerator <- round(numerator_vec, 6L)
    df$PAF_LOS <- round(paf_vec, 6L)
    df$denominator <- round(denom, 6L)

    paf_k <- sum(paf_vec)

    out[[path]] <- list(
      per_profile = df,
      PAF_k       = round(paf_k, 6L),
      denominator = round(denom, 6L)
    )
    message(sprintf(
      "'%s': PAF_k = %.4f | E[relative LOS] (denominator) = %.4f | %d profiles.",
      path, paf_k, denom, nrow(df)
    ))
  }

  return(out)
}


# -- Direct profile-specific YLD equations (associated / attributable) --------
#
# Shared machinery behind daly_calc_yld_associated() and
# daly_calc_yld_attributable(). Not exported; both public functions are thin
# wrappers so the two equations cannot drift out of sync with each other.
.daly_yld_direct <- function(
  mode,
  profiles_with_los,
  P_Lk_prime_tbl,
  yld_ref,
  DW_sepsis,
  pathogen_col,
  plk_col,
  incidence_col,
  facility_col,
  facility_name,
  facility_state_map,
  state_col,
  state_name,
  probability_col,
  dominant_class_col,
  losr_col,
  loss_col,
  out_col
) {
  # -- Input validation -------------------------------------------------------
  if (!is.list(profiles_with_los) || is.data.frame(profiles_with_los)) {
    stop(
      "profiles_with_los must be the named list returned by ",
      "daly_assign_rr_to_profiles() (called with los_r_col/los_s_col set)."
    )
  }
  if (!is.data.frame(P_Lk_prime_tbl)) {
    stop(
      "P_Lk_prime_tbl must be a data frame (the P_Lk_prime or ",
      "facility_level element from daly_calc_pathogen_fraction_nonfatal())."
    )
  }

  use_scalar_dw <- !is.null(DW_sepsis)
  if (use_scalar_dw) {
    if (!is.numeric(DW_sepsis) || length(DW_sepsis) != 1 || is.na(DW_sepsis)) {
      stop("DW_sepsis must be a single non-missing numeric value.")
    }
  } else {
    if (is.null(yld_ref) || !all(c("location_name", "DW_sepsis") %in% names(yld_ref))) {
      stop(
        "Provide either DW_sepsis as a numeric scalar, or yld_ref ",
        "with columns: 'location_name', 'DW_sepsis'."
      )
    }
    warning(
      "DW_L is being sourced from yld_ref (Proxy_YLD_per_case.csv-derived). ",
      "That table's 'DW_sepsis' column is an already duration-inclusive ",
      "YLD-per-case proxy (proxy_yld_per_case == yld_days_per_case / 365), ",
      "not a pure disability weight -- combining it with profile-specific ",
      "LOSR/LOSS here risks double-counting duration. Prefer supplying a ",
      "pure disability weight via the DW_sepsis scalar argument for this ",
      "direct-equation pathway.",
      call. = FALSE
    )
  }

  required_tbl_cols <- c(pathogen_col, plk_col, incidence_col)
  missing_tbl_cols <- setdiff(required_tbl_cols, names(P_Lk_prime_tbl))
  if (length(missing_tbl_cols) > 0L) {
    stop(sprintf(
      paste0(
        "Column(s) not found in P_Lk_prime_tbl: %s. incidence_col defaults ",
        "to 'N_NF_L', the non-fatal incidence column already returned by ",
        "daly_calc_pathogen_fraction_nonfatal()."
      ),
      paste(missing_tbl_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  # -- Optional single-facility restriction -----------------------------------
  if (!is.null(facility_name) && !is.null(facility_col) &&
    facility_col %in% names(P_Lk_prime_tbl)) {
    P_Lk_prime_tbl <- P_Lk_prime_tbl[P_Lk_prime_tbl[[facility_col]] == facility_name, , drop = FALSE]
    if (nrow(P_Lk_prime_tbl) == 0L) {
      stop(sprintf("No rows in P_Lk_prime_tbl for facility '%s'.", facility_name))
    }
  }

  n_row <- nrow(P_Lk_prime_tbl)

  # -- Resolve DW_L: one value per row of P_Lk_prime_tbl -----------------------
  has_facility_col <- !is.null(facility_col) && facility_col %in% names(P_Lk_prime_tbl)

  if (use_scalar_dw) {
    dw_vec <- rep(as.numeric(DW_sepsis), n_row)
  } else if (has_facility_col && !is.null(facility_state_map)) {
    if (!all(c(facility_col, state_col) %in% names(facility_state_map))) {
      stop(sprintf(
        "facility_state_map must have columns '%s' and '%s'.",
        facility_col, state_col
      ))
    }
    state_lookup <- stats::setNames(facility_state_map[[state_col]], facility_state_map[[facility_col]])
    dw_lookup <- stats::setNames(yld_ref$DW_sepsis, yld_ref$location_name)
    india_dw <- unname(dw_lookup[["India"]])
    fac_state <- unname(state_lookup[P_Lk_prime_tbl[[facility_col]]])
    dw_vec <- unname(dw_lookup[fac_state])
    n_missing_dw <- sum(is.na(dw_vec))
    if (n_missing_dw > 0L) {
      warning(sprintf(
        "%d row(s) had no state-level DW_sepsis match; using India-wide fallback.",
        n_missing_dw
      ))
      dw_vec[is.na(dw_vec)] <- india_dw
    }
  } else {
    loc <- if (is.null(state_name)) "India" else state_name
    dw_scalar <- yld_ref$DW_sepsis[yld_ref$location_name == loc]
    if (length(dw_scalar) == 0L) {
      stop(sprintf(
        "Location '%s' not found in yld_ref. Available: %s",
        loc, paste(utils::head(yld_ref$location_name, 10), collapse = ", ")
      ))
    }
    dw_vec <- rep(as.numeric(dw_scalar[1]), n_row)
  }

  # -- Per-row profile term: sum over RESISTANT profiles only ------------------
  # D_k = profiles with a dominant class (dominant_class != "all_susceptible").
  # The all-susceptible profile is excluded rather than included with an
  # implied LOSR==LOSS(->0) term, matching the pre-existing convention
  # elsewhere in this package that "associated"/"attributable" burden is a
  # partition over resistant profiles only (see daly_calc_paf_los(), whose
  # all-susceptible term already contributes exactly zero by construction).
  profile_term <- rep(NA_real_, n_row)
  n_unmatched_pathogen <- 0L

  for (i in seq_len(n_row)) {
    path <- as.character(P_Lk_prime_tbl[[pathogen_col]][i])
    prof_df <- profiles_with_los[[path]]

    if (is.null(prof_df)) {
      n_unmatched_pathogen <- n_unmatched_pathogen + 1L
      next
    }

    required_prof_cols <- c(probability_col, dominant_class_col, losr_col)
    if (mode == "attributable") required_prof_cols <- c(required_prof_cols, loss_col)
    missing_prof_cols <- setdiff(required_prof_cols, names(prof_df))
    if (length(missing_prof_cols) > 0L) {
      stop(sprintf(
        paste0(
          "Column(s) not found in profiles_with_los[['%s']]: %s. Call ",
          "daly_assign_rr_to_profiles() with los_r_col/los_s_col set to ",
          "produce them."
        ),
        path, paste(missing_prof_cols, collapse = ", ")
      ))
    }

    resist <- prof_df[prof_df[[dominant_class_col]] != "all_susceptible", , drop = FALSE]
    if (nrow(resist) == 0L) {
      profile_term[i] <- 0
      next
    }

    p <- resist[[probability_col]]
    losr <- resist[[losr_col]]
    n_na_losr <- sum(is.na(losr))
    if (n_na_losr > 0L) {
      warning(sprintf(
        paste0(
          "'%s': %d resistant profile(s) have no matched LOSR_years ",
          "(dominant class not in rr_table) -- excluded from the profile sum."
        ),
        path, n_na_losr
      ))
    }

    if (mode == "associated") {
      term_vec <- p * losr
    } else {
      loss <- resist[[loss_col]]
      n_na_loss <- sum(is.na(loss) & !is.na(losr))
      if (n_na_loss > 0L) {
        warning(sprintf(
          paste0(
            "'%s': %d resistant profile(s) have LOSR_years but no matched ",
            "LOSS_years -- excluded from the profile sum."
          ),
          path, n_na_loss
        ))
      }
      term_vec <- p * (losr - loss)
    }

    profile_term[i] <- sum(term_vec, na.rm = TRUE)
  }

  if (n_unmatched_pathogen > 0L) {
    warning(sprintf(
      "%d row(s) in P_Lk_prime_tbl had no matching entry in profiles_with_los; %s set to NA.",
      n_unmatched_pathogen, out_col
    ))
  }

  out <- P_Lk_prime_tbl
  out$I_L <- out[[incidence_col]]
  out$DW_L <- dw_vec
  out$profile_term <- round(profile_term, 6L)
  out[[out_col]] <- out$I_L * out[[plk_col]] * out$DW_L * out$profile_term

  message(sprintf(
    "%s computed: %d row(s), total %s = %.4f.",
    out_col, n_row, out_col, sum(out[[out_col]], na.rm = TRUE)
  ))

  out
}


#' Compute YLD Associated with Resistance (Direct, Profile-Specific)
#'
#' Primary, direct implementation of the profile-specific associated-YLD
#' equation:
#'
#'   YLD_assoc_k = I_L * P'_Lk * DW_L *
#'     sum_\{delta in D_k\} [R'_k_delta * LOSR_k,d*(delta)]
#'
#' where D_k is the set of RESISTANT resistance profiles for pathogen k (the
#' all-susceptible profile has no dominant class and is excluded -- see
#' "Details"). d*(delta) is the dominant drug class of profile delta (GBD max
#' rule, assigned by \code{daly_assign_rr_to_profiles()}).
#'
#' This does \strong{not} require computing a baseline/pooled YLD first:
#' unlike the pre-existing multiplicative pathway (a pooled baseline YLD
#' multiplied by an associated-burden fraction), this is the direct
#' equation, evaluated per pathogen (and per facility, if
#' \code{P_Lk_prime_tbl} is facility-level) from its own inputs.
#'
#' \strong{Component sourcing}:
#' \itemize{
#'   \item \strong{I_L} and \strong{P'_Lk}: read directly from
#'     \code{P_Lk_prime_tbl} -- the \code{P_Lk_prime} (pooled) or
#'     \code{facility_level} element already returned by
#'     \code{daly_calc_pathogen_fraction_nonfatal()}, which computes both the
#'     non-fatal incidence count and the pathogen fraction from the same
#'     observed cohort. Incidence is never back-calculated from deaths, CFR,
#'     or LOS.
#'   \item \strong{DW_L}: resolved exactly as in the package's pre-existing
#'     YLD weight mechanism -- either the \code{DW_sepsis} scalar, or a
#'     lookup in \code{yld_ref} (loaded from
#'     \code{inst/extdata/Proxy_YLD_per_case.csv}) by state/facility. See
#'     "Disability weight caveat" below.
#'   \item \strong{R'_k_delta}, \strong{d*(delta)}, and \strong{LOSR}: read
#'     from \code{profiles_with_los}, the list returned by
#'     \code{daly_assign_rr_to_profiles()} when called with
#'     \code{los_r_col}/\code{los_s_col} set.
#' }
#'
#' \strong{Why D_k excludes the all-susceptible profile}: this matches the
#' pre-existing convention elsewhere in this package that "associated"
#' burden is a partition over resistant profiles only (the removed
#' \code{daly_calc_fraction_associated_yld()}'s \code{Fraction_k} summed over
#' resistant profiles only; \code{daly_calc_paf_los()}'s all-susceptible term
#' already contributes exactly zero by construction, since RR=1 there).
#' \code{YLD_associated} therefore represents the burden occurring
#' specifically among pathogen-k patients with a resistant profile, not the
#' total burden across all pathogen-k patients.
#'
#' \strong{Disability weight caveat}: \code{Proxy_YLD_per_case.csv}'s
#' \code{DW_sepsis} column is empirically a duration-inclusive "YLD per case"
#' proxy (\code{proxy_yld_per_case == yld_days_per_case / 365}), not a pure
#' disability weight. Because this function already multiplies by
#' profile-specific \code{LOSR_years} separately, sourcing \code{DW_L} from
#' \code{yld_ref} risks double-counting duration; a warning is emitted when
#' \code{DW_sepsis} is not supplied directly. Prefer supplying a pure
#' disability weight via the \code{DW_sepsis} scalar for this pathway.
#'
#' @param profiles_with_los Named list (one entry per pathogen) from
#'   \code{daly_assign_rr_to_profiles(..., los_r_col = ..., los_s_col = ...)}.
#'   Each entry must have \code{probability_col}, \code{dominant_class_col},
#'   and \code{losr_col} columns.
#' @param P_Lk_prime_tbl Data frame: the \code{P_Lk_prime} (pooled) or
#'   \code{facility_level} element from
#'   \code{daly_calc_pathogen_fraction_nonfatal()}. Must contain
#'   \code{pathogen_col}, \code{plk_col}, and \code{incidence_col}.
#' @param yld_ref Data frame with columns \code{location_name} and
#'   \code{DW_sepsis}. Ignored when \code{DW_sepsis} is supplied.
#' @param DW_sepsis Numeric scalar or \code{NULL}. Disability weight used
#'   directly for every row when supplied (recommended -- see "Disability
#'   weight caveat").
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param plk_col Character. P'_Lk column in \code{P_Lk_prime_tbl}. Default
#'   \code{"P_Lk_prime"}.
#' @param incidence_col Character. I_L column in \code{P_Lk_prime_tbl}.
#'   Default \code{"N_NF_L"} (the column already produced by
#'   \code{daly_calc_pathogen_fraction_nonfatal()}).
#' @param facility_col Character or \code{NULL}. Facility identifier column
#'   in \code{P_Lk_prime_tbl}. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. Restrict to a single
#'   facility. Default \code{NULL}.
#' @param facility_state_map Data frame with \code{facility_col} and
#'   \code{state_col}, mapping each facility to a state. Used only when
#'   \code{DW_sepsis} is not supplied and \code{facility_col} is present.
#' @param state_col Character. Default \code{"state"}.
#' @param state_name Character or \code{NULL}. State for the pooled/no-facility
#'   DW lookup. \code{NULL} uses the "India" row. Ignored when
#'   \code{DW_sepsis} is supplied.
#' @param probability_col Character. R'_k_delta column in
#'   \code{profiles_with_los}. Default \code{"probability"}.
#' @param dominant_class_col Character. d*(delta) column. Default
#'   \code{"dominant_class"}.
#' @param losr_col Character. LOSR_k,d*(delta) column (years). Default
#'   \code{"LOSR_years"}.
#'
#' @return \code{P_Lk_prime_tbl} augmented with \code{I_L}, \code{DW_L},
#'   \code{profile_term} (= sum_delta R'_k_delta * LOSR_k_delta), and
#'   \code{YLD_associated}.
#' @export
daly_calc_yld_associated <- function(
  profiles_with_los,
  P_Lk_prime_tbl,
  yld_ref = NULL,
  DW_sepsis = NULL,
  pathogen_col = "pathogen",
  plk_col = "P_Lk_prime",
  incidence_col = "N_NF_L",
  facility_col = NULL,
  facility_name = NULL,
  facility_state_map = NULL,
  state_col = "state",
  state_name = NULL,
  probability_col = "probability",
  dominant_class_col = "dominant_class",
  losr_col = "LOSR_years"
) {
  .daly_yld_direct(
    mode = "associated",
    profiles_with_los = profiles_with_los,
    P_Lk_prime_tbl = P_Lk_prime_tbl,
    yld_ref = yld_ref,
    DW_sepsis = DW_sepsis,
    pathogen_col = pathogen_col,
    plk_col = plk_col,
    incidence_col = incidence_col,
    facility_col = facility_col,
    facility_name = facility_name,
    facility_state_map = facility_state_map,
    state_col = state_col,
    state_name = state_name,
    probability_col = probability_col,
    dominant_class_col = dominant_class_col,
    losr_col = losr_col,
    loss_col = NULL,
    out_col = "YLD_associated"
  )
}


#' Compute YLD Attributable to Resistance (Direct, Profile-Specific)
#'
#' Primary, direct implementation of the profile-specific attributable-YLD
#' equation:
#'
#'   YLD_attrib_k = I_L * P'_Lk * DW_L *
#'     sum_\{delta in D_k\} [R'_k_delta * (LOSR_k,d*(delta) - LOSS_k,d*(delta))]
#'
#' The excess \code{(LOSR - LOSS)} term makes this a genuine counterfactual:
#' how much MORE disability burden exists because these infections were
#' resistant, versus if they had been susceptible to the same (dominant)
#' class. D_k excludes the all-susceptible profile from the sum entirely
#' (rather than including a zero \code{LOSR - LOSS} term for it) -- see
#' \code{daly_calc_yld_associated()} for why.
#'
#' Obtains I_L, P'_Lk, DW_L, R'_k_delta, d*(delta), LOSR, and LOSS the same
#' way as \code{daly_calc_yld_associated()} -- see that function's
#' documentation for full component-by-component sourcing and the
#' disability-weight caveat. Does \strong{not} require a baseline YLD or a
#' PAF_LOS computed first; unlike the pre-existing multiplicative pathway
#' (\code{daly_calc_paf_los()} multiplied against a pooled baseline YLD),
#' this is the direct equation. \code{daly_calc_paf_los()} remains available
#' standalone for diagnostics / comparison against this direct calculation.
#'
#' @inheritParams daly_calc_yld_associated
#' @param loss_col Character. LOSS_k,d*(delta) column (years) in
#'   \code{profiles_with_los}. Default \code{"LOSS_years"}.
#'
#' @return \code{P_Lk_prime_tbl} augmented with \code{I_L}, \code{DW_L},
#'   \code{profile_term} (= sum_delta R'_k_delta * (LOSR_k_delta -
#'   LOSS_k_delta)), and \code{YLD_attributable}.
#' @export
daly_calc_yld_attributable <- function(
  profiles_with_los,
  P_Lk_prime_tbl,
  yld_ref = NULL,
  DW_sepsis = NULL,
  pathogen_col = "pathogen",
  plk_col = "P_Lk_prime",
  incidence_col = "N_NF_L",
  facility_col = NULL,
  facility_name = NULL,
  facility_state_map = NULL,
  state_col = "state",
  state_name = NULL,
  probability_col = "probability",
  dominant_class_col = "dominant_class",
  losr_col = "LOSR_years",
  loss_col = "LOSS_years"
) {
  .daly_yld_direct(
    mode = "attributable",
    profiles_with_los = profiles_with_los,
    P_Lk_prime_tbl = P_Lk_prime_tbl,
    yld_ref = yld_ref,
    DW_sepsis = DW_sepsis,
    pathogen_col = pathogen_col,
    plk_col = plk_col,
    incidence_col = incidence_col,
    facility_col = facility_col,
    facility_name = facility_name,
    facility_state_map = facility_state_map,
    state_col = state_col,
    state_name = state_name,
    probability_col = probability_col,
    dominant_class_col = dominant_class_col,
    losr_col = losr_col,
    loss_col = loss_col,
    out_col = "YLD_attributable"
  )
}
