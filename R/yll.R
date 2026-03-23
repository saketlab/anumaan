# yll.R
#
# YLL calculation pipeline functions.
#
# Functions in this file:
#   calculate_P_Lk          -- Fatal pathogen distribution P_{LK} (Step 4)
#   compute_base_yll_from_dl -- Base YLL block sum_x(D_x^L * e_x*) (Steps 1-3)
#
# Pipeline (Bhaswati Ganguli, DALY Methodology for AMR, 2026):
#
#   Step 1.  D_L = deaths by syndrome L (supplied by the caller as a scalar).
#            For 1000-incidence normalisation: D_L = 1000 x death_rate_L
#            where death_rate_L = (#deaths | L) / (#incidence_L)
#            and   #incidence_L  = #prevalence_L / (obs_period x avg_LOS).
#
#   Step 2.  Disaggregate D_L by age (x sex) using observed proportions from
#            the patient death cohort, yielding {D_x^L} for each age bin x
#            (same bins as the India life table).
#
#   Step 3.  Compute the base YLL block (Eq. 10):
#              sum_x  D_x^L  x  e_x*
#            where e_x* = expectation of life at age bin x (India life table).
#
# Aggregation levels (within the same function):
#   - Overall (all deaths for the syndrome)
#   - Per pathogen K  (using each pathogen's own age distribution)
#   - Per age x sex
#   - Per pathogen x age x sex
#   - Per facility   (using each facility's own age distribution)
#   - Per syndrome x pathogen
#
# This base quantity feeds into:
#   YLL_Kd  (Eq. 11)  = [sum_L base_YLL_L x P_KL] x R_Kd
#   YLL_K   (Eq. 12)  = sum_delta YLL_Kdelta
#   YLL_Kdelta  (Eq. 14)  = [sum_L base_YLL_L x P_KL] x Mort_PAF_Kdelta
#
# References:
#   Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026.
#   Antimicrobial Resistance Collaborators. Lancet. 2022.


# -- compute_base_yll_from_dl --------------------------------------------------

#' Compute Base YLL Block from a Scalar D_L (Steps 1-3)
#'
#' Implements the first part of the YLL pipeline. Takes \eqn{D_L} (deaths by
#' infectious syndrome, a single scalar), disaggregates it into age x sex
#' strata using observed proportions from the patient death cohort, and
#' multiplies by the India life expectancy for each stratum to yield the base
#' YLL block \eqn{\sum_x D_L^x \, e_*^x} (Equation 10).
#'
#' \strong{Formula:}
#' \deqn{
#'   \text{base\_YLL}_L = \sum_x D_L^x \, e_*^x
#' }
#' where
#' \deqn{D_L^x = D_L \times \hat{p}_x}
#' and \eqn{\hat{p}_x} is the observed proportion of deaths in age bin
#' \eqn{x} (within sex \eqn{s} when \code{use_sex = TRUE}), estimated from
#' \code{patient_data} (filtered to deaths, and optionally to
#' \code{syndrome_name}).
#'
#' \strong{D_L note:} For 1000-incidence normalisation:
#' \deqn{D_L = 1000 \times \text{death\_rate}_L,\quad
#'   \text{death\_rate}_L = \frac{\#\text{deaths}\mid L}{\#\text{incidence}_L}}
#' The caller supplies this scalar directly via \code{dl}.
#'
#' \strong{Per-pathogen / per-facility YLL:} age (x sex) proportions are
#' re-computed within each subgroup so that each subgroup's YLL reflects its
#' own age structure, scaled by the shared \code{dl}.
#'
#' @param dl Numeric scalar. \eqn{D_L}, the number of deaths for syndrome L
#'   (e.g. \code{45.2} deaths per 1000 incidences).
#' @param patient_data Data frame. Patient-level records (one row per
#'   patient x antibiotic test or per patient x pathogen). Used to derive
#'   observed age (x sex) proportions from the death cohort.
#' @param patient_col Character. Unique patient identifier column in
#'   \code{patient_data}.
#' @param outcome_col Character. Final-outcome column in \code{patient_data}.
#' @param death_value Character (scalar or vector). Value(s) in
#'   \code{outcome_col} indicating a fatal outcome. Default \code{"Death"}.
#' @param age_bin_col Character. Column in \code{patient_data} containing
#'   GBD-standard age bin labels (e.g. \code{"0-1"}, \code{"1-5"}, ...,
#'   \code{"85+"}).  Use \code{age_bin_map} to recode non-standard labels
#'   (default remaps \code{"<1"} -> \code{"0-1"}).
#' @param sex_col Character. Sex column in \code{patient_data}.
#' @param syndrome_col Character or \code{NULL}. Syndrome column in
#'   \code{patient_data} (e.g. \code{"infectious_syndrome"}).  Required when
#'   \code{syndrome_name} is not \code{NULL}.
#' @param syndrome_name Character or \code{NULL}. If supplied, \code{patient_data}
#'   is filtered to rows where \code{syndrome_col == syndrome_name} before
#'   computing proportions and outputs (e.g. \code{"Bloodstream infection"}).
#'   \code{NULL} uses all syndromes.
#' @param pathogen_col Character or \code{NULL}. Pathogen (organism) column.
#'   When supplied, per-pathogen age proportions are computed and the
#'   \code{per_pathogen} and \code{by_pathogen_age_sex} outputs are populated.
#' @param facility_col Character or \code{NULL}. Facility identifier column.
#'   When supplied, per-facility proportions are computed and \code{by_facility}
#'   is populated.
#' @param facility_name Character or \code{NULL}. If provided, filters
#'   \code{patient_data} to the specified facility before computation.
#' @param use_sex Logical. Whether to disaggregate by sex as well as age
#'   (uses sex-specific life expectancy). Default \code{TRUE}.
#' @param le_path Character. Path to the India life expectancy xlsx file.
#'   Defaults to the bundled \code{inst/extdata} copy.
#' @param male_value Character. Value in \code{sex_col} for males.
#'   Default \code{"Male"}.
#' @param female_value Character. Value in \code{sex_col} for females.
#'   Default \code{"Female"}.  All other values map to \code{"Combined"}.
#' @param age_bin_map Named character vector. Remaps non-standard age bin
#'   labels before joining to the life table.  Default \code{c("<1" = "0-1")}.
#' @param stratify_by Character vector or \code{NULL}. Additional columns
#'   from \code{patient_data} to include in the \code{stratified} output.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{total}}{Scalar: total base YLL = \eqn{\sum_x D_L^x e_*^x}.}
#'   \item{\code{by_age_sex}}{Data frame: base YLL by age bin x sex.}
#'   \item{\code{per_pathogen}}{Data frame: base YLL per pathogen K, computed
#'     using each pathogen's own death age distribution (only when
#'     \code{pathogen_col} is supplied).}
#'   \item{\code{by_pathogen_age_sex}}{Data frame: base YLL by pathogen x
#'     age bin x sex (only when \code{pathogen_col} is supplied).}
#'   \item{\code{by_facility}}{Data frame: base YLL by facility (only when
#'     \code{facility_col} is supplied).}
#'   \item{\code{by_syndrome_pathogen}}{Data frame: base YLL by syndrome x
#'     pathogen (only when both \code{syndrome_col} and \code{pathogen_col}
#'     are supplied).}
#'   \item{\code{stratified}}{Data frame: base YLL aggregated by
#'     \code{stratify_by} (only when \code{stratify_by} is supplied).}
#'   \item{\code{disaggregated_dl}}{Data frame: the full expanded table with
#'     columns \code{D_x_L} (\eqn{D_L^x}), \code{proportion} (\eqn{\hat{p}_x}),
#'     \code{life_expectancy} (\eqn{e_*^x}), and \code{yll_contribution}
#'     (\eqn{D_L^x \times e_*^x}) for every age x sex stratum. This is the
#'     row-level audit trail.}
#' }
#'
#' @export
#'
#' @references
#' Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026.
#'
#' Antimicrobial Resistance Collaborators. Global burden of bacterial
#' antimicrobial resistance in 2019. Lancet. 2022.
#'
#' @examples
#' \dontrun{
#' result <- compute_base_yll_from_dl(
#'   dl = 45.2,
#'   patient_data = bsi_data,
#'   outcome_col = "final_outcome",
#'   death_value = "Death",
#'   age_bin_col = "Age_bin",
#'   sex_col = "gender",
#'   syndrome_col = "infectious_syndrome",
#'   syndrome_name = "Bloodstream infection",
#'   pathogen_col = "organism_name",
#'   facility_col = "center_name",
#'   le_path = here::here(
#'     "anumaan", "inst", "extdata",
#'     "life_expectancy_all.xlsx"
#'   )
#' )
#'
#' result$total
#' result$by_age_sex
#' result$per_pathogen
#' result$by_syndrome_pathogen
#' result$disaggregated_dl
#' }
compute_base_yll_from_dl <- function(
  dl,
  patient_data,
  patient_col,
  outcome_col,
  death_value = "Death",
  age_bin_col,
  sex_col,
  syndrome_col = NULL,
  syndrome_name = NULL,
  pathogen_col = NULL,
  facility_col = NULL,
  facility_name = NULL,
  use_sex = TRUE,
  le_path = system.file("extdata", "life_table_india.csv",
    package = "anumaan"
  ),
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c("<1" = "0-1"),
  stratify_by = NULL
) {
  # -- 1. Input validation ---------------------------------------------------

  if (!is.numeric(dl) || length(dl) != 1L || is.na(dl) || dl < 0) {
    stop("'dl' must be a single non-negative numeric value (D_L).")
  }

  required_pt <- c(patient_col, outcome_col, age_bin_col, sex_col)
  if (!is.null(syndrome_col)) required_pt <- c(required_pt, syndrome_col)
  if (!is.null(pathogen_col)) required_pt <- c(required_pt, pathogen_col)
  if (!is.null(facility_col)) required_pt <- c(required_pt, facility_col)
  if (!is.null(stratify_by)) required_pt <- c(required_pt, stratify_by)

  missing_pt <- setdiff(required_pt, names(patient_data))
  if (length(missing_pt) > 0L) {
    stop(sprintf(
      "Missing column(s) in patient_data: %s",
      paste(missing_pt, collapse = ", ")
    ))
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("'syndrome_col' must be supplied when 'syndrome_name' is specified.")
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("'facility_col' must be supplied when 'facility_name' is specified.")
  }

  # -- 2. Load India life expectancy lookup ----------------------------------
  # CSV columns: age_bin_gbd, life_expectancy, sex ("Both"/"Male"/"Female")

  if (!file.exists(le_path)) {
    stop(sprintf("Life table CSV not found: %s", le_path))
  }

  le_lookup <- utils::read.csv(le_path, stringsAsFactors = FALSE)

  # -- 3. Filter patient_data to death cohort --------------------------------

  deaths_df <- patient_data %>%
    dplyr::filter(.data[[outcome_col]] %in% death_value)

  if (!is.null(syndrome_name)) {
    deaths_df <- deaths_df %>%
      dplyr::filter(.data[[syndrome_col]] == syndrome_name)
  }
  if (!is.null(facility_name)) {
    deaths_df <- deaths_df %>%
      dplyr::filter(.data[[facility_col]] == facility_name)
  }

  if (nrow(deaths_df) == 0L) {
    stop(paste0(
      "No fatal records found after filtering. ",
      "Check outcome_col, death_value, syndrome_name, and facility_name."
    ))
  }

  message(sprintf(
    "Death cohort: %d records%s%s | D_L = %.4f",
    nrow(deaths_df),
    if (!is.null(syndrome_name)) sprintf(" [%s]", syndrome_name) else "",
    if (!is.null(facility_name)) sprintf(" [%s]", facility_name) else "",
    dl
  ))

  # -- 4. Recode age bins; normalise sex -------------------------------------

  if (length(age_bin_map) > 0L) {
    deaths_df[[age_bin_col]] <- dplyr::recode(
      as.character(deaths_df[[age_bin_col]]), !!!age_bin_map
    )
  }

  # Map patient sex to CSV sex labels: "Male", "Female", "Both"
  sex_strat_col <- ".sex_norm"
  deaths_df <- deaths_df %>%
    dplyr::mutate(
      .sex_norm = if (use_sex) {
        dplyr::case_when(
          .data[[sex_col]] == male_value ~ "Male",
          .data[[sex_col]] == female_value ~ "Female",
          TRUE ~ "Both"
        )
      } else {
        "Both"
      }
    )

  # -- 4b. Deduplicate to unique patient x pathogen --------------------------
  # Raw data is long (one row per antibiotic test). Collapse to one row per
  # patient x pathogen so every patient counts exactly once per organism
  # when computing age/sex proportions and n_deaths.
  dedup_cols <- unique(c(
    patient_col, age_bin_col, sex_strat_col,
    pathogen_col, facility_col, syndrome_col
  ))
  dedup_cols <- dedup_cols[!is.null(dedup_cols) & dedup_cols %in% names(deaths_df)]

  deaths_df <- dplyr::distinct(
    deaths_df,
    dplyr::across(dplyr::all_of(dedup_cols))
  )

  message(sprintf(
    "After dedup: %d unique patient-pathogen records | %d unique patients",
    nrow(deaths_df),
    dplyr::n_distinct(deaths_df[[patient_col]])
  ))

  # -- Helper: compute proportions + YLL for a given subgroup of deaths ------
  #
  # For a data frame of deaths (already filtered to a subgroup), compute:
  #   1. age (x sex) proportions
  #   2. D_x^L = dl * proportion
  #   3. join life expectancy e_x*
  #   4. yll_contribution = D_x^L * e_x*
  # Returns a tidy data frame with one row per age (x sex) stratum.

  .compute_yll_strata <- function(df_deaths, label = "global") {
    counts <- df_deaths %>%
      dplyr::group_by(
        dplyr::across(dplyr::all_of(c(age_bin_col, sex_strat_col)))
      ) %>%
      dplyr::summarise(n_deaths = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        n_deaths_total = sum(n_deaths),
        proportion     = n_deaths / n_deaths_total,
        D_x_L          = dl * proportion
      )

    # Join life expectancy -- CSV uses age_bin_gbd and sex columns
    join_le <- stats::setNames(
      c("age_bin_gbd", "sex"),
      c(age_bin_col, sex_strat_col)
    )
    counts <- dplyr::left_join(counts, le_lookup, by = join_le)

    n_miss <- sum(is.na(counts$life_expectancy))
    if (n_miss > 0L) {
      warning(sprintf(
        "[%s] %d age/sex bin(s) had no life expectancy match; YLL set to NA.",
        label, n_miss
      ))
    }

    counts %>%
      dplyr::mutate(yll_contribution = D_x_L * life_expectancy)
  }

  # -- 5. Overall disaggregation (all deaths for the syndrome) ---------------

  overall_strata <- .compute_yll_strata(deaths_df, label = "overall")

  YLL_total <- sum(overall_strata$yll_contribution, na.rm = TRUE)

  by_age_sex <- overall_strata %>%
    dplyr::rename(sex = dplyr::all_of(sex_strat_col)) %>%
    dplyr::select(
      dplyr::all_of(c(age_bin_col, "sex")),
      proportion, D_x_L, life_expectancy, yll_contribution
    )

  # -- 6. Per-pathogen disaggregation ----------------------------------------
  #
  # For each pathogen K, re-compute age proportions using only that
  # pathogen's deaths. This reflects K's own age burden.

  per_pathogen <- NULL
  by_pathogen_age_sex <- NULL

  if (!is.null(pathogen_col)) {
    pathogens <- sort(unique(deaths_df[[pathogen_col]]))

    path_strata_list <- lapply(pathogens, function(k) {
      df_k <- deaths_df %>% dplyr::filter(.data[[pathogen_col]] == k)
      strata_k <- .compute_yll_strata(df_k, label = k)
      strata_k[[pathogen_col]] <- k
      strata_k
    })
    path_strata <- dplyr::bind_rows(path_strata_list)

    by_pathogen_age_sex <- path_strata %>%
      dplyr::rename(sex = dplyr::all_of(sex_strat_col)) %>%
      dplyr::select(
        dplyr::all_of(c(pathogen_col, age_bin_col, "sex")),
        proportion, D_x_L, life_expectancy, yll_contribution
      )

    per_pathogen <- by_pathogen_age_sex %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(pathogen_col))) %>%
      dplyr::summarise(
        YLL_base = sum(yll_contribution, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        deaths_df %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(pathogen_col))) %>%
          dplyr::summarise(n_deaths = dplyr::n(), .groups = "drop"),
        by = pathogen_col
      )
  }

  # -- 7. Per-facility disaggregation ----------------------------------------

  by_facility <- NULL

  if (!is.null(facility_col)) {
    facilities <- sort(unique(deaths_df[[facility_col]]))

    fac_strata_list <- lapply(facilities, function(f) {
      df_f <- deaths_df %>% dplyr::filter(.data[[facility_col]] == f)
      strata_f <- .compute_yll_strata(df_f, label = f)
      strata_f[[facility_col]] <- f
      strata_f
    })
    fac_strata <- dplyr::bind_rows(fac_strata_list)

    by_facility <- fac_strata %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(facility_col))) %>%
      dplyr::summarise(
        YLL_base = sum(yll_contribution, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        deaths_df %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(facility_col))) %>%
          dplyr::summarise(n_deaths = dplyr::n(), .groups = "drop"),
        by = facility_col
      )
  }

  # -- 8. By syndrome x pathogen ---------------------------------------------

  by_syndrome_pathogen <- NULL

  if (!is.null(syndrome_col) && !is.null(pathogen_col)) {
    syndromes <- sort(unique(deaths_df[[syndrome_col]]))

    syn_path_list <- lapply(syndromes, function(s) {
      df_s <- deaths_df %>% dplyr::filter(.data[[syndrome_col]] == s)
      pats <- sort(unique(df_s[[pathogen_col]]))

      lapply(pats, function(k) {
        df_sk <- df_s %>% dplyr::filter(.data[[pathogen_col]] == k)
        strata_sk <- .compute_yll_strata(df_sk,
          label = paste(s, k, sep = "|")
        )
        strata_sk[[syndrome_col]] <- s
        strata_sk[[pathogen_col]] <- k
        strata_sk
      })
    })
    syn_path_strata <- dplyr::bind_rows(unlist(syn_path_list, recursive = FALSE))

    by_syndrome_pathogen <- syn_path_strata %>%
      dplyr::group_by(
        dplyr::across(dplyr::all_of(c(syndrome_col, pathogen_col)))
      ) %>%
      dplyr::summarise(
        YLL_base = sum(yll_contribution, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # -- 9. User-defined stratification ---------------------------------------

  stratified <- NULL

  if (!is.null(stratify_by)) {
    strat_cols <- unique(c(stratify_by, pathogen_col))
    strat_cols <- strat_cols[!is.null(strat_cols)]

    strat_groups <- deaths_df %>%
      dplyr::select(dplyr::all_of(strat_cols)) %>%
      dplyr::distinct()

    strat_list <- lapply(seq_len(nrow(strat_groups)), function(i) {
      mask <- rep(TRUE, nrow(deaths_df))
      for (col in strat_cols) {
        mask <- mask & (deaths_df[[col]] == strat_groups[[col]][i])
      }
      df_sub <- deaths_df[mask, , drop = FALSE]
      if (nrow(df_sub) == 0L) {
        return(NULL)
      }
      strata_sub <- .compute_yll_strata(df_sub,
        label = paste(strat_groups[i, ], collapse = "|")
      )
      for (col in strat_cols) strata_sub[[col]] <- strat_groups[[col]][i]
      strata_sub
    })
    strat_all <- dplyr::bind_rows(strat_list)

    stratified <- strat_all %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(strat_cols))) %>%
      dplyr::summarise(
        YLL_base = sum(yll_contribution, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # -- 10. Summary message ---------------------------------------------------

  message(sprintf(
    "Base YLL (sum_x D_x^L * e_x*): %.4f years | D_L = %.4f | %d age-sex bins",
    YLL_total, dl, nrow(by_age_sex)
  ))

  # -- Return ----------------------------------------------------------------

  out <- list(
    total = YLL_total,
    by_age_sex = by_age_sex,
    disaggregated_dl = overall_strata %>%
      dplyr::rename(sex = dplyr::all_of(sex_strat_col))
  )
  if (!is.null(pathogen_col)) {
    out$per_pathogen <- per_pathogen
    out$by_pathogen_age_sex <- by_pathogen_age_sex
  }
  if (!is.null(facility_col)) out$by_facility <- by_facility
  if (!is.null(syndrome_col) && !is.null(pathogen_col)) out$by_syndrome_pathogen <- by_syndrome_pathogen
  if (!is.null(stratify_by)) out$stratified <- stratified

  out
}


# -- calculate_P_Lk ------------------------------------------------------------

#' Calculate Fatal Pathogen Distribution (P_\{LK\})
#'
#' Computes the fatal pathogen distribution for a given infectious syndrome
#' (L) using facility-level microbiology data. This quantity represents the
#' fractional contribution of each pathogen (K) to \strong{fatal} infection
#' cases, and is used in Step 4 of YLL estimation (Eq. 11).
#'
#' \strong{Unit of analysis:} the patient. Each fatal patient contributes a
#' total weight of 1, split equally across their valid pathogens:
#' \deqn{w_{r,k} = \frac{1}{m_r}}
#' where \eqn{m_r} is the number of distinct pathogens for patient \eqn{r}.
#' Monomicrobial patients (\eqn{m_r = 1}) receive weight 1; polymicrobial
#' patients have their death split equally (Approach 1).
#'
#' \strong{Pooled formula across facilities:}
#' \deqn{P_{LK}^{\text{pooled}} =
#'   \frac{\sum_f N^{F}_{f,L,K}}{\sum_f N^{F}_{f,L}}}
#'
#' @param data             Data frame of facility-level microbiology records.
#' @param syndrome_col     Character. Column containing infectious syndrome
#'   labels (L).
#' @param syndrome_name    Character. Syndrome to analyse (e.g.
#'   \code{"Bloodstream infection"}).
#' @param polymicrobial_col Character. Column flagging polymicrobial patients
#'   (1 = polymicrobial, 0 = monomicrobial).
#' @param patient_col      Character. Unique patient identifier column.
#' @param pathogen_col     Character. Pathogen (organism) column (K).
#' @param outcome_col      Character. Final patient outcome column.
#' @param death_value      Character. Value indicating a fatal outcome.
#'   Default \code{"Death"}.
#' @param specimen_col     Character or \code{NULL}. Specimen type column.
#'   When supplied, \code{specimen_name} must also be provided.
#' @param specimen_name    Character or \code{NULL}. Specimen value to
#'   restrict to (e.g. \code{"Blood"}).
#' @param glass_ref        Character vector of valid pathogen names, or a data
#'   frame with columns \code{specimen} and \code{pathogen}. When supplied,
#'   polymicrobial patients are restricted to pathogens on the GLASS list for
#'   the given specimen. Monomicrobial patients are never filtered.
#'   \code{NULL} skips GLASS filtering.
#' @param facility_col     Character or \code{NULL}. Facility identifier
#'   column. When supplied without \code{facility_name}, returns both
#'   facility-level and pooled P_LK.
#' @param facility_name    Character or \code{NULL}. Restricts computation to
#'   a single named facility.
#' @param pathogen_name    Character vector or \code{NULL}. Restricts output
#'   to specific pathogen(s).
#'
#' @return A named list:
#' \describe{
#'   \item{\code{P_Lk}}{Data frame of pooled P_LK across facilities. Columns:
#'     \code{pathogen_col}, \code{N_F_LK} (weighted fatal patient count for K),
#'     \code{N_F_L} (total fatal patients for L), \code{P_Lk} (fraction).}
#'   \item{\code{facility_level}}{Per-facility P_LK data frame (only when
#'     \code{facility_col} is supplied and \code{facility_name} is
#'     \code{NULL}).}
#' }
#'
#' @export
#'
#' @references
#' Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026.
#'
#' @examples
#' \dontrun{
#' pkL <- calculate_P_Lk(
#'   data              = bsi_data,
#'   syndrome_col      = "infectious_syndrome",
#'   syndrome_name     = "Bloodstream infection",
#'   polymicrobial_col = "is_polymicrobial",
#'   patient_col       = "PatientInformation_id",
#'   pathogen_col      = "organism_name",
#'   outcome_col       = "final_outcome",
#'   death_value       = "Death",
#'   facility_col      = "center_name"
#' )
#' pkL$P_Lk
#' pkL$facility_level
#' }
calculate_P_Lk <- function(data,
                           syndrome_col,
                           syndrome_name,
                           polymicrobial_col,
                           patient_col,
                           pathogen_col,
                           outcome_col,
                           death_value = "Death",
                           specimen_col = NULL,
                           specimen_name = NULL,
                           glass_ref = NULL,
                           facility_col = NULL,
                           facility_name = NULL,
                           pathogen_name = NULL) {
  # -- Input validation -------------------------------------------------------
  required_cols <- c(
    syndrome_col, polymicrobial_col,
    patient_col, pathogen_col, outcome_col
  )
  if (!is.null(specimen_col)) required_cols <- c(required_cols, specimen_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  if (!is.null(specimen_col) && is.null(specimen_name)) {
    stop("specimen_name must be provided when specimen_col is specified.")
  }

  # -- Step 1: Filter to syndrome + (optional specimen) + fatal outcome -------
  df <- data %>%
    dplyr::filter(
      .data[[syndrome_col]] == syndrome_name,
      .data[[outcome_col]] == death_value
    )

  if (!is.null(specimen_col)) {
    df <- df %>% dplyr::filter(.data[[specimen_col]] == specimen_name)
  }

  if (nrow(df) == 0) {
    stop(sprintf(
      "No fatal records remain after syndrome%s filtering (death_value = '%s').",
      if (!is.null(specimen_col)) "/specimen" else "",
      death_value
    ))
  }

  # -- Step 2: Optional single-facility restriction ---------------------------
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No fatal records found for facility '%s'.", facility_name))
    }
  }

  # -- Step 3: Optional pathogen filter --------------------------------------
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No fatal records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # -- Step 4: GLASS filter -- polymicrobial patients only --------------------
  # Monomicrobial patients are never filtered.
  if (!is.null(glass_ref)) {
    if (is.data.frame(glass_ref)) {
      valid_pathogens <- if (!is.null(specimen_name)) {
        glass_ref %>%
          dplyr::filter(.data[["specimen"]] == specimen_name) %>%
          dplyr::pull(.data[["pathogen"]])
      } else {
        unique(glass_ref[["pathogen"]])
      }
    } else {
      valid_pathogens <- glass_ref
    }

    df_mono <- df %>% dplyr::filter(.data[[polymicrobial_col]] == 0)
    df_poly <- df %>%
      dplyr::filter(
        .data[[polymicrobial_col]] == 1,
        .data[[pathogen_col]] %in% valid_pathogens
      )

    n_removed <- sum(df[[polymicrobial_col]] == 1) - nrow(df_poly)
    if (n_removed > 0) {
      message(sprintf(
        "GLASS filter: removed %d row(s) from polymicrobial patients%s.",
        n_removed,
        if (!is.null(specimen_name)) sprintf(" (not on GLASS list for '%s')", specimen_name) else ""
      ))
    }

    df <- dplyr::bind_rows(df_mono, df_poly)
    if (nrow(df) == 0) {
      stop("No fatal records remain after GLASS reference filtering.")
    }
  }

  # -- Step 5: Deduplicate to unique patient x pathogen, then weight --------
  # Raw data has one row per antibiotic tested. Collapse first so each
  # patient-pathogen pair contributes exactly once before weighting.
  # m_r = distinct valid pathogens for patient r.
  # Monomicrobial -> weight = 1; polymicrobial with m pathogens -> weight = 1/m.
  group_cols <- c(facility_col, patient_col) # c() drops NULL

  df <- df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_cols, pathogen_col)))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      m_r    = dplyr::n_distinct(.data[[pathogen_col]]),
      weight = 1 / m_r
    ) %>%
    dplyr::ungroup()

  # -- Step 6: N^F_LK = weighted patient count per (facility x) pathogen -----
  agg_cols <- c(facility_col, pathogen_col) # c() drops NULL

  N_LK <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) %>%
    dplyr::summarise(N_F_LK = sum(weight, na.rm = TRUE), .groups = "drop")

  # -- Step 7: N^F_L = unique fatal patients per (facility) ------------------
  fac_grp <- if (!is.null(facility_col)) facility_col else character(0)

  N_L <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(fac_grp))) %>%
    dplyr::summarise(
      N_F_L = dplyr::n_distinct(.data[[patient_col]]),
      .groups = "drop"
    )

  # -- Step 8: Compute P_LK --------------------------------------------------
  if (!is.null(facility_col)) {
    facility_level <- dplyr::left_join(N_LK, N_L, by = facility_col) %>%
      dplyr::mutate(P_Lk = N_F_LK / N_F_L)

    # Pooled: sum numerators / sum denominators across facilities
    pooled <- N_LK %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(N_F_LK = sum(N_F_LK), .groups = "drop") %>%
      dplyr::mutate(
        N_F_L = sum(N_L$N_F_L),
        P_Lk  = N_F_LK / N_F_L
      )

    message(sprintf(
      "P_LK computed: %d pathogen(s) across %d facility/facilities | sum = %.4f",
      dplyr::n_distinct(df[[pathogen_col]]),
      dplyr::n_distinct(df[[facility_col]]),
      sum(pooled$P_Lk)
    ))
    return(list(P_Lk = pooled, facility_level = facility_level))
  } else {
    N_total <- N_L$N_F_L
    result <- N_LK %>%
      dplyr::mutate(
        N_F_L = N_total,
        P_Lk  = N_F_LK / N_F_L
      )

    message(sprintf(
      "P_LK computed: %d pathogen(s) | N^F_L = %d fatal patient(s) | sum = %.4f",
      nrow(result), N_total, sum(result$P_Lk)
    ))
    return(list(P_Lk = result))
  }
}


# -- .build_class_resistance_wide ---------------------------------------------

#' Collapse drug-level data to antibiotic-class binary wide matrix
#'
#' Standardises antibiotic values (Intermediate -> R/S per GBD rules),
#' collapses to class level (class = 1 if ANY drug in class = R), and pivots
#' wide: one row per patient, one column per antibiotic class (0/1).
#' Column names are sanitised with make.names(); the original-to-safe name
#' mapping is stored as the "class_name_map" attribute.
#'
#' @keywords internal
.build_class_resistance_wide <- function(
  data,
  patient_id_col = "PatientInformation_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  untested_fill = 0L
) {
  abx_std <- data %>%
    dplyr::filter(
      !is.na(.data[[antibiotic_value_col]]),
      !is.na(.data[[antibiotic_class_col]])
    ) %>%
    dplyr::mutate(
      .abx_name = stringr::str_to_lower(stringr::str_trim(
        as.character(.data[[antibiotic_name_col]])
      )),
      .abx_val = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[antibiotic_value_col]])
      )),
      .abx_val = dplyr::case_when(
        .abx_name == "colistin" & .abx_val == "I" ~ "S",
        .abx_val == "I" ~ "R",
        .abx_val %in% c("R", "RESISTANT") ~ "R",
        .abx_val %in% c("S", "SUSCEPTIBLE", "SENSITIVE") ~ "S",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.abx_val))

  class_long <- abx_std %>%
    dplyr::group_by(
      .data[[patient_id_col]],
      .data[[antibiotic_class_col]]
    ) %>%
    dplyr::summarise(
      class_R = as.integer(any(.abx_val == "R")),
      .groups = "drop"
    )

  orig_classes <- sort(unique(class_long[[antibiotic_class_col]]))
  safe_classes <- make.names(orig_classes)
  class_name_map <- setNames(orig_classes, safe_classes)

  class_wide <- class_long %>%
    dplyr::mutate(
      .class_safe = make.names(.data[[antibiotic_class_col]])
    ) %>%
    tidyr::pivot_wider(
      id_cols     = !!rlang::sym(patient_id_col),
      names_from  = ".class_safe",
      values_from = "class_R",
      values_fill = untested_fill
    )

  attr(class_wide, "class_name_map") <- class_name_map
  return(class_wide)
}


# -- derive_infection_type_for_mortality ---------------------------------------

#' Derive Infection Type (HAI / CAI) for Mortality RR Model
#'
#' Variant of \code{derive_infection_type()} tailored for the mortality
#' relative-risk model. Processes all patients (dead and surviving), recognises
#' explicit HAI/CAI labels first, and falls back to date-gap derivation when
#' labels are ambiguous. Also performs a death-date data-quality check.
#'
#' @param data Data frame.
#' @param infection_type_col Character. Default \code{"type_of_infection"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_culture_col Character. Default \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_date_col Character. Default \code{"final_outcome_date"}.
#' @param death_value Character. Default \code{"Death"}.
#' @param hai_threshold_hours Numeric. Default \code{48}.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#'
#' @return \code{data} with column \code{infection_type_derived}
#'   (\code{"HAI"} / \code{"CAI"} / \code{"Not Known"}).
#' @export
derive_infection_type_for_mortality <- function(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  death_value = "Death",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
) {
  required <- c(
    infection_type_col, date_admission_col,
    date_culture_col, final_outcome_col
  )
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  has_outcome_date <- final_outcome_date_col %in% names(data)

  if (has_outcome_date) {
    is_dead <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[final_outcome_col]])
    )) == stringr::str_to_upper(death_value)

    missing_outcome_dt <- is.na(suppressWarnings(
      as.Date(as.character(data[[final_outcome_date_col]]))
    ))
    has_admit_or_culture <- !is.na(data[[date_admission_col]]) |
      !is.na(data[[date_culture_col]])
    dead_no_date <- is_dead & missing_outcome_dt & has_admit_or_culture

    if (any(dead_no_date)) {
      n_flag <- sum(dead_no_date)
      message(sprintf(
        paste0(
          "[Mortality] Data-quality flag: %d patient row(s) have ",
          "final_outcome='%s' but no outcome date recorded. ",
          "HAI/CAI derivation is unaffected."
        ),
        n_flag, death_value
      ))
      if (patient_id_col %in% names(data)) {
        flagged_ids <- unique(data[dead_no_date, patient_id_col, drop = TRUE])
        message(sprintf(
          "  Flagged patient IDs (first 20 of %d): %s",
          length(flagged_ids),
          paste(head(flagged_ids, 20L), collapse = ", ")
        ))
      }
      flagged_df <- data[dead_no_date, , drop = FALSE]
    } else {
      flagged_df <- data[0L, , drop = FALSE]
    }
  } else {
    message(sprintf(
      "[Mortality] Column '%s' not found; skipping death-date quality check.",
      final_outcome_date_col
    ))
    flagged_df <- data[0L, , drop = FALSE]
  }

  ambiguous_mask <- {
    inf_up <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[infection_type_col]])
    ))
    inf_up %in% c("NOT KNOWN", "NOT_KNOWN", "UNKNOWN", "NULL", "NA", "") |
      is.na(data[[infection_type_col]])
  }
  missing_admit <- is.na(data[[date_admission_col]])
  missing_culture <- is.na(data[[date_culture_col]])
  cannot_infer <- ambiguous_mask & (missing_admit | missing_culture)

  if (any(cannot_infer)) {
    message(sprintf(
      "[Mortality] Cannot infer HAI/CAI for %d row(s): assigned 'Not Known'.",
      sum(cannot_infer)
    ))
  }

  data <- data %>%
    dplyr::mutate(
      .inf_raw = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[infection_type_col]])
      )),
      .gap_h = as.numeric(difftime(
        as.Date(suppressWarnings(as.character(.data[[date_culture_col]]))),
        as.Date(suppressWarnings(as.character(.data[[date_admission_col]]))),
        units = "hours"
      )),
      .cannot_infer = (.inf_raw %in% c(
        "NOT KNOWN", "NOT_KNOWN",
        "UNKNOWN", "NULL", "NA", ""
      ) |
        is.na(.data[[infection_type_col]])),
      infection_type_derived = dplyr::case_when(
        .inf_raw %in% c(
          "HAI", "HOSPITAL ACQUIRED", "HOSPITAL-ACQUIRED",
          "HOSPITAL_ACQUIRED", "NOSOCOMIAL",
          "HOSPITAL ACQUIRED INFECTION",
          "HOSPITAL-ACQUIRED INFECTION"
        ) ~ "HAI",
        .inf_raw %in% c(
          "CAI", "COMMUNITY ACQUIRED",
          "COMMUNITY-ACQUIRED", "COMMUNITY_ACQUIRED",
          "COMMUNITY ACQUIRED INFECTION",
          "COMMUNITY-ACQUIRED INFECTION"
        ) ~ "CAI",
        .cannot_infer &
          !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours, "CAI", "HAI"),
        TRUE ~ "Not Known"
      )
    ) %>%
    dplyr::select(-".inf_raw", -".gap_h", -".cannot_infer")

  n_hai <- sum(data$infection_type_derived == "HAI", na.rm = TRUE)
  n_cai <- sum(data$infection_type_derived == "CAI", na.rm = TRUE)
  n_not_known <- sum(data$infection_type_derived == "Not Known", na.rm = TRUE)
  message(sprintf(
    "[Mortality] Infection type: HAI = %d | CAI = %d | Not Known = %d.",
    n_hai, n_cai, n_not_known
  ))

  attr(data, "missing_death_date_patients") <- flagged_df
  return(data)
}


# -- .derive_icu_binary --------------------------------------------------------

#' Derive ICU Binary Flag per Patient
#'
#' Collapses unit-type data to a patient-level binary ICU indicator using the
#' "ever in ICU" rule. Returns integer 0/1 or a 3-level factor when missing
#' data is high.
#'
#' @param data Data frame at drug-test level.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param unit_type_col Character. Default \code{"unit_type"}.
#' @param icu_values Character vector. ICU location labels.
#' @param missing_threshold Numeric. Proportion above which a 3-level factor
#'   is returned. Default \code{0.10}.
#'
#' @return Patient-level data frame with \code{patient_id_col} and \code{ICU}.
#' @keywords internal
.derive_icu_binary <- function(
  data,
  patient_id_col = "PatientInformation_id",
  unit_type_col = "unit_type",
  icu_values = c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU"),
  missing_threshold = 0.10
) {
  if (!unit_type_col %in% names(data)) {
    message(sprintf(
      "[ICU] Column '%s' not found -- ICU covariate will be omitted.",
      unit_type_col
    ))
    return(
      dplyr::distinct(data, !!rlang::sym(patient_id_col)) %>%
        dplyr::mutate(ICU = NA_integer_)
    )
  }

  icu_pattern <- paste(stringr::str_to_lower(icu_values), collapse = "|")

  per_patient <- data %>%
    dplyr::mutate(
      .unit_norm = stringr::str_to_lower(stringr::str_trim(
        as.character(.data[[unit_type_col]])
      )),
      .is_missing = is.na(.data[[unit_type_col]]) |
        .unit_norm %in% c("na", "", "null"),
      .is_icu = !.is_missing &
        grepl(icu_pattern, .unit_norm, fixed = FALSE)
    ) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::summarise(
      .any_icu     = any(.is_icu, na.rm = TRUE),
      .all_missing = all(.is_missing),
      .groups      = "drop"
    )

  prop_missing <- sum(per_patient$.all_missing) / nrow(per_patient)

  if (prop_missing > missing_threshold) {
    message(sprintf(
      "[ICU] %.1f%% of patients have entirely missing unit_type. Returning 3-level factor.",
      prop_missing * 100
    ))
    per_patient <- per_patient %>%
      dplyr::mutate(
        ICU = factor(
          dplyr::case_when(
            .all_missing ~ "Unknown",
            .any_icu ~ "ICU",
            TRUE ~ "Ward"
          ),
          levels = c("Ward", "ICU", "Unknown")
        )
      )
  } else {
    per_patient <- per_patient %>%
      dplyr::mutate(
        ICU = dplyr::case_when(
          .all_missing ~ NA_integer_,
          .any_icu ~ 1L,
          TRUE ~ 0L
        )
      )
  }

  n_icu <- sum(per_patient$ICU == "ICU" | per_patient$ICU == 1L, na.rm = TRUE)
  n_ward <- sum(per_patient$ICU == "Ward" | per_patient$ICU == 0L, na.rm = TRUE)
  n_unk <- sum(is.na(per_patient$ICU) | per_patient$ICU == "Unknown", na.rm = TRUE)
  message(sprintf(
    "[ICU] Per-patient: ICU = %d | Ward = %d | Unknown/NA = %d.",
    n_icu, n_ward, n_unk
  ))

  per_patient %>% dplyr::select(!!rlang::sym(patient_id_col), ICU)
}


# -- .encode_comorbidity_mortality ---------------------------------------------

#' Encode Comorbidity Column for Mortality Model
#'
#' Standardises a comorbidity column to numeric, binary (0/1), or ordinal
#' factor encoding for use in \code{fit_mortality_rr_logistic()}.
#'
#' @param data Patient-level data frame.
#' @param comorbidity_col Character. Default \code{"comorbidities"}.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#'
#' @return \code{data} with \code{comorbidity_encoded} added. Attribute
#'   \code{"comorbidity_encoding"} records the strategy used.
#' @keywords internal
.encode_comorbidity_mortality <- function(
  data,
  comorbidity_col = "comorbidities",
  patient_id_col = "PatientInformation_id"
) {
  if (!comorbidity_col %in% names(data)) {
    message(sprintf(
      "[Comorbidity] Column '%s' not found -- covariate will be omitted.",
      comorbidity_col
    ))
    data$comorbidity_encoded <- NA_real_
    attr(data, "comorbidity_encoding") <- "absent"
    return(data)
  }

  raw <- data[[comorbidity_col]]

  if (is.numeric(raw)) {
    data$comorbidity_encoded <- raw
    attr(data, "comorbidity_encoding") <- "numeric"
    message(sprintf(
      "[Comorbidity] Numeric index detected (range %.1f-%.1f). Used as-is.",
      min(raw, na.rm = TRUE), max(raw, na.rm = TRUE)
    ))
    return(data)
  }

  raw_up <- stringr::str_to_upper(stringr::str_trim(as.character(raw)))
  na_vals <- c(
    "NA", "NULL", "", "UNKNOWN", "NOT KNOWN", "NOT_KNOWN",
    "MISSING", "N/A", "NONE RECORDED"
  )
  raw_up[raw_up %in% na_vals | is.na(raw)] <- NA_character_

  binary_present <- c(
    "PRESENT", "YES", "1", "TRUE", "POSITIVE",
    "COMORBID", "COMORBIDITY PRESENT"
  )
  binary_none <- c(
    "NONE", "NO", "0", "FALSE", "ABSENT",
    "NO COMORBIDITY", "COMORBIDITY ABSENT"
  )
  non_na_vals <- raw_up[!is.na(raw_up)]

  if (length(non_na_vals) > 0L &&
    all(non_na_vals %in% c(binary_present, binary_none))) {
    data$comorbidity_encoded <- dplyr::case_when(
      raw_up %in% binary_present ~ 1L,
      raw_up %in% binary_none ~ 0L,
      TRUE ~ NA_integer_
    )
    attr(data, "comorbidity_encoding") <- "binary"
    message(sprintf("[Comorbidity] Binary encoding applied."))
    return(data)
  }

  ord_none <- c("NONE", "NO COMORBIDITY", "ABSENT", "0")
  ord_mild <- c("MILD", "LOW", "MINOR", "MINIMAL", "1")
  ord_moderate <- c("MODERATE", "MEDIUM", "MODERATE COMORBIDITY", "2")
  ord_severe <- c(
    "SEVERE", "HIGH", "MAJOR", "HEAVY", "SIGNIFICANT",
    "MULTIPLE", "3"
  )

  data$comorbidity_encoded <- factor(
    dplyr::case_when(
      raw_up %in% ord_none ~ "none",
      raw_up %in% ord_mild ~ "mild",
      raw_up %in% ord_moderate ~ "moderate",
      raw_up %in% ord_severe ~ "severe",
      TRUE ~ NA_character_
    ),
    levels = c("none", "mild", "moderate", "severe"), ordered = TRUE
  )
  attr(data, "comorbidity_encoding") <- "ordinal"
  message("[Comorbidity] Ordinal encoding applied.")
  return(data)
}


# -- .check_hai_icu_collinearity -----------------------------------------------

#' Check Collinearity Between HAI and ICU Covariates
#'
#' Computes the phi correlation for the HAI x ICU 2x2 table and warns when
#' it exceeds \code{phi_threshold}.
#'
#' @param df Patient-level data frame with integer 0/1 columns.
#' @param hai_col Character. Default \code{"HAI"}.
#' @param icu_col Character. Default \code{"ICU"}.
#' @param phi_threshold Numeric. Default \code{0.7}.
#'
#' @return Named list: \code{phi}, \code{tbl}, \code{warning_issued}.
#' @keywords internal
.check_hai_icu_collinearity <- function(df,
                                        hai_col = "HAI",
                                        icu_col = "ICU",
                                        phi_threshold = 0.7) {
  if (!hai_col %in% names(df) || !icu_col %in% names(df)) {
    return(list(phi = NA_real_, tbl = NULL, warning_issued = FALSE))
  }

  h <- as.integer(df[[hai_col]])
  ic <- as.integer(df[[icu_col]])
  ok <- !is.na(h) & !is.na(ic)
  h <- h[ok]
  ic <- ic[ok]

  if (length(h) < 10L) {
    message("[Collinearity] Too few observations to assess HAI/ICU correlation.")
    return(list(phi = NA_real_, tbl = NULL, warning_issued = FALSE))
  }

  tbl <- table(HAI = h, ICU = ic)
  message("[Collinearity] HAI x ICU cross-table:")
  print(tbl)

  if (any(rowSums(tbl) == 0L) || any(colSums(tbl) == 0L)) {
    warning("[Collinearity] Perfect separation detected between HAI and ICU.")
    return(list(phi = NA_real_, tbl = tbl, warning_issued = TRUE))
  }

  phi <- tryCatch(stats::cor(h, ic, method = "pearson"), error = function(e) NA_real_)
  warn_issued <- FALSE

  if (!is.na(phi) && abs(phi) >= phi_threshold) {
    warning(sprintf(
      "[Collinearity] HAI and ICU are highly correlated (phi = %.3f >= %.2f).",
      phi, phi_threshold
    ))
    warn_issued <- TRUE
  } else if (!is.na(phi)) {
    message(sprintf("[Collinearity] HAI/ICU phi = %.3f (OK).", phi))
  }

  list(phi = phi, tbl = tbl, warning_issued = warn_issued)
}


# -- fit_mortality_rr_logistic -------------------------------------------------

#' Fit Mixed Logistic Models to Estimate Mortality OR per Pathogen x Class
#'
#' For each pathogen x antibiotic class combination, fits a mixed-effects
#' logistic regression (facility random intercept) with resistance, age, sex,
#' and HAI as fixed effects. Returns a data frame of mortality odds ratios (OR)
#' and 95% Wald CIs. Use \code{convert_or_to_rr()} to convert ORs to RRs.
#'
#' @param data Data frame. Patient-level microbiology records.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param facility_col Character. Default \code{"center_name"}.
#' @param organism_col Character. Default \code{"organism_name"}.
#' @param syndrome_col Character. Default \code{"syndrome"}.
#' @param infection_type_col Character. Default \code{"type_of_infection"}.
#' @param antibiotic_class_col Character. Default \code{"antibiotic_class"}.
#' @param antibiotic_name_col Character. Default \code{"antibiotic_name"}.
#' @param antibiotic_value_col Character. Default \code{"antibiotic_value"}.
#' @param unit_type_col Character. Default \code{"unit_type"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_culture_col Character. Default \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_date_col Character. Default \code{"final_outcome_date"}.
#' @param age_col Character. Default \code{"Age"}.
#' @param sex_col Character. Default \code{"Gender"}.
#' @param comorbidity_col Character or \code{NULL}. Default \code{NULL}.
#' @param death_value Character. Default \code{"Death"}.
#' @param syndrome_name Character or \code{NULL}. Filter to one syndrome.
#' @param organism_name Character vector or \code{NULL}. Filter to pathogen(s).
#' @param facility_name Character or \code{NULL}. If provided, filters data
#'   to the specified facility before fitting.
#' @param hai_threshold_hours Numeric. Default \code{48}.
#' @param icu_values Character vector. ICU location labels.
#' @param phi_threshold Numeric. HAI/ICU collinearity threshold. Default \code{0.7}.
#' @param min_n Integer. Minimum patients per model. Default \code{10L}.
#' @param min_deaths Integer. Minimum deaths per model. Default \code{5L}.
#'
#' @return Data frame with columns: \code{pathogen}, \code{antibiotic_class},
#'   \code{OR_death}, \code{CI_lower}, \code{CI_upper}, \code{n_patients},
#'   \code{n_deaths}, \code{convergence_warning}, \code{syndrome_scope},
#'   \code{comorbidity_encoding}.
#' @export
fit_mortality_rr_logistic <- function(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "type_of_infection",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  unit_type_col = "unit_type",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  age_col = "Age",
  sex_col = "Gender",
  comorbidity_col = NULL,
  death_value = "Death",
  syndrome_name = NULL,
  organism_name = NULL,
  facility_name = NULL,
  hai_threshold_hours = 48,
  icu_values = c(
    "ICU", "Intensive Care", "Critical Care",
    "PICU", "NICU"
  ),
  phi_threshold = 0.7,
  min_n = 10L,
  min_deaths = 5L
) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required. Install with: install.packages('lme4')")
  }

  required <- c(
    patient_id_col, facility_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_culture_col,
    final_outcome_col, age_col, sex_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  use_syndrome_covariate <- is.null(syndrome_name) && syndrome_col %in% names(data)
  missing_req <- setdiff(required, names(data))
  if (length(missing_req) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_req, collapse = ", ")
    ))
  }

  data <- derive_infection_type_for_mortality(
    data,
    infection_type_col     = infection_type_col,
    date_admission_col     = date_admission_col,
    date_culture_col       = date_culture_col,
    final_outcome_col      = final_outcome_col,
    final_outcome_date_col = final_outcome_date_col,
    death_value            = death_value,
    hai_threshold_hours    = hai_threshold_hours,
    patient_id_col         = patient_id_col
  )

  df <- data
  if (!is.null(syndrome_name)) {
    df <- dplyr::filter(df, .data[[syndrome_col]] == syndrome_name)
  }
  if (!is.null(facility_name)) {
    df <- dplyr::filter(df, .data[[facility_col]] == facility_name)
  }

  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(df[[organism_col]]))
  }

  icu_tbl <- .derive_icu_binary(df,
    patient_id_col = patient_id_col,
    unit_type_col = unit_type_col,
    icu_values = icu_values
  )
  use_icu <- unit_type_col %in% names(df) && !all(is.na(icu_tbl$ICU))

  patient_covars <- df %>%
    dplyr::mutate(
      death = as.integer(
        stringr::str_to_upper(stringr::str_trim(
          as.character(.data[[final_outcome_col]])
        )) == stringr::str_to_upper(death_value)
      ),
      HAI = as.integer(infection_type_derived == "HAI")
    ) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      !!rlang::sym(patient_id_col), !!rlang::sym(facility_col),
      !!rlang::sym(organism_col), !!rlang::sym(age_col),
      !!rlang::sym(sex_col), death, HAI, infection_type_derived
    ) %>%
    dplyr::mutate(
      !!rlang::sym(age_col) := suppressWarnings(as.numeric(.data[[age_col]])),
      !!rlang::sym(sex_col) := factor(.data[[sex_col]])
    )

  if (use_syndrome_covariate) {
    syndrome_raw <- df %>%
      dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!rlang::sym(patient_id_col), !!rlang::sym(syndrome_col)) %>%
      dplyr::mutate(!!rlang::sym(syndrome_col) := factor(.data[[syndrome_col]]))
    patient_covars <- dplyr::left_join(patient_covars, syndrome_raw,
      by = patient_id_col
    )
    message(sprintf(
      "[Syndrome] '%s' added as covariate (%d levels).",
      syndrome_col, nlevels(patient_covars[[syndrome_col]])
    ))
  }

  use_comorbidity <- FALSE
  comorbidity_encoding <- "absent"

  if (!is.null(comorbidity_col) && comorbidity_col %in% names(df)) {
    comorbid_raw <- df %>%
      dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!rlang::sym(patient_id_col), !!rlang::sym(comorbidity_col))
    comorbid_enc <- .encode_comorbidity_mortality(comorbid_raw,
      comorbidity_col = comorbidity_col,
      patient_id_col  = patient_id_col
    )
    comorbidity_encoding <- attr(comorbid_enc, "comorbidity_encoding")
    patient_covars <- dplyr::left_join(
      patient_covars,
      dplyr::select(comorbid_enc, !!rlang::sym(patient_id_col), comorbidity_encoded),
      by = patient_id_col
    )
    use_comorbidity <- comorbidity_encoding != "absent" &&
      !all(is.na(patient_covars$comorbidity_encoded))
  } else {
    patient_covars$comorbidity_encoded <- NA_real_
  }

  patient_covars <- dplyr::left_join(patient_covars, icu_tbl, by = patient_id_col)

  all_or <- list()

  for (path in pathogens) {
    path_ids <- df %>%
      dplyr::filter(stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
        stringr::str_to_lower(stringr::str_trim(path))) %>%
      dplyr::distinct(!!rlang::sym(patient_id_col)) %>%
      dplyr::pull(!!rlang::sym(patient_id_col))

    if (length(path_ids) == 0L) {
      message(sprintf("'%s': no patients, skipping.", path))
      next
    }

    path_df_raw <- df %>%
      dplyr::filter(stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
        stringr::str_to_lower(stringr::str_trim(path)))

    resist_wide <- .build_class_resistance_wide(
      path_df_raw,
      patient_id_col = patient_id_col,
      antibiotic_class_col = antibiotic_class_col,
      antibiotic_name_col = antibiotic_name_col,
      antibiotic_value_col = antibiotic_value_col
    )
    class_name_map <- attr(resist_wide, "class_name_map")
    class_safe <- setdiff(names(resist_wide), patient_id_col)

    model_data <- patient_covars %>%
      dplyr::filter(.data[[patient_id_col]] %in% path_ids) %>%
      dplyr::inner_join(resist_wide, by = patient_id_col) %>%
      dplyr::filter(
        infection_type_derived != "Not Known",
        !is.na(.data[[age_col]]), !is.na(.data[[sex_col]]),
        !is.na(HAI)
      )
    if (use_syndrome_covariate) {
      model_data <- dplyr::filter(model_data, !is.na(.data[[syndrome_col]]))
    }

    if (nrow(model_data) < min_n) {
      message(sprintf("'%s': %d patients < min_n=%d, skipping.", path, nrow(model_data), min_n))
      next
    }
    if (sum(model_data$death, na.rm = TRUE) < min_deaths) {
      message(sprintf(
        "'%s': %d deaths < min_deaths=%d, skipping.",
        path, sum(model_data$death, na.rm = TRUE), min_deaths
      ))
      next
    }

    if (use_icu) {
      invisible(.check_hai_icu_collinearity(model_data, phi_threshold = phi_threshold))
    }

    n_facilities <- length(unique(model_data[[facility_col]]))

    for (cls in class_safe) {
      orig_class <- class_name_map[[cls]]
      sub <- model_data %>% dplyr::filter(!is.na(.data[[cls]]))

      if (length(unique(stats::na.omit(sub[[cls]]))) < 2L) {
        message(sprintf("'%s' | '%s': no resistance variation, skipping.", path, orig_class))
        next
      }
      if (nrow(sub) < min_n) next
      if (sum(sub$death) < min_deaths || sum(sub$death == 0L) < min_deaths) next

      fixed_terms <- c(
        sprintf("`%s`", cls), sprintf("`%s`", age_col),
        sprintf("`%s`", sex_col), "HAI"
      )
      if (use_icu && length(unique(stats::na.omit(as.character(sub$ICU)))) > 1L) {
        fixed_terms <- c(fixed_terms, "ICU")
      }
      if (use_comorbidity && !all(is.na(sub$comorbidity_encoded)) &&
        length(unique(stats::na.omit(as.character(sub$comorbidity_encoded)))) > 1L) {
        fixed_terms <- c(fixed_terms, "comorbidity_encoded")
      }
      if (use_syndrome_covariate &&
        length(unique(stats::na.omit(as.character(sub[[syndrome_col]])))) > 1L) {
        fixed_terms <- c(fixed_terms, sprintf("`%s`", syndrome_col))
      }

      needed_cols <- unique(gsub("`", "", c("death", fixed_terms)))
      sub_fit <- sub %>%
        dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(needed_cols))))

      if (nrow(sub_fit) < min_n) next
      n_deaths_fit <- sum(sub_fit$death, na.rm = TRUE)
      n_surv_fit <- sum(sub_fit$death == 0L, na.rm = TRUE)
      if (n_deaths_fit < min_deaths || n_surv_fit < min_deaths) next

      epv <- n_deaths_fit / length(fixed_terms)
      if (epv < 10) {
        message(sprintf(
          "'%s' | '%s': EPV = %.1f < 10 -- estimates may be unreliable.",
          path, orig_class, epv
        ))
      }

      ct_sep <- table(sub_fit[[cls]], sub_fit$death)
      if (any(ct_sep == 0L)) {
        message(sprintf(
          "'%s' | '%s': zero cell (complete separation), skipping.",
          path, orig_class
        ))
        next
      }

      re_term <- if (n_facilities > 1L) sprintf("(1 | `%s`)", facility_col) else NULL
      fmla <- stats::as.formula(sprintf(
        "death ~ %s%s",
        paste(fixed_terms, collapse = " + "),
        if (!is.null(re_term)) paste0(" + ", re_term) else ""
      ))

      conv_warn <- FALSE
      model <- tryCatch(
        {
          fit <- suppressWarnings(lme4::glmer(
            fmla,
            data = sub_fit, family = stats::binomial,
            control = lme4::glmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 2e5)
            )
          ))
          if (lme4::isSingular(fit)) {
            message(sprintf("'%s' | '%s': singular fit.", path, orig_class))
            conv_warn <- TRUE
          }
          conv_msgs <- fit@optinfo$conv$lme4$messages
          if (!is.null(conv_msgs) && length(conv_msgs) > 0L) {
            message(sprintf(
              "'%s' | '%s': convergence warning -- %s",
              path, orig_class, paste(conv_msgs, collapse = "; ")
            ))
            conv_warn <- TRUE
          }
          fit
        },
        error = function(e) {
          message(sprintf(
            "'%s' | '%s': model failed -- %s",
            path, orig_class, conditionMessage(e)
          ))
          NULL
        }
      )
      if (is.null(model)) next

      coefs <- lme4::fixef(model)
      vcmat <- stats::vcov(model)
      pat <- paste0("^`?", gsub(".", "\\.", cls, fixed = TRUE), "`?$")
      cname <- grep(pat, names(coefs), value = TRUE)
      if (length(cname) == 0L) next

      beta <- coefs[[cname[1L]]]
      se <- sqrt(vcmat[cname[1L], cname[1L]])

      if (requireNamespace("car", quietly = TRUE)) {
        vif_vals <- tryCatch(car::vif(model), error = function(e) NULL)
        if (!is.null(vif_vals)) {
          vif_num <- if (is.matrix(vif_vals)) vif_vals[, "GVIF"] else as.numeric(vif_vals)
          high_vif <- names(vif_num)[vif_num > 10]
          if (length(high_vif) > 0) {
            message(sprintf(
              "'%s' | '%s': high VIF (>10) for: %s",
              path, orig_class, paste(high_vif, collapse = ", ")
            ))
          }
        }
      }

      all_or[[length(all_or) + 1L]] <- data.frame(
        pathogen             = path,
        antibiotic_class     = orig_class,
        OR_death             = round(exp(beta), 4L),
        CI_lower             = round(exp(beta - 1.96 * se), 4L),
        CI_upper             = round(exp(beta + 1.96 * se), 4L),
        n_patients           = nrow(sub_fit),
        n_deaths             = n_deaths_fit,
        convergence_warning  = conv_warn,
        syndrome_scope       = if (!is.null(syndrome_name)) syndrome_name else if (use_syndrome_covariate) "covariate" else "all",
        comorbidity_encoding = comorbidity_encoding,
        stringsAsFactors     = FALSE
      )
    }

    n_cls_fitted <- sum(sapply(all_or, function(x) x$pathogen == path))
    message(sprintf("'%s': %d class OR(s) estimated.", path, n_cls_fitted))
  }

  if (length(all_or) == 0L) {
    warning("No mortality OR values could be estimated.")
    return(data.frame())
  }
  dplyr::bind_rows(all_or)
}


# -- compute_p0 ----------------------------------------------------------------

#' Compute Baseline Mortality Rate Among Fully Susceptible Patients (p0)
#'
#' Computes the baseline death probability \eqn{p_0} among patients who are
#' fully susceptible -- i.e. \strong{all} antibiotic test results are
#' \code{"S"}. Any patient with at least one \code{"R"} result is excluded.
#' Used as the denominator correction for OR -> RR conversion (Zhang & Yu 1998).
#'
#' @param data               Data frame. One row per patient x antibiotic test.
#' @param patient_col        Character. Unique patient identifier column.
#' @param antibiotic_value_col Character. Antibiotic result column
#'   (\code{"S"}, \code{"R"}, \code{"I"}).
#' @param outcome_col        Character. Final outcome column.
#' @param death_value        Character. Death indicator. Default \code{"Death"}.
#' @param resistant_value    Character. Resistance indicator used to exclude
#'   patients. Default \code{"R"}.
#' @param syndrome_col       Character or \code{NULL}.
#' @param syndrome_name      Character or \code{NULL}.
#' @param facility_col       Character or \code{NULL}. Facility identifier
#'   column. Required when \code{facility_name} is specified.
#' @param facility_name      Character or \code{NULL}. If provided, filters
#'   data to the specified facility before computing p0.
#'
#' @return Named list: \code{p0} (scalar), \code{n_susceptible},
#'   \code{n_susceptible_deaths}, \code{summary} (one-row data frame).
#' @export
#'
#' @references
#' Zhang J, Yu KF. What's the relative risk? JAMA. 1998;280(19):1690-1.
compute_p0 <- function(data,
                       patient_col,
                       antibiotic_value_col,
                       outcome_col,
                       death_value = "Death",
                       resistant_value = "R",
                       syndrome_col = NULL,
                       syndrome_name = NULL,
                       facility_col = NULL,
                       facility_name = NULL) {
  required_cols <- c(patient_col, antibiotic_value_col, outcome_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("syndrome_col must be supplied when syndrome_name is specified.")
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be supplied when facility_name is specified.")
  }

  df <- data
  if (!is.null(syndrome_name)) {
    df <- df %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
  }
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
  }
  if (nrow(df) == 0) stop("No records remain after filtering.")

  # Flag patients with at least one R result -- they are excluded
  patient_susceptibility <- df %>%
    dplyr::group_by(.data[[patient_col]]) %>%
    dplyr::summarise(
      has_resistant = any(.data[[antibiotic_value_col]] == resistant_value,
        na.rm = TRUE
      ),
      outcome = dplyr::first(.data[[outcome_col]]),
      .groups = "drop"
    ) %>%
    dplyr::filter(!has_resistant)

  n_susceptible <- nrow(patient_susceptibility)
  n_susceptible_deaths <- sum(patient_susceptibility$outcome == death_value,
    na.rm = TRUE
  )

  if (n_susceptible == 0) {
    stop("No fully susceptible patients found. Check antibiotic_value_col and resistant_value.")
  }

  p0 <- n_susceptible_deaths / n_susceptible

  message(sprintf(
    "p0 = %.4f | susceptible patients: %d | susceptible deaths: %d%s%s",
    p0, n_susceptible, n_susceptible_deaths,
    if (!is.null(syndrome_name)) sprintf(" [%s]", syndrome_name) else "",
    if (!is.null(facility_name)) sprintf(" [%s]", facility_name) else ""
  ))

  list(
    p0 = p0,
    n_susceptible = n_susceptible,
    n_susceptible_deaths = n_susceptible_deaths,
    summary = data.frame(
      syndrome             = if (!is.null(syndrome_name)) syndrome_name else "all",
      facility_name        = if (!is.null(facility_name)) facility_name else NA_character_,
      p0                   = p0,
      n_susceptible        = n_susceptible,
      n_susceptible_deaths = n_susceptible_deaths,
      stringsAsFactors     = FALSE
    )
  )
}


# -- convert_or_to_rr ---------------------------------------------------------

#' Convert Odds Ratios to Relative Risks Using Baseline Mortality (p0)
#'
#' Applies the Zhang & Yu (1998) formula to convert odds ratios from logistic
#' regression to approximate relative risks:
#' \deqn{RR = \frac{OR}{(1 - p_0) + p_0 \times OR}}
#' Applied to the point estimate and both CI bounds. Works directly on the
#' output of \code{fit_mortality_rr_logistic()}.
#'
#' @param or_data      Data frame with OR columns.
#' @param p0           Numeric scalar. From \code{compute_p0()$p0}.
#' @param or_col       Character. Default \code{"OR_death"}.
#' @param ci_lower_col Character. Default \code{"CI_lower"}.
#' @param ci_upper_col Character. Default \code{"CI_upper"}.
#'
#' @return \code{or_data} with columns \code{RR_death}, \code{RR_lower},
#'   \code{RR_upper} appended.
#' @export
#'
#' @references
#' Zhang J, Yu KF. What's the relative risk? JAMA. 1998;280(19):1690-1.
convert_or_to_rr <- function(or_data,
                             p0,
                             or_col = "OR_death",
                             ci_lower_col = "CI_lower",
                             ci_upper_col = "CI_upper") {
  if (!is.numeric(p0) || length(p0) != 1L || is.na(p0) || p0 < 0 || p0 > 1) {
    stop("'p0' must be a single numeric value between 0 and 1.")
  }

  missing_cols <- setdiff(c(or_col, ci_lower_col, ci_upper_col), names(or_data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in or_data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  .to_rr <- function(or) or / ((1 - p0) + (p0 * or))

  result <- or_data %>%
    dplyr::mutate(
      RR_death = round(.to_rr(.data[[or_col]]), 4L),
      RR_lower = round(.to_rr(.data[[ci_lower_col]]), 4L),
      RR_upper = round(.to_rr(.data[[ci_upper_col]]), 4L)
    )

  message(sprintf(
    "OR -> RR conversion: p0 = %.4f | %d rows | RR range [%.3f, %.3f]",
    p0, nrow(result),
    min(result$RR_death, na.rm = TRUE),
    max(result$RR_death, na.rm = TRUE)
  ))

  result
}


# -- assign_rr_to_profiles ----------------------------------------------------

#' Assign Per-Class RR to Resistance Profiles (Max Rule)
#'
#' For each resistance profile delta, determines the profile-level RR using
#' the GBD max rule: RR_kd = max over resistant classes of RR_kc.
#' The CI reported for each profile is that of its dominant (max-RR) class.
#' Also used for mortality ORs when called with \code{rr_col = "OR_death"}.
#'
#' @param profiles_output Named list from \code{compute_resistance_profiles()}.
#' @param rr_table Data frame from \code{fit_mortality_rr_logistic()} or a
#'   LOS RR table. Must have \code{pathogen_col}, \code{class_col},
#'   \code{rr_col}, and optionally \code{CI_lower} / \code{CI_upper}.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param rr_col Character. Default \code{"RR_LOS"}. Use \code{"OR_death"}
#'   for mortality.
#' @param fallback_rr Numeric. RR for classes with no match. Default \code{1}.
#'
#' @return Named list (one per pathogen): profiles data frame augmented with
#'   \code{RR_LOS_profile}, \code{dominant_class}, and CI columns.
#' @export
assign_rr_to_profiles <- function(
  profiles_output,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rr_col = "RR_LOS",
  fallback_rr = 1
) {
  if (!is.list(profiles_output) ||
    !all(sapply(
      profiles_output,
      function(x) all(c("profiles", "classes") %in% names(x))
    ))) {
    stop("profiles_output must be the list returned by compute_resistance_profiles().")
  }

  missing_rr <- setdiff(c(pathogen_col, class_col, rr_col), names(rr_table))
  if (length(missing_rr) > 0L) {
    stop(sprintf(
      "Column(s) not found in rr_table: %s",
      paste(missing_rr, collapse = ", ")
    ))
  }

  has_ci <- all(c("CI_lower", "CI_upper") %in% names(rr_table))
  out <- list()

  for (path in names(profiles_output)) {
    prof_list <- profiles_output[[path]]
    profiles <- prof_list$profiles
    classes <- prof_list$classes

    rr_k <- rr_table[rr_table[[pathogen_col]] == path, ]
    rr_lookup <- setNames(rr_k[[rr_col]], rr_k[[class_col]])
    ci_lo_lookup <- if (has_ci) setNames(rr_k$CI_lower, rr_k[[class_col]]) else NULL
    ci_hi_lookup <- if (has_ci) setNames(rr_k$CI_upper, rr_k[[class_col]]) else NULL

    n_prof <- nrow(profiles)
    rr_profile <- numeric(n_prof)
    dom_class <- character(n_prof)
    ci_lo_prof <- if (has_ci) numeric(n_prof) else NULL
    ci_hi_prof <- if (has_ci) numeric(n_prof) else NULL

    for (i in seq_len(n_prof)) {
      resist_cls <- classes[as.integer(profiles[i, classes]) == 1L]
      if (length(resist_cls) == 0L) {
        rr_profile[i] <- 1.0
        dom_class[i] <- "all_susceptible"
        if (has_ci) {
          ci_lo_prof[i] <- 1.0
          ci_hi_prof[i] <- 1.0
        }
        next
      }
      rrs_c <- sapply(resist_cls, function(cc) {
        ifelse(cc %in% names(rr_lookup), rr_lookup[[cc]], fallback_rr)
      })
      rrs_c[is.na(rrs_c) | is.nan(rrs_c)] <- fallback_rr
      max_idx <- which.max(rrs_c)
      rr_profile[i] <- rrs_c[max_idx]
      dom_class[i] <- resist_cls[max_idx]
      if (has_ci) {
        dc <- dom_class[i]
        ci_lo_prof[i] <- ifelse(dc %in% names(ci_lo_lookup), ci_lo_lookup[[dc]], fallback_rr)
        ci_hi_prof[i] <- ifelse(dc %in% names(ci_hi_lookup), ci_hi_lookup[[dc]], fallback_rr)
      }
    }

    # Name the profile-level RR column after the source column
    # e.g. rr_col="RR_LOS"   -> "RR_LOS_profile"
    #      rr_col="RR_death" -> "RR_death_profile"
    rr_profile_col_name <- paste0(rr_col, "_profile")
    profiles[[rr_profile_col_name]] <- round(rr_profile, 4L)
    profiles$dominant_class <- dom_class
    if (has_ci) {
      profiles$CI_lower_profile <- round(ci_lo_prof, 4L)
      profiles$CI_upper_profile <- round(ci_hi_prof, 4L)
    }
    out[[path]] <- profiles
    message(sprintf(
      "'%s': RR assigned to %d profiles. Max = %.4f (%s).",
      path, n_prof, max(rr_profile),
      dom_class[which.max(rr_profile)]
    ))
  }
  return(out)
}


# -- filter_profiles_to_rr_classes --------------------------------------------

#' Filter Profiles to Classes with Actual RR Estimates
#'
#' Drops profiles where every resistant class is absent from the RR table
#' (they would receive fallback_rr = 1 and silently distort the PAF), then
#' re-normalises the surviving profile probabilities to sum to 1.
#'
#' @param profiles_with_rr Named list from \code{assign_rr_to_profiles()}.
#' @param rr_table Data frame. Must have \code{pathogen_col} and
#'   \code{class_col}.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param probability_col Character. Default \code{"probability"}.
#' @param fallback_rr Numeric. Default \code{1}.
#'
#' @return Named list with unmatched profiles removed and probabilities
#'   re-normalised.
#' @export
filter_profiles_to_rr_classes <- function(
  profiles_with_rr,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  probability_col = "probability",
  fallback_rr = 1
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }
  if (!all(c(pathogen_col, class_col) %in% names(rr_table))) {
    stop(sprintf(
      "Columns '%s' and '%s' must be present in rr_table.",
      pathogen_col, class_col
    ))
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]
    rr_classes <- unique(rr_table[rr_table[[pathogen_col]] == path, class_col])
    rr_classes <- rr_classes[!is.na(rr_classes)]

    non_class_cols <- c(
      "profile", probability_col, "RR_LOS_profile",
      "dominant_class", "CI_lower_profile", "CI_upper_profile",
      "numerator", "PAF_LOS", "denominator"
    )
    class_cols <- setdiff(names(df), non_class_cols)
    matchable <- intersect(class_cols, rr_classes)

    if (length(matchable) == 0L) {
      warning(sprintf("'%s': no classes overlap with RR table -- all profiles dropped.", path))
      next
    }

    all_class_cols <- setdiff(names(df), non_class_cols)
    is_all_s <- if (length(all_class_cols) > 0L) {
      rowSums(df[, all_class_cols, drop = FALSE] == 1L, na.rm = TRUE) == 0L
    } else {
      rep(TRUE, nrow(df))
    }

    has_matchable_r <- if (length(matchable) == 1L) {
      df[[matchable]] == 1L
    } else {
      rowSums(df[, matchable, drop = FALSE] == 1L) > 0L
    }

    keep <- is_all_s | has_matchable_r
    n_before <- nrow(df)
    df_kept <- df[keep, , drop = FALSE]
    n_dropped <- n_before - nrow(df_kept)

    if (nrow(df_kept) == 0L) {
      warning(sprintf("'%s': all profiles dropped -- skipping.", path))
      next
    }

    prob_sum <- sum(df_kept[[probability_col]], na.rm = TRUE)
    if (prob_sum <= 0) {
      warning(sprintf("'%s': probability sum is zero -- skipping.", path))
      next
    }
    df_kept[[probability_col]] <- df_kept[[probability_col]] / prob_sum
    out[[path]] <- df_kept
    message(sprintf(
      "'%s': %d -> %d profiles kept (%d dropped). Prob re-normalised.",
      path, n_before, nrow(df_kept), n_dropped
    ))
  }
  return(out)
}


# -- compute_paf_rr_mortality -------------------------------------------------

#' Compute Mortality PAF per Resistance Profile (Levin formula)
#'
#' Computes PAF_kd_mortality for each pathogen x resistance profile using the
#' GBD multi-exposure Levin formula with mortality OR substituted for LOS RR:
#'
#' \deqn{\text{PAF}_{kd} =
#'   \frac{R'_{K\delta}(OR_{K\delta} - 1)}
#'        {1 + \sum_\delta R'_{K\delta}(OR_{K\delta} - 1)}}
#'
#' \strong{Pipeline:}
#' \preformatted{
#'   or_table          <- fit_mortality_rr_logistic(data, ...)
#'   p0_res            <- compute_p0(data, ...)
#'   rr_table          <- convert_or_to_rr(or_table, p0 = p0_res$p0)
#'   profiles_with_or  <- assign_rr_to_profiles(profiles_output,
#'                            rr_table = or_table, rr_col = "OR_death")
#'   paf_mort          <- compute_paf_rr_mortality(profiles_with_or)
#' }
#'
#' @param profiles_with_rr Named list from \code{assign_rr_to_profiles()}
#'   called with \code{rr_col = "OR_death"}.
#' @param probability_col Character. Default \code{"probability"}.
#' @param rr_profile_col Character. Default \code{"RR_LOS_profile"}.
#' @param profile_col Character. Default \code{"profile"}.
#' @param facility_col Character or \code{NULL}. Facility identifier column.
#'   Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. If provided, stored in
#'   the output for provenance tracking. Default \code{NULL}.
#'
#' @return Named list (one per pathogen) with \code{per_profile},
#'   \code{PAF_k_mort}, and \code{denominator_mort}.
#' @export
compute_paf_rr_mortality <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile",
  facility_col = NULL,
  facility_name = NULL
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]
    for (col in c(probability_col, rr_profile_col, profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf("Column '%s' not found in profiles for '%s'.", col, path))
      }
    }

    p <- df[[probability_col]]
    or <- df[[rr_profile_col]]
    numerator_vec <- p * (or - 1.0)
    denom <- 1.0 + sum(numerator_vec, na.rm = TRUE)

    if (!is.finite(denom) || denom <= 0) {
      warning(sprintf("'%s': PAF denominator = %.6g -- skipping.", path, denom))
      next
    }

    df$numerator_mort <- round(numerator_vec, 6L)
    df$PAF_mortality <- round(numerator_vec / denom, 6L)
    df$denominator_mort <- round(denom, 6L)

    paf_k_mort <- sum(df$PAF_mortality, na.rm = TRUE)
    out[[path]] <- list(
      per_profile      = df,
      PAF_k_mort       = round(paf_k_mort, 6L),
      denominator_mort = round(denom, 6L),
      facility_name    = if (!is.null(facility_name)) facility_name else NA_character_
    )
    message(sprintf(
      "'%s': PAF_k_mortality = %.4f | E[OR] = %.4f | %d profiles.",
      path, paf_k_mort, denom, nrow(df)
    ))
  }
  return(out)
}


# -- compute_fatal_resistance_prevalence ---------------------------------------

#' Compute Per-Profile Fatal Prevalence of Resistance R_Kdelta* (Eq. 13)
#'
#' For each **resistant** profile delta of pathogen K, computes the fatal
#' prevalence of resistance:
#'
#' \deqn{
#'   R^*_{K\delta} =
#'     \frac{R'_{K\delta} \cdot RR_{Kd^*}}
#'          {\left(1 - \textstyle\sum_\delta R'_{K\delta}\right)
#'           + \textstyle\sum_\delta R'_{K\delta} \cdot RR_{Kd^*}}
#' }
#'
#' where the sums \eqn{\sum_\delta} run over \strong{resistant profiles only}
#' (profiles with at least one resistant antibiotic class, i.e.
#' \code{dominant_class != "all_susceptible"}).  The term
#' \eqn{(1 - \sum_\delta R'_{K\delta})} is the susceptible fraction and must
#' be > 0 for the formula to give R*_Kdelta < 1.
#'
#' \strong{Important:} pass the output of \code{assign_rr_to_profiles()}
#' \emph{directly} -- do \emph{not} call \code{filter_profiles_to_rr_classes()}
#' first, because renormalisation would set \eqn{\sum R'_{K\delta} = 1} and
#' collapse the denominator, producing R*_Kdelta = 1 for every pathogen.
#'
#' @param profiles_with_rr Named list from
#'   \code{assign_rr_to_profiles(rr_col = "RR_death")}, one element per
#'   pathogen.  Each element must contain \code{probability_col},
#'   \code{rr_profile_col}, and \code{dominant_class_col}.
#' @param probability_col     Character. Profile prevalence column R'_Kdelta.
#'   Default \code{"probability"}.
#' @param rr_profile_col      Character. Dominant-class converted RR.
#'   Default \code{"RR_LOS_profile"}.
#' @param dominant_class_col  Character. Column identifying the dominant class
#'   (or \code{"all_susceptible"}).  Default \code{"dominant_class"}.
#' @param facility_col Character or \code{NULL}. Facility identifier column.
#'   Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. If provided, stored in
#'   the output for provenance tracking. Default \code{NULL}.
#'
#' @return Named list (one per pathogen).  Each element:
#' \describe{
#'   \item{\code{per_profile}}{Data frame of resistant profiles augmented with
#'     \code{R_star} (R*_Kdelta per profile) and \code{numerator_delta}.}
#'   \item{\code{R_K_star}}{Scalar: \eqn{\sum_\delta R^*_{K\delta}} =
#'     total fatal resistance prevalence for pathogen K.}
#'   \item{\code{sum_r_prime}}{Scalar: \eqn{\sum_\delta R'_{K\delta}}
#'     (resistant profiles only).}
#'   \item{\code{susceptible_fraction}}{Scalar:
#'     \eqn{1 - \sum_\delta R'_{K\delta}}.}
#'   \item{\code{denominator}}{Scalar: shared denominator for all profiles.}
#'   \item{\code{n_resistant_profiles}}{Integer: number of resistant profiles.}
#' }
#' @export
#'
#' @references
#' Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026. Eq. 13.
compute_fatal_resistance_prevalence <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  dominant_class_col = "dominant_class",
  facility_col = NULL,
  facility_name = NULL
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  rows <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]
    for (col in c(probability_col, rr_profile_col, dominant_class_col)) {
      if (!col %in% names(df)) {
        stop(sprintf("Column '%s' not found in profiles for '%s'.", col, path))
      }
    }

    # -- Keep only RESISTANT profiles (dominant_class != "all_susceptible") -
    df_res <- df[df[[dominant_class_col]] != "all_susceptible", , drop = FALSE]

    if (nrow(df_res) == 0L) {
      warning(sprintf("'%s': no resistant profiles found -- skipping.", path))
      next
    }

    r_prime <- df_res[[probability_col]] # R'_Kdelta
    rr_dominant <- df_res[[rr_profile_col]] # RR_Kd*(delta)

    # -- Numerator: Sigma_delta (R'_Kdelta x RR_Kd*(delta))  -- summed over all resistant delta -
    numerator <- sum(r_prime * rr_dominant, na.rm = TRUE)

    # -- Denominator: (1 - Sigma_delta R'_Kdelta) + Sigma_delta (R'_Kdelta x RR_Kd*(delta)) ----------
    sum_r_prime <- sum(r_prime, na.rm = TRUE)
    susceptible_frac <- 1.0 - sum_r_prime
    denominator <- susceptible_frac + numerator

    if (!is.finite(denominator) || denominator <= 0) {
      warning(sprintf("'%s': denominator = %.6g -- skipping.", path, denominator))
      next
    }

    # -- R*_K = numerator / denominator  (scalar per pathogen K) -----------
    R_K_star <- numerator / denominator

    message(sprintf(
      "'%s': Sigma(R'xRR)=%.4f | susc_frac=%.4f | denom=%.4f | R*_K=%.4f | %d resistant profiles.",
      path, numerator, susceptible_frac, denominator, R_K_star, nrow(df_res)
    ))

    rows[[path]] <- data.frame(
      pathogen             = path,
      facility_name        = if (!is.null(facility_name)) facility_name else NA_character_,
      sum_r_prime          = round(sum_r_prime, 6L),
      numerator            = round(numerator, 6L),
      susceptible_fraction = round(susceptible_frac, 6L),
      denominator          = round(denominator, 6L),
      R_K_star             = round(R_K_star, 6L),
      n_resistant_profiles = nrow(df_res),
      stringsAsFactors     = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out)
}


# -- compute_yll_associated ----------------------------------------------------

#' Compute Associated YLL per Pathogen (Eq. 11-12)
#'
#' Implements YLL -- Associated (Eq. 12):
#'
#' \deqn{
#'   YLL_K = \left[\sum_L \left(\sum_x D_L^x e_x^*\right) \times P_{KL}\right]
#'           \times R_K^*
#' }
#'
#' where \eqn{R_K^* = \sum_\delta R^*_{K\delta}} (sum of per-profile fatal
#' resistance prevalences from \code{compute_fatal_resistance_prevalence}).
#' For a single syndrome (BSI): YLL_K = base_YLL x P_KL x R_K*.
#'
#' Also returns per-profile YLL_Kdelta = base_YLL x P_KL x R*_Kdelta.
#'
#' @param base_yll Numeric scalar or list from
#'   \code{compute_base_yll_from_dl()} (uses \code{$total}).
#' @param P_Lk Named list from \code{calculate_P_Lk()} or the P_Lk data
#'   frame directly (column \code{pathogen_col} and \code{p_lk_col}).
#'   \code{NULL} applies equal weight 1/K.
#' @param fatal_resistance_prevalence Named list from
#'   \code{compute_fatal_resistance_prevalence()}.
#' @param pathogen_col Character. Pathogen column in P_Lk. Default
#'   \code{"pathogen"}.
#' @param p_lk_col Character. P_Lk value column. Default \code{"P_Lk"}.
#' @param profiles_with_rr Named list or \code{NULL}. Output from
#'   \code{assign_rr_to_profiles()} for per-profile / per-class breakdowns.
#' @param rr_profile_col Character. Dominant-class RR column in profile data.
#'   Default \code{"RR_death_profile"}.
#' @param probability_col Character. Profile prevalence column.
#'   Default \code{"probability"}.
#' @param dominant_class_col Character. Column identifying the dominant
#'   antibiotic class. Default \code{"dominant_class"}.
#' @param patient_data Data frame or \code{NULL}. Patient-level records for
#'   per-patient YLL output.
#' @param patient_col Character or \code{NULL}. Patient identifier column in
#'   \code{patient_data}.
#' @param outcome_col Character or \code{NULL}. Final outcome column in
#'   \code{patient_data}.
#' @param death_value Character. Value(s) indicating a fatal outcome.
#'   Default \code{"Death"}.
#' @param age_bin_col Character or \code{NULL}. Age bin column in
#'   \code{patient_data}.
#' @param sex_col Character or \code{NULL}. Sex column in \code{patient_data}.
#' @param syndrome_col Character or \code{NULL}. Syndrome column.
#' @param syndrome_name Character or \code{NULL}. If supplied, filters
#'   \code{patient_data} to this syndrome.
#' @param patient_pathogen_col Character or \code{NULL}. Pathogen column name
#'   in \code{patient_data} (e.g. \code{"organism_name"}).
#' @param le_path Character or \code{NULL}. Path to the life expectancy file.
#' @param male_value Character. Value in \code{sex_col} for males.
#'   Default \code{"Male"}.
#' @param female_value Character. Value in \code{sex_col} for females.
#'   Default \code{"Female"}.
#' @param age_bin_map Named character vector remapping non-standard age bin
#'   labels. Default \code{c("<1" = "0-1")}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{by_pathogen}}{Data frame: \code{pathogen}, \code{P_Lk},
#'     \code{sum_r_prime}, \code{susceptible_fraction}, \code{denominator},
#'     \code{R_K_star}, \code{base_yll}, \code{YLL_K}.}
#'   \item{\code{by_pathogen_profile}}{Data frame: one row per pathogen x
#'     resistant profile with \code{P_Lk}, \code{R_star},
#'     \code{YLL_Kdelta}.}
#'   \item{\code{total_yll}}{Scalar: sum of YLL_K.}
#'   \item{\code{base_yll_used}}{Scalar base YLL.}
#' }
#' @export
#'
#' @references
#' Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026. Eq. 12.
compute_yll_associated <- function(
  base_yll,
  P_Lk,
  fatal_resistance_prevalence,
  pathogen_col = "pathogen",
  p_lk_col = "P_Lk",
  # Optional: supply profiles_with_rr to get by_profile / by_class output
  profiles_with_rr = NULL,
  rr_profile_col = "RR_death_profile",
  probability_col = "probability",
  dominant_class_col = "dominant_class",
  # Optional: supply patient-level data to get by_patient output
  patient_data = NULL,
  patient_col = NULL,
  outcome_col = NULL,
  death_value = "Death",
  age_bin_col = NULL,
  sex_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  patient_pathogen_col = NULL, # pathogen column in patient_data (e.g. "organism_name")
  le_path = NULL,
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c("<1" = "0-1")
) {
  # -- 1. Extract scalar base YLL ---------------------------------------------
  if (is.list(base_yll) && "total" %in% names(base_yll)) {
    base_yll_val <- base_yll$total
  } else if (is.numeric(base_yll) && length(base_yll) == 1L) {
    base_yll_val <- base_yll
  } else {
    stop("'base_yll' must be a scalar or the list from compute_base_yll_from_dl().")
  }
  if (!is.finite(base_yll_val) || base_yll_val < 0) {
    stop(sprintf("base_yll_val must be finite and non-negative; got %.6g.", base_yll_val))
  }

  # -- 2. Validate fatal_resistance_prevalence --------------------------------
  if (!is.data.frame(fatal_resistance_prevalence)) {
    stop("'fatal_resistance_prevalence' must be the data frame from compute_fatal_resistance_prevalence().")
  }
  required_cols <- c("pathogen", "R_K_star")
  missing_cols <- setdiff(required_cols, names(fatal_resistance_prevalence))
  if (length(missing_cols)) {
    stop(sprintf(
      "'fatal_resistance_prevalence' missing columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  fr_df <- fatal_resistance_prevalence
  all_pathogens <- fr_df$pathogen

  # -- 3. Build P_Lk lookup ---------------------------------------------------
  if (!is.null(P_Lk)) {
    p_df <- if (!is.data.frame(P_Lk) && is.list(P_Lk) && "P_Lk" %in% names(P_Lk)) P_Lk$P_Lk else P_Lk
    if (!is.data.frame(p_df)) {
      stop("Cannot extract a data frame from 'P_Lk'. Provide the list from calculate_P_Lk().")
    }
    if (!pathogen_col %in% names(p_df)) {
      stop(sprintf("Column '%s' not found in P_Lk data frame.", pathogen_col))
    }
    if (!p_lk_col %in% names(p_df)) {
      stop(sprintf("Column '%s' not found in P_Lk data frame.", p_lk_col))
    }
    p_lookup <- stats::setNames(p_df[[p_lk_col]], p_df[[pathogen_col]])
  } else {
    k <- length(all_pathogens)
    p_lookup <- stats::setNames(rep(1.0 / k, k), all_pathogens)
    message(sprintf("P_Lk is NULL -- equal weighting (1 / %d = %.4f).", k, 1 / k))
  }

  # -- 4. Compute YLL_K = base_yll x P_KL x R*_K  (scalar per pathogen) -----
  p_lk_vals <- ifelse(
    all_pathogens %in% names(p_lookup),
    as.numeric(p_lookup[all_pathogens]),
    0.0
  )
  missing_p <- all_pathogens[!all_pathogens %in% names(p_lookup)]
  if (length(missing_p)) {
    warning(sprintf("Pathogens not in P_Lk (set to 0): %s", paste(missing_p, collapse = ", ")))
  }

  yll_k_vec <- base_yll_val * p_lk_vals * fr_df$R_K_star

  by_pathogen <- data.frame(
    pathogen             = all_pathogens,
    P_Lk                 = round(p_lk_vals, 6L),
    sum_r_prime          = round(fr_df$sum_r_prime, 6L),
    susceptible_fraction = round(fr_df$susceptible_fraction, 6L),
    numerator            = round(fr_df$numerator, 6L),
    denominator          = round(fr_df$denominator, 6L),
    R_K_star             = round(fr_df$R_K_star, 6L),
    base_yll             = round(base_yll_val, 4L),
    YLL_K                = round(yll_k_vec, 4L),
    stringsAsFactors     = FALSE
  )

  for (i in seq_len(nrow(by_pathogen))) {
    message(sprintf(
      "'%s': P_Lk=%.4f | R*_K=%.4f | YLL_K=%.4f",
      by_pathogen$pathogen[i], by_pathogen$P_Lk[i],
      by_pathogen$R_K_star[i], by_pathogen$YLL_K[i]
    ))
  }

  total_yll <- sum(by_pathogen$YLL_K, na.rm = TRUE)
  message(sprintf("Total associated YLL = %.4f years", total_yll))

  # -- 5. Per-patient YLL (optional) -----------------------------------------
  by_patient <- NULL

  if (!is.null(patient_data) && !is.null(patient_col) &&
    !is.null(outcome_col) && !is.null(age_bin_col) &&
    !is.null(sex_col) && !is.null(le_path) &&
    !is.null(pathogen_col)) {
    if (!file.exists(le_path)) {
      warning(sprintf("le_path not found ('%s'); skipping by_patient.", le_path))
    } else {
      # Load life table
      if (grepl("\\.xlsx$", le_path, ignore.case = TRUE)) {
        le_tbl <- readxl::read_excel(le_path)
      } else {
        le_tbl <- utils::read.csv(le_path, stringsAsFactors = FALSE)
      }

      le_tbl <- le_tbl %>%
        dplyr::select(age_bin_gbd, sex, life_expectancy) %>%
        dplyr::distinct()

      # Filter to deaths (and syndrome if supplied); recode age bins; normalise sex
      pt <- patient_data %>%
        dplyr::filter(.data[[outcome_col]] %in% death_value)
      if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
        pt <- pt %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
      }
      pt <- pt %>%
        dplyr::mutate(
          !!age_bin_col := dplyr::recode(
            as.character(.data[[age_bin_col]]), !!!age_bin_map
          ),
          .sex_norm = dplyr::case_when(
            .data[[sex_col]] == male_value ~ "Male",
            .data[[sex_col]] == female_value ~ "Female",
            TRUE ~ "Both"
          )
        )

      # Resolve the pathogen column name in patient_data
      pt_path_col <- if (!is.null(patient_pathogen_col)) patient_pathogen_col else pathogen_col

      if (!pt_path_col %in% names(pt)) {
        stop(sprintf(
          "patient_pathogen_col '%s' not found in patient_data. ",
          pt_path_col,
          "Set patient_pathogen_col to the pathogen column name in your patient data."
        ))
      }

      # Deduplicate to one row per patient x pathogen (drop antibiotic-level columns)
      dedup_cols <- intersect(
        c(patient_col, pt_path_col, age_bin_col, ".sex_norm"),
        names(pt)
      )
      pt <- pt %>%
        dplyr::select(dplyr::all_of(dedup_cols)) %>%
        dplyr::distinct()

      # Join life expectancy (sex-specific, fall back to "Both")
      pt <- pt %>%
        dplyr::left_join(
          le_tbl %>% dplyr::filter(sex != "Both") %>%
            dplyr::rename(life_expectancy_years = life_expectancy),
          by = c(stats::setNames("age_bin_gbd", age_bin_col),
            ".sex_norm" = "sex"
          )
        ) %>%
        dplyr::left_join(
          le_tbl %>% dplyr::filter(sex == "Both") %>%
            dplyr::select(age_bin_gbd, life_expectancy) %>%
            dplyr::rename(le_both = life_expectancy),
          by = stats::setNames("age_bin_gbd", age_bin_col)
        ) %>%
        dplyr::mutate(
          life_expectancy_years = dplyr::coalesce(
            life_expectancy_years, le_both
          )
        ) %>%
        dplyr::select(-le_both)

      # Join R*_K and P_Lk per pathogen using pt_path_col
      rk_lookup <- stats::setNames(fr_df$R_K_star, fr_df$pathogen)
      pt$R_K_star <- rk_lookup[pt[[pt_path_col]]]
      pt$P_Lk_val <- as.numeric(p_lookup[pt[[pt_path_col]]])

      pt$YLL_associated <- round(
        pt$life_expectancy_years * pt$P_Lk_val * pt$R_K_star, 4L
      )
      pt$life_expectancy_years <- round(pt$life_expectancy_years, 4L)

      by_patient <- pt

      message(sprintf(
        "by_patient: %d records | %d unique patients | mean life_expectancy = %.2f yrs",
        nrow(by_patient),
        dplyr::n_distinct(by_patient[[patient_col]]),
        mean(by_patient$life_expectancy_years, na.rm = TRUE)
      ))
    }
  }

  # -- 6. Per-profile and per-class YLL (optional) ---------------------------
  by_profile <- NULL
  by_class <- NULL

  if (!is.null(profiles_with_rr)) {
    profile_rows <- list()

    for (path in names(profiles_with_rr)) {
      df <- profiles_with_rr[[path]]

      for (col in c(probability_col, rr_profile_col, dominant_class_col)) {
        if (!col %in% names(df)) {
          stop(sprintf("Column '%s' not found in profiles for '%s'.", col, path))
        }
      }

      df_res <- df[df[[dominant_class_col]] != "all_susceptible", , drop = FALSE]
      if (nrow(df_res) == 0L) next

      r_prime <- df_res[[probability_col]]
      rr_dominant <- df_res[[rr_profile_col]]

      # Shared denominator (same as compute_fatal_resistance_prevalence)
      numerator_k <- sum(r_prime * rr_dominant, na.rm = TRUE)
      denom_k <- (1.0 - sum(r_prime, na.rm = TRUE)) + numerator_k
      if (!is.finite(denom_k) || denom_k <= 0) next

      # Per-profile R*_Kdelta and YLL_Kdelta
      R_star_delta <- (r_prime * rr_dominant) / denom_k

      p_lk_val <- if (path %in% names(p_lookup)) as.numeric(p_lookup[[path]]) else 0.0
      yll_kdelta <- base_yll_val * p_lk_val * R_star_delta

      profile_rows[[length(profile_rows) + 1L]] <- data.frame(
        pathogen = path,
        profile = df_res[[dominant_class_col]],
        R_prime = round(r_prime, 6L),
        RR_dominant = round(rr_dominant, 6L),
        R_star_delta = round(R_star_delta, 6L),
        YLL_Kdelta = round(yll_kdelta, 4L),
        stringsAsFactors = FALSE
      )
    }

    if (length(profile_rows) > 0L) {
      by_profile <- do.call(rbind, profile_rows)
      rownames(by_profile) <- NULL

      # Aggregate to dominant class level (matches yll_attr$by_profile structure)
      by_class <- by_profile %>%
        dplyr::group_by(pathogen, profile) %>%
        dplyr::summarise(
          YLL_Kdelta = sum(YLL_Kdelta, na.rm = TRUE),
          .groups = "drop"
        )
    }
  }

  list(
    by_pathogen   = by_pathogen,
    by_profile    = by_profile,
    by_class      = by_class,
    by_patient    = by_patient,
    total_yll     = round(total_yll, 4L),
    base_yll_used = round(base_yll_val, 4L)
  )
}


# -- compute_yll_attributable --------------------------------------------------

#' Compute Attributable YLL per Pathogen x Resistance Profile (Eq. 14-15)
#'
#' Implements YLL -- Attributable (Bhaswati Ganguli, DALY Methodology, 2026):
#'
#' \strong{Eq. 15 -- Mortality PAF per profile delta:}
#' \deqn{
#'   \text{MortPAF}_{K\delta} =
#'     \frac{R'_{K\delta} \cdot (RR_{Kd^*}(\delta) - 1)}
#'          {1 + \sum_{\delta'} R'_{K\delta'} \cdot (RR_{Kd^*}(\delta') - 1)}
#' }
#'
#' Numerator is \strong{per profile} delta; denominator is \strong{shared} across
#' all resistant profiles of pathogen K. The sum \eqn{\sum_{\delta'}} runs
#' over resistant profiles only (\code{dominant_class != "all_susceptible"}).
#'
#' \strong{Eq. 14 -- Attributable YLL per profile:}
#' \deqn{
#'   YLL_{K\delta} = \text{BaseYLL} \times P_{KL} \times \text{MortPAF}_{K\delta}
#' }
#'
#' \strong{Attributable YLL per pathogen:}
#' \deqn{YLL_K = \sum_\delta YLL_{K\delta}}
#'
#' @param profiles_with_rr Named list from
#'   \code{assign_rr_to_profiles(rr_col = "RR_death")}, one element per
#'   pathogen. Each element must have \code{probability_col},
#'   \code{rr_profile_col}, and \code{dominant_class_col}.
#' @param base_yll Numeric scalar or list from
#'   \code{compute_base_yll_from_dl()} (uses \code{$total}).
#' @param P_Lk Data frame (pathogen column + P_Lk column) or list from
#'   \code{calculate_P_Lk()}. \code{NULL} applies equal weighting 1/K.
#' @param pathogen_col Character. Pathogen column in \code{P_Lk}.
#'   Default \code{"pathogen"}.
#' @param p_lk_col Character. P_Lk value column. Default \code{"P_Lk"}.
#' @param probability_col Character. Profile prevalence R'_Kdelta column.
#'   Default \code{"probability"}.
#' @param rr_profile_col Character. Dominant-class mortality RR column.
#'   Default \code{"RR_death_profile"}.
#' @param dominant_class_col Character. Dominant class column.
#'   Default \code{"dominant_class"}.
#' @param patient_data Data frame or \code{NULL}. Patient-level records for
#'   per-patient attributable YLL output.
#' @param patient_col Character or \code{NULL}. Patient identifier column in
#'   \code{patient_data}.
#' @param outcome_col Character or \code{NULL}. Final outcome column in
#'   \code{patient_data}.
#' @param death_value Character. Value(s) indicating a fatal outcome.
#'   Default \code{"Death"}.
#' @param age_bin_col Character or \code{NULL}. Age bin column in
#'   \code{patient_data}.
#' @param sex_col Character or \code{NULL}. Sex column in \code{patient_data}.
#' @param syndrome_col Character or \code{NULL}. Syndrome column.
#' @param syndrome_name Character or \code{NULL}. If supplied, filters
#'   \code{patient_data} to this syndrome.
#' @param patient_pathogen_col Character or \code{NULL}. Pathogen column name
#'   in \code{patient_data} (e.g. \code{"organism_name"}).
#' @param le_path Character or \code{NULL}. Path to the life expectancy file.
#' @param male_value Character. Value in \code{sex_col} for males.
#'   Default \code{"Male"}.
#' @param female_value Character. Value in \code{sex_col} for females.
#'   Default \code{"Female"}.
#' @param age_bin_map Named character vector remapping non-standard age bin
#'   labels. Default \code{c("<1" = "0-1")}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{by_profile}}{Data frame: one row per pathogen x resistant
#'     profile with \code{R_prime}, \code{RR_dominant}, \code{excess_risk},
#'     \code{numerator_delta}, \code{denominator}, \code{MortPAF},
#'     \code{P_Lk}, \code{base_yll}, \code{YLL_Kdelta}.}
#'   \item{\code{by_pathogen}}{Data frame: one row per pathogen with
#'     \code{P_Lk}, \code{denominator}, \code{sum_MortPAF}, \code{YLL_K}.}
#'   \item{\code{total_yll}}{Scalar: \eqn{\sum_K YLL_K}.}
#'   \item{\code{base_yll_used}}{Scalar base YLL used.}
#' }
#' @export
#'
#' @references
#' Bhaswati Ganguli. DALY Methodology for AMR. March 2026. Eq. 14-15.
compute_yll_attributable <- function(
  profiles_with_rr,
  base_yll,
  P_Lk,
  pathogen_col = "pathogen",
  p_lk_col = "P_Lk",
  probability_col = "probability",
  rr_profile_col = "RR_death_profile",
  dominant_class_col = "dominant_class",
  # Optional: supply patient-level data to get by_patient output
  patient_data = NULL,
  patient_col = NULL,
  outcome_col = NULL,
  death_value = "Death",
  age_bin_col = NULL,
  sex_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  patient_pathogen_col = NULL,
  le_path = NULL,
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c("<1" = "0-1")
) {
  # -- 1. Extract scalar base YLL ---------------------------------------------
  if (is.list(base_yll) && "total" %in% names(base_yll)) {
    base_yll_val <- base_yll$total
  } else if (is.numeric(base_yll) && length(base_yll) == 1L) {
    base_yll_val <- base_yll
  } else {
    stop("'base_yll' must be a scalar or the list from compute_base_yll_from_dl().")
  }
  if (!is.finite(base_yll_val) || base_yll_val < 0) {
    stop(sprintf("base_yll_val must be finite and non-negative; got %.6g.", base_yll_val))
  }

  # -- 2. Validate profiles_with_rr ------------------------------------------
  if (!is.list(profiles_with_rr)) {
    stop("'profiles_with_rr' must be the list from assign_rr_to_profiles().")
  }

  # -- 3. Build P_Lk lookup ---------------------------------------------------
  all_pathogens <- names(profiles_with_rr)

  if (!is.null(P_Lk)) {
    p_df <- if (!is.data.frame(P_Lk) && is.list(P_Lk) && "P_Lk" %in% names(P_Lk)) P_Lk$P_Lk else P_Lk
    if (!is.data.frame(p_df)) {
      stop("Cannot extract a data frame from 'P_Lk'.")
    }
    if (!pathogen_col %in% names(p_df)) {
      stop(sprintf("Column '%s' not found in P_Lk.", pathogen_col))
    }
    if (!p_lk_col %in% names(p_df)) {
      stop(sprintf("Column '%s' not found in P_Lk.", p_lk_col))
    }
    p_lookup <- stats::setNames(p_df[[p_lk_col]], p_df[[pathogen_col]])
  } else {
    k <- length(all_pathogens)
    p_lookup <- stats::setNames(rep(1.0 / k, k), all_pathogens)
    message(sprintf("P_Lk is NULL -- equal weighting (1 / %d = %.4f).", k, 1 / k))
  }

  # -- 4. Compute MortPAF and YLL per pathogen x profile ---------------------
  profile_rows <- list()
  pathogen_rows <- list()

  for (path in all_pathogens) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col, dominant_class_col)) {
      if (!col %in% names(df)) {
        stop(sprintf("Column '%s' not found for '%s'.", col, path))
      }
    }

    # Resistant profiles only
    df_res <- df[df[[dominant_class_col]] != "all_susceptible", , drop = FALSE]

    if (nrow(df_res) == 0L) {
      warning(sprintf("'%s': no resistant profiles -- skipping.", path))
      next
    }

    r_prime <- df_res[[probability_col]] # R'_Kdelta  (per profile)
    rr_dominant <- df_res[[rr_profile_col]] # RR_Kd*(delta) (per profile)
    excess <- rr_dominant - 1.0 # (RR_Kd* - 1) per profile

    # -- Eq. 15 denominator (shared for all delta of pathogen K) ---------------
    # 1 + Sigma_delta R'_Kdelta x (RR_Kd*(delta) - 1)
    denominator <- 1.0 + sum(r_prime * excess, na.rm = TRUE)

    if (!is.finite(denominator) || denominator <= 0) {
      warning(sprintf("'%s': denominator = %.6g -- skipping.", path, denominator))
      next
    }

    # -- Eq. 15 numerator (per profile delta) ----------------------------------
    # R'_Kdelta x (RR_Kd*(delta) - 1)
    numerator_delta <- r_prime * excess

    # -- MortPAF_Kdelta = numerator_delta / denominator ------------------------
    mort_paf_delta <- numerator_delta / denominator

    # -- P_KL --------------------------------------------------------------
    p_lk_val <- if (path %in% names(p_lookup)) {
      as.numeric(p_lookup[[path]])
    } else {
      warning(sprintf("'%s' not in P_Lk -- set to 0.", path))
      0.0
    }

    # -- Eq. 14: YLL_Kdelta = BaseYLL x P_KL x MortPAF_Kdelta --------------------
    yll_kdelta <- base_yll_val * p_lk_val * mort_paf_delta

    # -- YLL_K = Sigma_delta YLL_Kdelta ------------------------------------------------
    yll_k <- sum(yll_kdelta, na.rm = TRUE)

    message(sprintf(
      "'%s': P_Lk=%.4f | denom=%.4f | SigmaMortPAF=%.4f | YLL_K_attr=%.4f",
      path, p_lk_val, denominator,
      sum(mort_paf_delta, na.rm = TRUE), yll_k
    ))

    profile_rows[[path]] <- data.frame(
      pathogen = path,
      profile = df_res[[dominant_class_col]],
      R_prime = round(r_prime, 6L),
      RR_dominant = round(rr_dominant, 6L),
      excess_risk = round(excess, 6L),
      numerator_delta = round(numerator_delta, 6L),
      denominator = round(denominator, 6L),
      MortPAF = round(mort_paf_delta, 6L),
      P_Lk = round(p_lk_val, 6L),
      base_yll = round(base_yll_val, 4L),
      YLL_Kdelta = round(yll_kdelta, 4L),
      stringsAsFactors = FALSE
    )

    pathogen_rows[[path]] <- data.frame(
      pathogen = path,
      P_Lk = round(p_lk_val, 6L),
      denominator = round(denominator, 6L),
      sum_MortPAF = round(sum(mort_paf_delta, na.rm = TRUE), 6L),
      base_yll = round(base_yll_val, 4L),
      YLL_K = round(yll_k, 4L),
      stringsAsFactors = FALSE
    )
  }

  by_profile <- do.call(rbind, profile_rows)
  rownames(by_profile) <- NULL
  by_pathogen <- do.call(rbind, pathogen_rows)
  rownames(by_pathogen) <- NULL

  total_yll <- sum(by_pathogen$YLL_K, na.rm = TRUE)
  message(sprintf("Total attributable YLL = %.4f years", total_yll))

  # -- Per-patient attributable YLL (optional) --------------------------------
  by_patient <- NULL

  if (!is.null(patient_data) && !is.null(patient_col) &&
    !is.null(outcome_col) && !is.null(age_bin_col) &&
    !is.null(sex_col) && !is.null(le_path)) {
    if (!file.exists(le_path)) {
      warning(sprintf("le_path not found ('%s'); skipping by_patient.", le_path))
    } else {
      # Load life table
      if (grepl("\\.xlsx$", le_path, ignore.case = TRUE)) {
        le_tbl <- readxl::read_excel(le_path)
      } else {
        le_tbl <- utils::read.csv(le_path, stringsAsFactors = FALSE)
      }

      le_tbl <- le_tbl %>%
        dplyr::select(age_bin_gbd, sex, life_expectancy) %>%
        dplyr::distinct()

      # Filter to deaths + optional syndrome; recode age bins; normalise sex
      pt <- patient_data %>%
        dplyr::filter(.data[[outcome_col]] %in% death_value)
      if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
        pt <- pt %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
      }
      pt <- pt %>%
        dplyr::mutate(
          !!age_bin_col := dplyr::recode(
            as.character(.data[[age_bin_col]]), !!!age_bin_map
          ),
          .sex_norm = dplyr::case_when(
            .data[[sex_col]] == male_value ~ "Male",
            .data[[sex_col]] == female_value ~ "Female",
            TRUE ~ "Both"
          )
        )

      pt_path_col <- if (!is.null(patient_pathogen_col)) patient_pathogen_col else pathogen_col

      if (!pt_path_col %in% names(pt)) {
        stop(sprintf(
          "patient_pathogen_col '%s' not found in patient_data.", pt_path_col
        ))
      }

      # Select and deduplicate to one row per patient x pathogen
      dedup_cols <- intersect(
        c(patient_col, pt_path_col, age_bin_col, ".sex_norm"),
        names(pt)
      )
      pt <- pt %>%
        dplyr::select(dplyr::all_of(dedup_cols)) %>%
        dplyr::distinct()

      # Join life expectancy (sex-specific, fall back to "Both")
      pt <- pt %>%
        dplyr::left_join(
          le_tbl %>% dplyr::filter(sex != "Both") %>%
            dplyr::rename(life_expectancy_years = life_expectancy),
          by = c(stats::setNames("age_bin_gbd", age_bin_col),
            ".sex_norm" = "sex"
          )
        ) %>%
        dplyr::left_join(
          le_tbl %>% dplyr::filter(sex == "Both") %>%
            dplyr::select(age_bin_gbd, life_expectancy) %>%
            dplyr::rename(le_both = life_expectancy),
          by = stats::setNames("age_bin_gbd", age_bin_col)
        ) %>%
        dplyr::mutate(
          life_expectancy_years = dplyr::coalesce(life_expectancy_years, le_both)
        ) %>%
        dplyr::select(-le_both)

      # Join sum_MortPAF and P_Lk per pathogen
      paf_lookup <- stats::setNames(by_pathogen$sum_MortPAF, by_pathogen$pathogen)
      pt$sum_MortPAF <- paf_lookup[pt[[pt_path_col]]]
      pt$P_Lk_val <- as.numeric(p_lookup[pt[[pt_path_col]]])

      pt$YLL_attributable <- round(
        pt$life_expectancy_years * pt$P_Lk_val * pt$sum_MortPAF, 4L
      )
      pt$life_expectancy_years <- round(pt$life_expectancy_years, 4L)

      by_patient <- pt

      message(sprintf(
        "by_patient (attributable): %d records | %d unique patients | mean YLL_attr = %.2f yrs",
        nrow(by_patient),
        dplyr::n_distinct(by_patient[[patient_col]]),
        mean(by_patient$YLL_attributable, na.rm = TRUE)
      ))
    }
  }

  list(
    by_profile    = by_profile,
    by_pathogen   = by_pathogen,
    by_patient    = by_patient,
    total_yll     = round(total_yll, 4L),
    base_yll_used = round(base_yll_val, 4L)
  )
}
