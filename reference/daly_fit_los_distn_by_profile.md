# Estimate Profile-Specific Median LOS (Section 5.3)

For each resistance profile \\\delta\\ of each pathogen \\k\\, fits
parametric distributions (default: Weibull, log-normal, gamma) to the
LOS of the patients with that profile, \\\\LOS_i : \delta_i =
\delta\\\\, selects the distribution with the lowest AIC, and takes the
**theoretical median** of the selected distribution: \$\$LOS\_{k\delta}
= Median\_{f\_\delta}(LOS)\$\$ The same procedure applied to the fully
susceptible reference profile \\\delta_S\\ gives \\LOS\_{k\delta_S}\\,
and `RR_LOS = LOS_{k delta} / LOS_{k delta_S}`.

## Usage

``` r
daly_fit_los_distn_by_profile(
  data,
  profile_data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  profile_id_col = patient_id_col,
  profile_col = "profile_delta",
  profile_pathogen_col = NULL,
  ref_profile = NULL,
  syndrome_name = NULL,
  organism_name = NULL,
  facility_name = NULL,
  distributions = c("weibull", "lnorm", "gamma"),
  max_los = 365,
  min_n = 10L
)
```

## Arguments

- data:

  Data frame used for LOS (patient, facility, organism, dates, final
  outcome).

- profile_data:

  Data frame with one profile per patient (or per patient x pathogen):
  columns `profile_id_col`, `profile_col`, and optionally
  `profile_pathogen_col`.

- patient_id_col:

  Character. Patient identifier in `data`.

- facility_col:

  Character. Facility / centre column in `data`.

- organism_col:

  Character. Organism column in `data`; used to assign pathogens when
  `profile_pathogen_col` is `NULL`.

- syndrome_col:

  Character. Syndrome column (used only with `syndrome_name`).

- date_admission_col:

  Character. Admission date column.

- date_discharge_col:

  Character. Final outcome / discharge date column.

- final_outcome_col:

  Character. Final outcome column.

- final_outcome_value:

  Character. Value indicating discharge.

- profile_id_col:

  Character. ID column in `profile_data` matching `patient_id_col`.
  Default `patient_id_col`.

- profile_col:

  Character. Profile label column in `profile_data` (e.g. `"RRSS"`).
  Default `"profile_delta"`.

- profile_pathogen_col:

  Character or `NULL`. Pathogen column in `profile_data`, when profiles
  are pathogen-specific. If `NULL`, pathogens are taken from
  `organism_col` in `data`.

- ref_profile:

  Character or `NULL`. Label of the fully susceptible reference profile.
  If `NULL`, the label consisting only of `"S"` characters (e.g.
  `"SSSS"`) is used per pathogen.

- syndrome_name:

  Character or `NULL`. Restrict to this syndrome.

- organism_name:

  Character vector or `NULL`. Restrict to these pathogens.

- facility_name:

  Character or `NULL`. Restrict to this facility.

- distributions:

  Character vector of candidate distributions. Default
  `c("weibull", "lnorm", "gamma")`.

- max_los:

  Numeric. Maximum LOS in days. Default `365`.

- min_n:

  Integer. Minimum LOS observations required to fit a profile. Default
  `10L`.

## Value

Data frame with one row per pathogen x profile: `pathogen`,
`profile_delta`, `is_reference`, `n_patients`, `n_los` (LOS
observations), `median_LOS` (days), `median_LOS_years`, `best_dist`,
`aic`, `median_LOS_ref`, `RR_LOS`, `status` (`"fitted"`,
`"insufficient_n"`, `"fit_failed"`) and `syndrome_scope`. Profiles that
cannot be fitted are kept with `NA` estimates.

## Details

**Profiles are not derived here.** Each patient's resistance profile
must be supplied in `profile_data` (e.g. the `profile_delta` labels from
the `profile_*` pipeline). Patients without a profile, or with
conflicting profiles, are excluded and reported in a message.

LOS is computed by
[`daly_compute_patient_los()`](https://saketlab.github.io/anumaan/reference/daly_compute_patient_los.md)
(final outcome date - admission date, discharged patients only, 0 \< LOS
\<= `max_los`). Theoretical medians: Weibull \\\lambda (\ln 2)^{1/k}\\;
log-normal \\e^{\mu}\\; gamma `qgamma(0.5, shape, rate)`.

## Two options for profile-level LOS

The package offers two ways to obtain LOS for each resistance profile:

1.  **Profile-specific (this function; Section 5.3).** Fits the LOS of
    the patients who actually have each profile \\\delta\\ and compares
    it with the fully susceptible profile \\\delta_S\\. Requires
    patient-level profiles (e.g. from the `profile_*` pipeline) and
    enough patients per profile (`min_n`).

2.  **Per-class with max rule.**
    [`daly_fit_los_rr_distribution()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr_distribution.md)
    fits resistant vs susceptible LOS separately for each antibiotic
    class, then
    [`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
    gives each profile the estimate of its dominant (highest-RR)
    resistant class. Needs no patient-level profiles and pools data
    across profiles, so it is a fallback or sensitivity analysis when
    profiles are sparse; its reference group is "susceptible to that
    class", not fully susceptible.

Note: the YLD functions
([`daly_calc_yld_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_associated.md),
[`daly_calc_yld_attributable()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_attributable.md))
currently read the `LOSR_years`/`LOSS_years` columns produced by option
2; the output of this function is not yet joined into that format.

## See also

[`daly_fit_los_rr_distribution`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr_distribution.md)
and
[`daly_assign_rr_to_profiles`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
(option 2)
