# Calculate Fatal Pathogen Distribution (P\_{LK})

Computes the fatal pathogen distribution for a given infectious syndrome
(L) using facility-level microbiology data. This quantity represents the
fractional contribution of each pathogen (K) to **fatal** infection
cases, and is used in Step 4 of YLL estimation (Eq. 11).

## Usage

``` r
calculate_P_Lk(
  data,
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
  pathogen_name = NULL
)
```

## Arguments

- data:

  Data frame of facility-level microbiology records.

- syndrome_col:

  Character. Column containing infectious syndrome labels (L).

- syndrome_name:

  Character. Syndrome to analyse (e.g. `"Bloodstream infection"`).

- polymicrobial_col:

  Character. Column flagging polymicrobial patients (1 = polymicrobial,
  0 = monomicrobial).

- patient_col:

  Character. Unique patient identifier column.

- pathogen_col:

  Character. Pathogen (organism) column (K).

- outcome_col:

  Character. Final patient outcome column.

- death_value:

  Character. Value indicating a fatal outcome. Default `"Death"`.

- specimen_col:

  Character or `NULL`. Specimen type column. When supplied,
  `specimen_name` must also be provided.

- specimen_name:

  Character or `NULL`. Specimen value to restrict to (e.g. `"Blood"`).

- glass_ref:

  Character vector of valid pathogen names, or a data frame with columns
  `specimen` and `pathogen`. When supplied, polymicrobial patients are
  restricted to pathogens on the GLASS list for the given specimen.
  Monomicrobial patients are never filtered. `NULL` skips GLASS
  filtering.

- facility_col:

  Character or `NULL`. Facility identifier column. When supplied without
  `facility_name`, returns both facility-level and pooled P_LK.

- facility_name:

  Character or `NULL`. Restricts computation to a single named facility.

- pathogen_name:

  Character vector or `NULL`. Restricts output to specific pathogen(s).

## Value

A named list:

- `P_Lk`:

  Data frame of pooled P_LK across facilities. Columns: `pathogen_col`,
  `N_F_LK` (weighted fatal patient count for K), `N_F_L` (total fatal
  patients for L), `P_Lk` (fraction).

- `facility_level`:

  Per-facility P_LK data frame (only when `facility_col` is supplied and
  `facility_name` is `NULL`).

## Details

**Unit of analysis:** the patient. Each fatal patient contributes a
total weight of 1, split equally across their valid pathogens:
\$\$w\_{r,k} = \frac{1}{m_r}\$\$ where \\m_r\\ is the number of distinct
pathogens for patient \\r\\. Monomicrobial patients (\\m_r = 1\\)
receive weight 1; polymicrobial patients have their death split equally
(Approach 1).

**Pooled formula across facilities:** \$\$P\_{LK}^{\text{pooled}} =
\frac{\sum_f N^{F}\_{f,L,K}}{\sum_f N^{F}\_{f,L}}\$\$

## References

Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026.

## Examples

``` r
if (FALSE) { # \dontrun{
pkL <- calculate_P_Lk(
  data              = bsi_data,
  syndrome_col      = "infectious_syndrome",
  syndrome_name     = "Bloodstream infection",
  polymicrobial_col = "is_polymicrobial",
  patient_col       = "PatientInformation_id",
  pathogen_col      = "organism_name",
  outcome_col       = "final_outcome",
  death_value       = "Death",
  facility_col      = "center_name"
)
pkL$P_Lk
pkL$facility_level
} # }
```
