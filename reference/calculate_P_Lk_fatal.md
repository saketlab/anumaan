# Calculate fatal pathogen distribution (P\_{Lk})

Computes the fatal pathogen distribution for a given infectious syndrome
(L) using facility-level microbiology data. This quantity represents the
fractional contribution of each pathogen (k) to **fatal** infection
cases, and is used in YLL / mortality burden estimation.

## Usage

``` r
calculate_P_Lk_fatal(
  data,
  syndrome_col,
  syndrome_name,
  specimen_col = NULL,
  specimen_name = NULL,
  polymicrobial_col,
  patient_col,
  pathogen_col,
  outcome_col,
  death_value = "Death",
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

  Character. Syndrome to analyse.

- specimen_col:

  Character or NULL. Column containing specimen type. If `NULL`, no
  specimen filter is applied.

- specimen_name:

  Character or NULL. Specimen to restrict to (e.g., `"Blood"`). Required
  when `specimen_col` is provided.

- polymicrobial_col:

  Character. Column flagging polymicrobial patients (1 = polymicrobial,
  0 = monomicrobial).

- patient_col:

  Character. Unique patient identifier column.

- pathogen_col:

  Character. Pathogen (organism) column (k).

- outcome_col:

  Character. Final patient outcome column.

- death_value:

  Character. Value indicating a fatal outcome. Default `"Death"`.

- glass_ref:

  Character vector of valid pathogen names, or a data frame with columns
  `specimen` and `pathogen`. Applied to polymicrobial patients only.
  When `specimen_name` is `NULL` and `glass_ref` is a data frame, all
  pathogens in the reference are used regardless of specimen. `NULL`
  skips GLASS filtering.

- facility_col:

  Character or NULL. Facility identifier column. When provided without
  `facility_name`, returns both facility-level and pooled P_LK.

- facility_name:

  Character or NULL. Restricts to a single facility.

- pathogen_name:

  Character vector or NULL. Filter to specific pathogen(s).

## Value

A list:

- P_Lk_fatal:

  Pooled P_LK data frame: `pathogen_col`, `N_F_LK`, `N_F_L`,
  `P_Lk_fatal`.

- facility_level:

  Per-facility P_LK (only when `facility_col` is supplied and
  `facility_name` is NULL).

## Details

The unit of analysis is the **patient**. Each patient contributes total
weight 1, distributed equally across their valid pathogens. For a
patient with \\m_r\\ valid pathogens, each pathogen receives weight
\\1/m_r\\.

For **polymicrobial patients** (`polymicrobial_col == 1`), only
pathogens listed in the GLASS reference (`glass_ref`) for the given
specimen type are retained before weighting. Monomicrobial patients
(`polymicrobial_col == 0`) are never filtered.

The pooled formula across facilities is: \$\$P\_{LK}^{\text{pooled}} =
\frac{\sum_f N^{F}\_{f,L,K}}{\sum_f N^{F}\_{f,L}}\$\$
