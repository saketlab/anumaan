# Calculate case fatality ratio by syndrome and pathogen (CFR\_{Lk})

Computes the case fatality ratio (CFR) for each pathogen (k) within a
specified infectious syndrome (L) using facility-level microbiology
data.

## Usage

``` r
calculate_cfr_lk(
  data,
  syndrome_col,
  syndrome_name,
  specimen_col,
  specimen_name,
  polymicrobial_col,
  patient_col,
  pathogen_col,
  outcome_col,
  death_value = "Died",
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

  Character. Column containing specimen type.

- specimen_name:

  Character. Specimen to restrict to (e.g., `"Blood"`).

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

  Character. Value indicating death. Default `"Died"`.

- glass_ref:

  Character vector of valid pathogen names, or a data frame with columns
  `specimen` and `pathogen`. Applied to polymicrobial patients only.
  `NULL` skips GLASS filtering.

- facility_col:

  Character or NULL. Facility identifier column.

- facility_name:

  Character or NULL. Restricts to a single facility.

- pathogen_name:

  Character vector or NULL. Filter to specific pathogen(s).

## Value

A list:

- cfr_table:

  Pooled CFR: `pathogen_col`, `weighted_deaths`, `weighted_total`,
  `CFR_LK`.

- facility_level:

  Per-facility CFR (only when `facility_col` is supplied and
  `facility_name` is NULL).

## Details

The unit of analysis is the **patient**. Each patient contributes total
weight 1, distributed equally across their valid pathogens via a
fractional weight \\1/m_r\\ (where \\m_r\\ is the number of distinct
valid pathogens for patient r). This ensures polymicrobial patients are
not double-counted.

For **polymicrobial patients** (`polymicrobial_col == 1`), only
pathogens listed in the GLASS reference (`glass_ref`) for the given
specimen type are retained before weighting. Monomicrobial patients
(`polymicrobial_col == 0`) are never filtered.

When `facility_col` is supplied and `facility_name` is NULL, both
per-facility and pooled CFR are returned. The pooled CFR sums weighted
deaths and totals across facilities before dividing.
