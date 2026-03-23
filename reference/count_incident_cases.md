# Count incident cases by syndrome from facility data

Counts the number of unique patients per infectious syndrome directly
from facility-level data. This is the direct-count approach to incidence
– no CFR, no CR_L adjustment, no pathogen weighting. Use this when you
want raw facility-reported case counts rather than the formula-derived
estimate from
[`calculate_incidence_L()`](https://saketlab.github.io/anumaan/reference/calculate_incidence_L.md).

## Usage

``` r
count_incident_cases(
  data,
  syndrome_col,
  syndrome_name,
  patient_col,
  facility_col = NULL,
  facility_name = NULL,
  pathogen_col = NULL,
  pathogen_name = NULL
)
```

## Arguments

- data:

  Data frame of facility-level records.

- syndrome_col:

  Character. Column containing infectious syndrome labels.

- syndrome_name:

  Character. Syndrome to count cases for (e.g.,
  `"Bloodstream infections"`).

- patient_col:

  Character. Unique patient identifier column. Cases are counted as
  distinct patients, not rows.

- facility_col:

  Character or NULL. Facility identifier column. When provided without
  `facility_name`, counts are broken down per facility. Default `NULL`.

- facility_name:

  Character or NULL. If provided, restricts the count to that facility
  only. Default `NULL`.

- pathogen_col:

  Character or NULL. Pathogen identifier column. Required when
  `pathogen_name` is specified. Default `NULL`.

- pathogen_name:

  Character or NULL. If provided, restricts the count to the specified
  pathogen(s). Default `NULL`.

## Value

Data frame with columns: `syndrome_col`, `n_cases` (unique patient
count), and `facility_col` if supplied.

## Details

Results are returned facility-wise when `facility_col` is supplied and
`facility_name` is NULL. When `facility_name` is specified, only that
facility is returned. When no facility information is provided, a single
pooled count across all records is returned.
