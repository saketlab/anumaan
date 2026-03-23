# Identify top N pathogens by occurrence

Ranks pathogens by the number of records (rows) in the dataset,
optionally filtered to a specific syndrome, specimen, facility, or
outcome. Returns the top N pathogens overall or broken down per
facility.

## Usage

``` r
get_top_pathogens(
  data,
  pathogen_col,
  n = 5L,
  syndrome_col = NULL,
  syndrome_name = NULL,
  specimen_col = NULL,
  specimen_name = NULL,
  outcome_col = NULL,
  outcome_name = NULL,
  facility_col = NULL,
  facility_name = NULL
)
```

## Arguments

- data:

  Data frame of facility-level records.

- pathogen_col:

  Character. Pathogen column.

- n:

  Integer. Number of top pathogens to return. Default 5.

- syndrome_col:

  Character or NULL. Filter to this syndrome column.

- syndrome_name:

  Character or NULL. Syndrome value to filter to.

- specimen_col:

  Character or NULL. Filter to this specimen column.

- specimen_name:

  Character or NULL. Specimen value to filter to.

- outcome_col:

  Character or NULL. Filter to this outcome column.

- outcome_name:

  Character or NULL. Outcome value to filter to (e.g., `"Died"` or
  `"Discharged"`).

- facility_col:

  Character or NULL. Facility identifier column. When provided without
  `facility_name`, returns top N per facility.

- facility_name:

  Character or NULL. If provided, restricts to that facility only before
  ranking.

## Value

Data frame with columns: `pathogen_col`, `n_records`, `rank` (1 = most
common), and `facility_col` if supplied.

## Details

Use this to decide which pathogens to focus on before calling
`calculate_P_Lk_prime_BSI()`,
[`calculate_cfr_lk()`](https://saketlab.github.io/anumaan/reference/calculate_cfr_lk.md),
or
[`calculate_YLD()`](https://saketlab.github.io/anumaan/reference/calculate_YLD.md).
