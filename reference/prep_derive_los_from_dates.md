# Derive Length of Stay from Date Columns

Computes LOS when missing (or always when `overwrite = TRUE`). Adds
`los_method` and `los_confidence` audit columns. For AIIMS data, uses
`unit_admission_date` + `unit_duration_days` as a fallback when primary
date columns are absent.

## Usage

``` r
prep_derive_los_from_dates(
  data,
  admission_col = "admission_date",
  outcome_col = "outcome_date",
  los_col = "los_days",
  unit_admission_col = "unit_admission_date",
  unit_duration_col = "unit_duration_days",
  overwrite = FALSE
)
```

## Arguments

- data:

  Data frame.

- admission_col:

  Character. Admission date column. Default "admission_date".

- outcome_col:

  Character. Outcome/discharge date column. Default "outcome_date".

- los_col:

  Character. Target LOS column to populate. Default "los_days".

- unit_admission_col:

  Character. AIIMS unit admission date fallback. Default
  "unit_admission_date".

- unit_duration_col:

  Character. AIIMS unit duration fallback. Default "unit_duration_days".

- overwrite:

  Logical. Recalculate even when `los_col` is already present. Default
  FALSE.

## Value

Data frame with `los_col`, `los_method`, and `los_confidence` populated
where possible.

## Details

Priority order:

1.  `outcome_col` - `admission_col`

2.  `unit_admission_date` + `unit_duration_days` (AIIMS fallback)

3.  Already-present `los_col` value retained as-is
