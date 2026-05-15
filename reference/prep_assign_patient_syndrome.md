# Assign One Syndrome Per Patient or Event

Selects the highest-priority syndrome per patient using the infectious
syndrome hierarchy. Matching is case-insensitive (uses normalized
names).

## Usage

``` r
prep_assign_patient_syndrome(
  data,
  patient_col = "patient_id",
  syndrome_col = "syndrome",
  score_col = "icd_score",
  hierarchy_ref = NULL,
  hierarchy_syndrome_col = "infectious_syndrome",
  hierarchy_rank_col = "rank",
  keep_all_candidates = TRUE,
  collapse_unspecified_respiratory = FALSE,
  keep_only_burden_syndromes = FALSE
)
```

## Arguments

- data:

  Data frame with one row per syndrome candidate.

- patient_col:

  Character. Patient or event ID column. Default `"patient_id"`.

- syndrome_col:

  Character. Syndrome column. Default `"syndrome"`.

- score_col:

  Character or NULL. ICD match score column for secondary priority.
  Default `"icd_score"`.

- hierarchy_ref:

  Data frame or NULL. Syndrome hierarchy (from
  `load_syndrome_hierarchy()`). If NULL, loads automatically.

- hierarchy_syndrome_col:

  Character. Syndrome column in hierarchy. Default
  `"infectious_syndrome"`.

- hierarchy_rank_col:

  Character. Rank column in hierarchy. Default `"rank"`.

- keep_all_candidates:

  Logical. TRUE -\> return all rows with `syndrome_selected` flag; FALSE
  -\> return one row per `patient_col` (selected only). Default TRUE.

- collapse_unspecified_respiratory:

  Logical. Map "Other unspecified respiratory site infections" -\>
  "Lower respiratory infections" before selection. Default FALSE.

- keep_only_burden_syndromes:

  Logical. Retain only syndromes where
  `contributes_to_amr_burden == TRUE`. Default FALSE.

## Value

Data frame. When `keep_all_candidates = FALSE`, includes
`contributes_to_amr_burden` column from the hierarchy.

## Details

Priority order per patient:

1.  Hierarchy rank (lower = higher priority)

2.  ICD match score (higher = higher priority), when `score_col` present

3.  Alphabetical tie-break

This function absorbs the logic of the former
`infer_patient_syndrome_long()`: use
`collapse_unspecified_respiratory = TRUE` and
`keep_only_burden_syndromes = TRUE` to replicate that behaviour.
