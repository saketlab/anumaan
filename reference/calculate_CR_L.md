# Calculate the CFR adjustment factor (CR_L)

Computes CR_L, the factor that adjusts a hospital-derived CFR to account
for infection cases managed outside the inpatient setting. The
adjustment type for each syndrome is looked up from the `adjustment_ref`
table (loaded from `inst/extdata/adjustment_for_CFR`).

## Usage

``` r
calculate_CR_L(
  data,
  syndrome_col,
  syndrome_name,
  patient_col,
  visit_type_col,
  inpatient_value = "Inpatient",
  outpatient_value = "Outpatient",
  adjustment_ref,
  facility_col = NULL,
  facility_name = NULL
)
```

## Arguments

- data:

  Data frame of facility-level records.

- syndrome_col:

  Character. Column containing infectious syndrome labels.

- syndrome_name:

  Character. Syndrome to compute CR_L for.

- patient_col:

  Character. Unique patient identifier column.

- visit_type_col:

  Character. Column indicating visit type per record (inpatient /
  outpatient).

- inpatient_value:

  Character. Value in `visit_type_col` that denotes an inpatient visit.
  Default `"Inpatient"`.

- outpatient_value:

  Character. Value in `visit_type_col` that denotes an outpatient visit.
  Default `"Outpatient"`.

- adjustment_ref:

  Data frame with columns `infectious_syndrome` and
  `adjustment_factor_on_CFR`. Load from
  `inst/extdata/adjustment_for_CFR`.

- facility_col:

  Character or NULL. Facility identifier column. When provided (and
  `facility_name` is NULL), CR_L is returned per facility.

- facility_name:

  Character or NULL. If provided, restricts to that facility only.

## Value

Data frame with columns `syndrome` (= `syndrome_name`),
`adjustment_type`, `CR_L`, and (when `facility_col` is supplied)
`facility_col`. Additional columns (`n_inpatient`, `n_total`, etc.) give
the raw counts used.

## Details

Three adjustment types are supported:

- None:

  CR_L = 1. The hospital CFR applies directly (e.g., BSI, Meningitis,
  hospital-acquired infections).

- Inpatient ratio:

  CR_L = (patients with \\\ge\\ 1 inpatient visit) / (all patients).
  Used when community cases are captured partly in outpatient data
  (e.g., community-acquired LRI, UTI).

- Outpatient to inpatient ratio:

  CR_L = (patients with \\\ge\\ 1 outpatient AND \\\ge\\ 1 inpatient
  visit) / (patients with \\\ge\\ 1 outpatient visit). Used for
  syndromes where OP-to-IP transition captures disease severity (e.g.,
  STI, Skin, Eye, Oral, Bone/joint infections).
