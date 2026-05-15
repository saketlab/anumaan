# Derive HAI/CAI Infection Type

Infers Community-Acquired (CAI) vs Hospital-Acquired (HAI) infection
using the admission-to-culture date gap (default cutoff: 2 days / 48
hours).

## Usage

``` r
prep_derive_hai_cai(
  data,
  infection_type_col = "infection_type",
  admission_col = "date_of_admission",
  culture_col = "date_of_culture",
  hai_cutoff = 2,
  overwrite = FALSE
)
```

## Arguments

- data:

  Data frame.

- infection_type_col:

  Character. Infection type column. Default "infection_type".

- admission_col:

  Character. Admission date column. Default "date_of_admission".

- culture_col:

  Character. Culture date column. Default "date_of_culture".

- hai_cutoff:

  Numeric. Days after admission to classify as HAI. Default 2.

- overwrite:

  Logical. Recalculate even when already present. Default FALSE.

## Value

Data frame with `infection_type` enriched.
