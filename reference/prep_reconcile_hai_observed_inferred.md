# Reconcile Observed and Inferred HAI/CAI Classification

Where both an observed (raw) infection type and an inferred
(date-derived) type exist, flags rows where they disagree. Retains
observed classification by default and marks discordant records.

## Usage

``` r
prep_reconcile_hai_observed_inferred(
  data,
  observed_col = "infection_type_observed",
  inferred_col = "infection_type_inferred",
  output_col = "infection_type"
)
```

## Arguments

- data:

  Data frame.

- observed_col:

  Character. Observed infection type column. Default
  "infection_type_observed".

- inferred_col:

  Character. Inferred infection type column. Default
  "infection_type_inferred".

- output_col:

  Character. Reconciled output column. Default "infection_type".

## Value

Data frame with reconciled `infection_type` and `hai_discordant` logical
flag.
