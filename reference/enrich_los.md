# Enrich Length of Stay

Derives Length of Stay (LOS) from admission and outcome dates.

## Usage

``` r
enrich_los(
  data,
  los_col = "Length_of_stay",
  admission_col = "date_of_admission",
  outcome_col = "date_of_final_outcome",
  overwrite = FALSE
)
```

## Arguments

- data:

  Data frame

- los_col:

  Character. LOS column name. Default "Length_of_stay".

- admission_col:

  Character. Admission date. Default "date_of_admission".

- outcome_col:

  Character. Outcome date. Default "date_of_final_outcome".

- overwrite:

  Logical. Recalculate even if present. Default FALSE.

## Value

Data frame with LOS enriched

## Examples

``` r
if (FALSE) { # \dontrun{
data_enriched <- enrich_los(data)
} # }
```
