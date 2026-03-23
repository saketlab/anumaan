# Enrich Age

Derives Age when missing using DOB and culture date. Uses multi-path
logic with confidence scoring.

## Usage

``` r
enrich_age(
  data,
  age_col = "Age",
  dob_col = "DOB",
  date_col = "date_of_culture",
  overwrite = FALSE
)
```

## Arguments

- data:

  Data frame

- age_col:

  Character. Age column. Default "Age".

- dob_col:

  Character. Date of birth column. Default "DOB".

- date_col:

  Character. Reference date for age calculation. Default
  "date_of_culture".

- overwrite:

  Logical. If TRUE, recalculates age even if present. Default FALSE
  (only fill missing).

## Value

Data frame with Age column enriched

## Examples

``` r
if (FALSE) { # \dontrun{
data_enriched <- enrich_age(data)
} # }
```
