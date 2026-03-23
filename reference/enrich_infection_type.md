# Enrich Infection Type

Infers Community-Acquired (CAI) vs Hospital-Acquired (HAI) infection
using admission-culture date gap.

## Usage

``` r
enrich_infection_type(
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

  Data frame

- infection_type_col:

  Character. Infection type column. Default "infection_type".

- admission_col:

  Character. Admission date. Default "date_of_admission".

- culture_col:

  Character. Culture date. Default "date_of_culture".

- hai_cutoff:

  Numeric. Days after admission to classify as HAI. Default 2 (48
  hours).

- overwrite:

  Logical. Recalculate even if present. Default FALSE.

## Value

Data frame with infection_type enriched

## Examples

``` r
if (FALSE) { # \dontrun{
# Default 2-day cutoff
data_enriched <- enrich_infection_type(data)

# 3-day cutoff
data_enriched <- enrich_infection_type(data, hai_cutoff = 3)
} # }
```
