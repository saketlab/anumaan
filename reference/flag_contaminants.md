# Flag Contaminant Organisms

Identifies likely contaminant organisms using multi-path logic.

## Usage

``` r
flag_contaminants(
  data,
  method = "auto",
  organism_col = "organism_normalized",
  specimen_col = "specimen_type"
)
```

## Arguments

- data:

  Data frame

- method:

  Character. "auto" (try all methods), "device_based", "heuristic",
  "provided". Default "auto".

- organism_col:

  Character. Normalized organism column.

- specimen_col:

  Character. Specimen type column.

## Value

Data frame with is_contaminant, contaminant_confidence,
contaminant_method columns
