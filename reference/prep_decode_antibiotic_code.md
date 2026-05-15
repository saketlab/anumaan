# Decode Antibiotic Short Codes to Full Names

Converts antibiotic short codes (e.g. AMK, AMP, AMPSUL) to full
antibiotic names using a caller-supplied code map CSV.

## Usage

``` r
prep_decode_antibiotic_code(
  data,
  code_col = "antibiotic_name_raw",
  map_path = NULL,
  output_col = "antibiotic_name_std"
)
```

## Arguments

- data:

  Data frame with a column of antibiotic codes.

- code_col:

  Character. Column containing antibiotic codes. Default
  "antibiotic_name_raw".

- map_path:

  Character. Path to the code map CSV file. Must have columns `code` and
  `antibiotic_name_full`. If NULL, falls back to
  `inst/extdata/antibiotic_code_map.csv` if it exists.

- output_col:

  Character. Name of the output column for decoded names. Default
  "antibiotic_name_std".

## Value

Data frame with `output_col` column added.

## Details

The code map must have columns `code` and `antibiotic_name_full`. This
function is dataset-agnostic: the code map and the column containing the
codes are supplied by the caller (typically from the analysis layer).
