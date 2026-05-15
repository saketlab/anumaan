# Build and Validate a Column Map Against a Dataset

Validates a user-supplied column map against the actual columns present
in a data frame. Reports which mappings are present, which are absent,
and which data columns are not covered.

## Usage

``` r
prep_build_column_map(data, column_map = NULL, custom_map = NULL)
```

## Arguments

- data:

  Data frame with raw column names.

- column_map:

  Named character vector: standard_name = raw_name. Pass a
  dataset-specific map from your analysis layer. If NULL, generic alias
  detection is attempted.

- custom_map:

  Named character vector. Additional raw -\> standard overrides applied
  on top of `column_map`. Optional.

## Value

Named character vector with only mappings whose raw names exist in data.

## Details

Dataset-specific maps (ICMR stewardship, surveillance, AIIMS, etc.)
should be defined in the analysis layer and passed in via `column_map`.
If no map is supplied, falls back to auto-detection using
`default_column_mappings` aliases.
