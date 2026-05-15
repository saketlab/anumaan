# Standardize Column Names to Package Convention

Convenience wrapper around
[`prep_build_column_map()`](https://saketlab.github.io/anumaan/reference/prep_build_column_map.md) +
[`prep_apply_column_map()`](https://saketlab.github.io/anumaan/reference/prep_apply_column_map.md).
Builds an exact-match map from `mapping` aliases, optionally extends it
with fuzzy matching, then applies the rename in one call.

## Usage

``` r
prep_standardize_column_names(
  data,
  mapping = default_column_mappings,
  fuzzy_match = TRUE,
  fuzzy_threshold = 0.3,
  interactive = FALSE
)
```

## Arguments

- data:

  A data frame with raw column names.

- mapping:

  Named list of standard_name -\> aliases vectors. Default uses
  `default_column_mappings`.

- fuzzy_match:

  Logical. Enable fuzzy matching for unresolved standards. Default TRUE.

- fuzzy_threshold:

  Numeric. Max Jaro-Winkler distance (0-1). Default 0.3.

- interactive:

  Logical. Prompt user to confirm fuzzy matches. Default FALSE.

## Value

List: `data` (renamed), `mapping_log`, `unmapped`.

## Details

For more control (e.g. dataset-specific maps, post-rename assertions)
use the staged API directly:
[`prep_build_column_map()`](https://saketlab.github.io/anumaan/reference/prep_build_column_map.md)
-\>
[`prep_apply_column_map()`](https://saketlab.github.io/anumaan/reference/prep_apply_column_map.md)
-\>
[`prep_assert_standard_names()`](https://saketlab.github.io/anumaan/reference/prep_assert_standard_names.md).
