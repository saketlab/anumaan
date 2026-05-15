# Assert Standard Names Are Present

Confirms that a column mapping succeeded for critical standard columns.
Stops on failure when strict = TRUE.

## Usage

``` r
prep_assert_standard_names(data, required_standard_names, strict = TRUE)
```

## Arguments

- data:

  Data frame after column mapping.

- required_standard_names:

  Character vector of standard column names that must be present.

- strict:

  Logical. Stop on failure. Default TRUE.

## Value

Invisibly returns data. Stops or warns on missing columns.
