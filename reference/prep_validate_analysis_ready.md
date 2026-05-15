# Validate Analysis-Ready Dataset

Final gate check before burden estimation. Confirms that the minimum
standard AMR schema fields are present and have acceptable completeness.

## Usage

``` r
prep_validate_analysis_ready(
  data,
  min_rows = 10L,
  max_missing_pct = 30,
  stop_on_failure = FALSE
)
```

## Arguments

- data:

  Analysis-ready data frame.

- min_rows:

  Integer. Minimum acceptable row count. Default 10.

- max_missing_pct:

  Numeric. Maximum percent missing allowed for core fields. Default 30.

- stop_on_failure:

  Logical. Stop if validation fails. Default FALSE.

## Value

Invisibly returns a list: passes (logical), issues (character vector).
