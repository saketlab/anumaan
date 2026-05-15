# Log Data Source Provenance

Records metadata about the source of data entering the pipeline. Returns
a provenance metadata list attached as an attribute.

## Usage

``` r
prep_log_source(
  data,
  study_type = "generic",
  centre_name = NULL,
  file_path = NULL,
  sheet_name = NULL
)
```

## Arguments

- data:

  Data frame.

- study_type:

  Character. One of "stewardship", "surveillance", "aiims_icu_bsi",
  "generic".

- centre_name:

  Character. Name of the contributing centre.

- file_path:

  Character. Path to the source file (optional).

- sheet_name:

  Character. Sheet name within Excel file (optional).

## Value

Data frame with `.provenance` attribute attached.
