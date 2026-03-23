# Generate Preprocessing Report

Creates a comprehensive report documenting all transformations applied
during preprocessing, including column mappings, data quality metrics,
and transformation summaries.

## Usage

``` r
generate_preprocessing_report(
  raw_data,
  processed_data,
  preprocessing_log,
  config,
  output_file = NULL,
  include_plots = TRUE
)
```

## Arguments

- raw_data:

  Original input data frame (before preprocessing)

- processed_data:

  Final output data frame (after preprocessing)

- preprocessing_log:

  List containing logs from each preprocessing step

- config:

  AMR configuration object used

- output_file:

  Optional. Path to save report. If NULL, returns list. Supports .html,
  .pdf, .txt, .rds formats.

- include_plots:

  Logical. Include data quality plots. Default TRUE.

## Value

If output_file is NULL, returns a list with report components.
Otherwise, saves report to file and returns file path.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate report
report <- generate_preprocessing_report(
  raw_data = my_raw_data,
  processed_data = result$data,
  preprocessing_log = result$log,
  config = result$config,
  output_file = "preprocessing_report.html"
)
} # }
```
