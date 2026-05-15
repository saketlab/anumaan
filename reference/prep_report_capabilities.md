# Report Preprocessing Capabilities

Prints a human-readable summary of which preprocessing steps can be run
given the current data columns.

## Usage

``` r
prep_report_capabilities(data)
```

## Arguments

- data:

  Data frame after column mapping. Or pass the result of
  [`detect_preprocessing_capabilities()`](https://saketlab.github.io/anumaan/reference/detect_preprocessing_capabilities.md)
  directly.

## Value

Invisibly returns the capabilities vector.
