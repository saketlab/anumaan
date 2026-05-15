# Classify Readmission Type

Applies readmission classification rules to a data frame that already
has a `readmission_class` column, re-standardizing the values.

## Usage

``` r
prep_classify_readmission(data, readmission_col = "readmission_class")
```

## Arguments

- data:

  Data frame.

- readmission_col:

  Character. Readmission column to standardize. Default
  "readmission_class".

## Value

Data frame with standardized readmission classification.

## Details

Useful when data arrives with non-standard labels that need mapping to
the controlled vocabulary used by
[`prep_flag_readmission()`](https://saketlab.github.io/anumaan/reference/prep_flag_readmission.md).
