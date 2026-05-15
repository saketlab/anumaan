# Parse Age Bin Labels to Breaks

Parses age-bin labels into numeric breakpoints for interval binning.

## Usage

``` r
parse_age_bin_labels(labels)
```

## Arguments

- labels:

  Character vector of age-bin labels, such as `"<1"`, `"1-5"`, or
  `"85+"`.

## Value

A list with `breaks` (numeric breakpoints) and `labels` (original
labels).
