# Detect Preprocessing Capabilities

Inspects a data frame (after column mapping) and returns a named logical
vector indicating which preprocessing steps can be executed given the
columns that are present.

## Usage

``` r
detect_preprocessing_capabilities(data)
```

## Arguments

- data:

  Data frame after column mapping.

## Value

Named logical vector. TRUE means the step can run; FALSE means required
inputs are absent and the step will be skipped with a warning.

## Details

This is used by
[`run_preprocess()`](https://saketlab.github.io/anumaan/reference/run_preprocess.md)
to decide which layers to run and which to skip gracefully.
