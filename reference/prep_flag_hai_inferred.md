# Flag HAI Inferred vs Observed

Adds an `infection_type_src` column indicating whether the infection
type was present in the raw data ("observed"), derived from dates
("inferred"), or could not be determined ("unknown").

## Usage

``` r
prep_flag_hai_inferred(
  data,
  infection_type_col = "infection_type",
  infection_type_method_col = "infection_type_method"
)
```

## Arguments

- data:

  Data frame.

- infection_type_col:

  Character. Raw/observed infection type column. Default
  "infection_type".

- infection_type_method_col:

  Character. Method column added by
  [`prep_derive_hai_cai()`](https://saketlab.github.io/anumaan/reference/prep_derive_hai_cai.md).
  Default "infection_type_method".

## Value

Data frame with `infection_type_src` column added.
