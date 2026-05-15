# Split Polymicrobial Episodes by Strategy

Applies a chosen strategy for handling polymicrobial episodes:

- fractional:

  Keep all rows; `poly_weight` is applied downstream during DALY
  calculation.

- exclude:

  Drop polymicrobial episodes from the main analysis frame.

- separate:

  Retain polymicrobial episodes in a separate data frame stored in
  `attr(result, "poly_data")` for sensitivity analyses.

## Usage

``` r
prep_split_poly_episode(
  data,
  strategy = c("fractional", "exclude", "separate"),
  polymicrobial_col = "is_polymicrobial"
)
```

## Arguments

- data:

  Data frame with `is_polymicrobial` column (output of
  [`prep_flag_polymicrobial()`](https://saketlab.github.io/anumaan/reference/prep_flag_polymicrobial.md)).

- strategy:

  Character. "fractional" (default), "exclude", or "separate".

- polymicrobial_col:

  Character. Polymicrobial flag column (0/1). Default
  "is_polymicrobial".

## Value

Data frame with polymicrobial rows handled per strategy. When
`strategy = "separate"`, the poly subset is in
`attr(result, "poly_data")`.
