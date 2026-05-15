# Compute Polymicrobial Weights

Calculates proportional weights for polymicrobial infections.
Monomicrobial patients always receive weight = 1.0. For polymicrobial
patients, weight is computed per organism within each episode using one
of three methods.

## Usage

``` r
prep_compute_poly_weights(
  data,
  episode_col = "episode_id",
  organism_col = "organism_normalized",
  polymicrobial_col = "is_polymicrobial",
  method = "monomicrobial_proportion",
  weight_map = NULL,
  facility_col = NULL,
  facility_name = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)

compute_polymicrobial_weight(
  data,
  episode_col = "episode_id",
  organism_col = "organism_normalized",
  polymicrobial_col = "is_polymicrobial",
  method = "monomicrobial_proportion",
  weight_map = NULL,
  facility_col = NULL,
  facility_name = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame with `episode_id`, `is_polymicrobial` (0/1), and organism
  columns (output of
  [`prep_flag_polymicrobial()`](https://saketlab.github.io/anumaan/reference/prep_flag_polymicrobial.md)).

- episode_col:

  Character. Episode ID column. Default "episode_id".

- organism_col:

  Character. Organism column. Default "organism_normalized".

- polymicrobial_col:

  Character. Polymicrobial flag column (0/1). Default
  "is_polymicrobial".

- method:

  Character. One of "monomicrobial_proportion" (default), "equal", or
  "manual".

- weight_map:

  Named numeric vector. Custom organism weights when
  `method = "manual"`.

- facility_col:

  Character or NULL. Scope monomicrobial reference pool per facility.

- facility_name:

  Character or NULL. Filter to this facility.

- syndrome_col:

  Character or NULL. Scope reference pool per syndrome.

- syndrome_name:

  Character or NULL. Filter to this syndrome.

## Value

Data frame with `poly_weight`, `weight_method`, and `weight_confidence`
columns.
