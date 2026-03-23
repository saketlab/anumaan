# Compute Polymicrobial Weights

Calculates proportional weights for polymicrobial infections.
Monomicrobial patients always receive weight = 1.0. For polymicrobial
patients, weight is computed per organism within each episode using one
of three methods.

## Usage

``` r
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
  [`flag_polymicrobial()`](https://saketlab.github.io/anumaan/reference/flag_polymicrobial.md)).

- episode_col:

  Character. Episode ID column. Default `"episode_id"`.

- organism_col:

  Character. Organism column. Default `"organism_normalized"`.

- polymicrobial_col:

  Character. Polymicrobial flag column (0/1). Default
  `"is_polymicrobial"`.

- method:

  Character. Weighting method: `"monomicrobial_proportion"` (default),
  `"equal"`, or `"manual"`.

- weight_map:

  Named numeric vector. Custom organism weights when
  `method = "manual"`.

- facility_col:

  Character or `NULL`. Facility column. When supplied, monomicrobial
  proportions are computed per-facility so each facility's local
  organism distribution is used as reference.

- facility_name:

  Character or `NULL`. When supplied together with `facility_col`, data
  are first filtered to that facility.

- syndrome_col:

  Character or `NULL`. Syndrome column. When supplied, monomicrobial
  proportions are computed per-syndrome.

- syndrome_name:

  Character or `NULL`. When supplied together with `syndrome_col`, data
  are first filtered to that syndrome.

## Value

Data frame with `polymicrobial_weight` column (range 0-1), plus
`weight_method` and `weight_confidence` audit columns. `episode_id` is
removed (internal use only).

## Details

When `facility_col` or `syndrome_col` are supplied, the monomicrobial
reference pool (method `"monomicrobial_proportion"`) is computed within
each stratum, so the reference distribution is local to each facility /
syndrome rather than global.

## Examples

``` r
if (FALSE) { # \dontrun{
# Global monomicrobial proportion weights
data_weighted <- compute_polymicrobial_weight(data_flagged)

# Per-facility reference distribution
data_weighted <- compute_polymicrobial_weight(
  data_flagged,
  facility_col = "center_name"
)

# Filter to one facility + one syndrome, then weight
data_weighted <- compute_polymicrobial_weight(
  data_flagged,
  facility_col  = "center_name",
  facility_name = "PGIMER",
  syndrome_col  = "infectious_syndrome",
  syndrome_name = "Bloodstream infection"
)
} # }
```
