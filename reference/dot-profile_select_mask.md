# (Internal) Select a Reproducible, Eligibility-Preserving Mask

Stratified masking over a `Y` matrix (any `D`), by group x class x
observed R/S value; for every group x class cell, at most enough cells
are masked to keep tested/resistant/susceptible counts at or above the
`min_*_after_mask` floors; only observed (non-`NA`) cells are
candidates; no event is ever left with zero observed classes.

## Usage

``` r
.profile_select_mask(
  Y,
  ids,
  fraction_to_mask,
  max_mask_per_cell,
  min_tested_after_mask,
  min_resistant_after_mask,
  min_susceptible_after_mask,
  seed
)
```
