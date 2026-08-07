# Check pairwise co-testing overlap per hospital x class-pair

Reports the full 2x2 co-testing breakdown (n_RR/n_RS/n_SR/n_SS), not
just the total co-tested count. A pair can clear `min_pairwise_cotested`
while showing no variation at all (e.g. 100 co-tested events that are
all resistant-resistant) – such a pair carries no information about
residual cross-class correlation, so `sufficient` additionally requires
all four cells to have at least `min_pairwise_cell` observations.

## Usage

``` r
.check_pairwise_cotesting(
  event_class_data,
  class_cols,
  upper_re_col,
  min_pairwise_cotested = 30L,
  min_pairwise_cell = 1L
)
```
