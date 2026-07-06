# Compute Pairwise Co-resistance Using Pearson Back-calculation

For each antibiotic-class pair (A, B) and each pathogen \[x stratum\],
this function:

1.  Computes the Pearson correlation \\\rho\_{AB}\\ between binary
    resistance indicators across co-tested isolates.

2.  Back-calculates pairwise prevalence using the GBD formula: \$\$P(A
    \cap B) = P(A) P(B) + \rho\_{AB} \sqrt{P(A)(1-P(A)) \cdot
    P(B)(1-P(B))}\$\$ where \\P(A)\\ and \\P(B)\\ come from `marginals`
    (which may include externally overridden values).

3.  Caps the result to \\\min(P(A), P(B))\\ and floors it to 0. Cells
    below the floor use the independence product as fallback.

## Usage

``` r
compute_pairwise_from_data(
  data_wide,
  marginals,
  col_map,
  panel_map,
  stratify_by = NULL,
  outcome_col = NULL,
  min_co_tested = 10L
)
```

## Arguments

- data_wide:

  Tibble. Output of `preprocess_for_profiles()$data_wide`.

- marginals:

  Tibble. Output of
  [`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md).

- col_map:

  Named list. `col_map_resolved`.

- panel_map:

  Named list. Same panel as used in preprocessing.

- stratify_by:

  Character vector or `NULL`. Must match what was used in
  [`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md).
  Default `NULL`.

- outcome_col:

  Character or `NULL`. Default `NULL`.

- min_co_tested:

  Integer. Pairs with fewer co-tested isolates are reported with
  `method = "independence_fallback"`. Default `10L`.

## Value

Tibble with one row per pathogen x class-pair \[x stratum\]:

- `pathogen`, `antibiotic_class_1`, `antibiotic_class_2`:

  Keys.

- `n_co_tested`:

  Isolates tested for both classes.

- `rho`:

  Pearson correlation between binary resistance indicators.

- `pairwise_prevalence`:

  Back-calculated P(AnB), capped and floored.

- `method`:

  `"pearson_back_calc"` or `"independence_fallback"`.
