# Compute Pairwise Co-resistance Matrices per Pathogen

For every pathogen and every pair of antibiotic classes \\(c_i, c_j)\\
that were both tested: \$\$ T\_{k,i,j} = \sum_e \mathbf{1}(c_i \text{
tested} \land c_j \text{ tested}) \$\$ \$\$ R\_{k,i,j} = \sum_e
\mathbf{1}(R\_{e,k,c_i}=1 \land R\_{e,k,c_j}=1) \$\$ \$\$
\text{Prev}\_{k,i,j} = R\_{k,i,j} \\/\\ T\_{k,i,j} \$\$

## Usage

``` r
compute_pairwise_coresistance(
  marginal_output,
  pathogen_col = "organism_name",
  org_group_col = "org_group",
  isolate_col = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  min_co_tested = 10,
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL
)
```

## Arguments

- marginal_output:

  The list returned by
  [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md).
  The `$class_long` element is used as input.

- pathogen_col:

  Character. Must match the column name used in Step 1. Default
  `"organism_name"`.

- org_group_col:

  Character. Default `"org_group"`.

- isolate_col:

  Character. Default `"isolate_id"`.

- antibiotic_class_col:

  Character. Default `"antibiotic_class"`.

- min_co_tested:

  Integer. Pairwise cells with fewer than `min_co_tested` co-tested
  isolates are set to `NA` in the prevalence matrix. Default `10`.

- facility_col:

  Character or `NULL`. Column identifying the facility/site. When
  provided together with `facility_name`, filters `class_long` to the
  specified facility before building matrices. For this to work,
  [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md)
  must have been called with the same `facility_col` argument so that
  the column is present in `$class_long`. Both must be supplied together
  or both left `NULL`. Default `NULL`.

- facility_name:

  Character or `NULL`. Facility value to retain. Default `NULL`.

- outcome_col:

  Character or `NULL`. Column containing patient outcomes. When provided
  together with `outcome_value`, filters `class_long` to isolates with
  the specified outcome before building matrices.
  [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md)
  must have been called with the same `outcome_col` so that the column
  is present in `$class_long`. Both must be supplied together or both
  left `NULL`. Default `NULL`.

- outcome_value:

  Character or `NULL`. Outcome value to retain (e.g. `"discharged"`,
  `"dead"`). Default `NULL`.

## Value

Named list. One entry per pathogen (keyed by pathogen name), each
containing:

- `prevalence`:

  Square symmetric matrix of pairwise co-resistance rates. Diagonal and
  cells with `n < min_co_tested` are `NA`.

- `T_matrix`:

  Integer matrix of co-tested isolate counts.

- `R_matrix`:

  Integer matrix of co-resistant isolate counts.

- `classes`:

  Character vector of class names (row/column order).

## Details

Matrices are computed for **all** tested classes – no filtering by
marginal resistance or GBD core list is applied here.

## Examples

``` r
set.seed(1)
ast_class <- data.frame(
  isolate_id       = rep(sprintf("ISO%02d", 1:80), each = 3),
  organism_name    = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 120),
  org_group        = "Enterobacterales",
  antibiotic_class = rep(c("Carbapenems", "3GC", "Fluoroquinolones"), times = 80),
  antibiotic_value = sample(c("S", "R"), 240, replace = TRUE)
)
marg <- compute_marginal_resistance(ast_class)
#> Collapsed to class level: 80 isolates | 2 pathogens | 3 classes.
#> No classes with marginal_resistance <= 0.
#> min_n_tested = 30: all 6 pathogen-class combination(s) passed the threshold.
#> class_long after min_n_tested filter: 240 isolate-class row(s) retained.
marg$marginal
#> # A tibble: 6 × 6
#>   organism_name         org_group        antibiotic_class n_tested n_resistant
#>   <chr>                 <chr>            <chr>               <int>       <int>
#> 1 Escherichia coli      Enterobacterales 3GC                    40          22
#> 2 Escherichia coli      Enterobacterales Fluoroquinolones       40          20
#> 3 Escherichia coli      Enterobacterales Carbapenems            40          16
#> 4 Klebsiella pneumoniae Enterobacterales Fluoroquinolones       40          24
#> 5 Klebsiella pneumoniae Enterobacterales 3GC                    40          19
#> 6 Klebsiella pneumoniae Enterobacterales Carbapenems            40          18
#> # ℹ 1 more variable: marginal_resistance <dbl>
co_res <- compute_pairwise_coresistance(marg)
#> 'Escherichia coli': 3 classes, 40 isolates -- co-resistance matrix built.
#> 'Klebsiella pneumoniae': 3 classes, 40 isolates -- co-resistance matrix built.
co_res[["Klebsiella pneumoniae"]]$prevalence
#>                    3GC Carbapenems Fluoroquinolones
#> 3GC                 NA       0.225            0.325
#> Carbapenems      0.225          NA            0.225
#> Fluoroquinolones 0.325       0.225               NA
```
