# Compute Resistance Profile Probabilities per Pathogen

For each pathogen \\k\\ with \\n_k\\ antibiotic classes, enumerates all
\\2^{n_k}\\ binary resistance profiles \\\delta \in \\0,1\\^{n_k}\\ and
estimates their probabilities by solving a simplex-constrained weighted
least-squares Quadratic Programme (GBD equation 7.5.1.3):

## Usage

``` r
compute_resistance_profiles(
  marginal_output,
  coresistance_output,
  pathogens = NULL,
  top_n_pathogens = NULL,
  exclude_near_zero = TRUE,
  top_n_classes = NULL,
  sigma_sq = 1,
  ridge = 1e-08,
  pathogen_col = "organism_name",
  antibiotic_class_col = "antibiotic_class",
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL,
  n_cores = 1L
)
```

## Arguments

- marginal_output:

  List returned by
  [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md).

- coresistance_output:

  List returned by
  [`compute_pairwise_coresistance()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_coresistance.md).

- pathogens:

  Character vector. Pathogen(s) to process. `NULL` (default) processes
  every pathogen in `marginal_output`.

- top_n_pathogens:

  Integer or `NULL` (default). When set, only the top `top_n_pathogens`
  pathogens ranked by total isolates tested (sum of `n_tested` across
  all antibiotic classes, descending) are processed. Applied after the
  `pathogens` argument filter. Useful for focusing on the most data-rich
  pathogens – e.g. `top_n_pathogens = 5` runs profiles for only the 5
  most-tested pathogens.

- exclude_near_zero:

  Logical. If `TRUE` (default), antibiotic classes that appear in
  `marginal_output$near_zero` for a given pathogen are excluded from
  profile enumeration.

- top_n_classes:

  Integer or `NULL` (default). When set, only the top `top_n_classes`
  antibiotic classes ranked by `n_tested` (descending) are kept per
  pathogen before profile enumeration. Useful for capping the
  combinatorial explosion (2^n profiles) for pathogens tested against
  many drug classes – e.g. `top_n_classes = 5` gives at most 32
  profiles. Applied after `exclude_near_zero`.

- sigma_sq:

  Positive numeric. Assumed variance for each constraint (uniform).
  Default `1`.

- ridge:

  Positive numeric. Ridge term added to the QP Hessian for numerical
  stability. Default `1e-8`.

- pathogen_col:

  Character. Column name for pathogens. Must match the column used in
  Steps 1-2. Default `"organism_name"`.

- antibiotic_class_col:

  Character. Column name for antibiotic classes. Default
  `"antibiotic_class"`.

- facility_col:

  Character or `NULL`. Column identifying the facility/site. When
  provided together with `facility_name`, filters
  `marginal_output$marginal` and `marginal_output$near_zero` to the
  specified facility before profile enumeration. For this to work,
  [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md)
  must have been called with the same `facility_col` argument. Both must
  be supplied together or both left `NULL`. Default `NULL`.

- facility_name:

  Character or `NULL`. Facility value to retain. Default `NULL`.

- outcome_col:

  Character or `NULL`. Column containing patient outcomes. When provided
  together with `outcome_value`, filters `marginal_output$marginal` and
  `marginal_output$near_zero` to the specified outcome before profile
  enumeration.
  [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md)
  must have been called with the same `outcome_col` argument. Both must
  be supplied together or both left `NULL`. Default `NULL`.

- outcome_value:

  Character or `NULL`. Outcome value to retain (e.g. `"discharged"`,
  `"dead"`). Default `NULL`.

- n_cores:

  Integer. Number of CPU cores for parallel computation. Default `1L`
  (sequential).

## Value

Named list, one entry per pathogen, each a list with:

- `profiles`:

  Data frame: `profile` (character label, e.g.\\ `"RSR"`), `probability`
  (\\\hat{p}\_\delta\\), and one binary (0/1) integer column per
  antibiotic class indicating whether that class is resistant (`1`) or
  susceptible (`0`) in each profile.

- `classes`:

  Character vector of antibiotic class names used (alphabetical; bit 0 =
  classes\[1\]).

- `n_classes`:

  Integer. Number of classes used.

- `constraint_residuals`:

  Named numeric vector of \\m_i^\top \hat{p} - v_i\\ for each
  constraint. Small absolute values indicate good constraint
  satisfaction.

## Details

\$\$ \hat{p} = \arg\min\_{p \in \Delta\_{2^n}} \sum\_{i=1}^{m}
\frac{(m_i^\top p - v_i)^2}{\sigma_i^2} \$\$

where \\m = n(n+1)/2\\ data-derived linear constraints encode **n
marginal** resistance rates and **n(n-1)/2 pairwise** co-resistance
rates, and \\\Delta\\ is the standard probability simplex.

### Constraint rows in M

- Marginal (rows 1 to n):

  Row \\d\\: \\M\_{d,\delta} = 1\\ iff \\\delta_d = 1\\ (i.e., class
  \\d\\ is resistant in profile \\\delta\\). Constraint: \\\sum\_\delta
  M\_{d,\delta}\\p\_\delta = \hat{r}\_{kd}\\.

- Pairwise (rows n+1 to m):

  Row \\(d_1,d_2)\\: \\M\_{d_1 d_2,\delta} = 1\\ iff \\\delta\_{d_1} = 1
  \land \delta\_{d_2} = 1\\. Constraint: \\\sum\_\delta M\_{d_1
  d_2,\delta}\\p\_\delta = \hat{r}\_{k,d_1 d_2}\\. When a pairwise
  estimate is unavailable (too few co-tested isolates), the product of
  marginals (independence assumption) is used as fallback.

The QP is solved via
[`quadprog::solve.QP`](https://rdrr.io/pkg/quadprog/man/solve.QP.html).
A small ridge term (`ridge`) is added to the Hessian to guarantee strict
positive-definiteness. On solver failure the pathogen gets a uniform
distribution over all profiles.

### Performance Notes

This function has been optimized for speed with vectorized profile
generation and label creation (10-100x faster than previous versions).
However, computational complexity is still exponential in the number of
classes:

- **n \<= 14**: Fast (seconds to minutes)

- **n = 15**: Moderate (minutes)

- **n \>= 16**: Slow and memory-intensive (use `top_n_classes`)

For large datasets with many pathogens, consider using `top_n_classes`
to limit each pathogen to its most-tested classes (e.g.,
`top_n_classes = 12`).

## Examples

``` r
if (FALSE) { # \dontrun{
marg <- compute_marginal_resistance(amr_clean)
co_res <- compute_pairwise_coresistance(marg)
rp <- compute_resistance_profiles(marg, co_res)

# All profiles and probabilities for K. pneumoniae
rp[["Klebsiella pneumoniae"]]$profiles

# Constraint residuals (quality check)
rp[["Klebsiella pneumoniae"]]$constraint_residuals

# Single pathogen
rp_kp <- compute_resistance_profiles(
  marg, co_res,
  pathogens = "Klebsiella pneumoniae"
)
} # }
```
