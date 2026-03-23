# Convert Odds Ratios to Relative Risks Using Baseline Mortality (p0)

Applies the Zhang & Yu (1998) formula to convert odds ratios from
logistic regression to approximate relative risks: \$\$RR =
\frac{OR}{(1 - p_0) + p_0 \times OR}\$\$ Applied to the point estimate
and both CI bounds. Works directly on the output of
[`fit_mortality_rr_logistic()`](https://saketlab.github.io/anumaan/reference/fit_mortality_rr_logistic.md).

## Usage

``` r
convert_or_to_rr(
  or_data,
  p0,
  or_col = "OR_death",
  ci_lower_col = "CI_lower",
  ci_upper_col = "CI_upper"
)
```

## Arguments

- or_data:

  Data frame with OR columns.

- p0:

  Numeric scalar. From `compute_p0()$p0`.

- or_col:

  Character. Default `"OR_death"`.

- ci_lower_col:

  Character. Default `"CI_lower"`.

- ci_upper_col:

  Character. Default `"CI_upper"`.

## Value

`or_data` with columns `RR_death`, `RR_lower`, `RR_upper` appended.

## References

Zhang J, Yu KF. What's the relative risk? JAMA. 1998;280(19):1690-1.
