# Compute Baseline Mortality Rate Among Fully Susceptible Patients (p0)

Computes the baseline death probability \\p_0\\ among patients who are
fully susceptible – i.e. **all** antibiotic test results are `"S"`. Any
patient with at least one `"R"` result is excluded. Used as the
denominator correction for OR -\> RR conversion (Zhang & Yu 1998).

## Usage

``` r
compute_p0(
  data,
  patient_col,
  antibiotic_value_col,
  outcome_col,
  death_value = "Death",
  resistant_value = "R",
  syndrome_col = NULL,
  syndrome_name = NULL,
  facility_col = NULL,
  facility_name = NULL
)
```

## Arguments

- data:

  Data frame. One row per patient x antibiotic test.

- patient_col:

  Character. Unique patient identifier column.

- antibiotic_value_col:

  Character. Antibiotic result column (`"S"`, `"R"`, `"I"`).

- outcome_col:

  Character. Final outcome column.

- death_value:

  Character. Death indicator. Default `"Death"`.

- resistant_value:

  Character. Resistance indicator used to exclude patients. Default
  `"R"`.

- syndrome_col:

  Character or `NULL`.

- syndrome_name:

  Character or `NULL`.

- facility_col:

  Character or `NULL`. Facility identifier column. Required when
  `facility_name` is specified.

- facility_name:

  Character or `NULL`. If provided, filters data to the specified
  facility before computing p0.

## Value

Named list: `p0` (scalar), `n_susceptible`, `n_susceptible_deaths`,
`summary` (one-row data frame).

## References

Zhang J, Yu KF. What's the relative risk? JAMA. 1998;280(19):1690-1.
