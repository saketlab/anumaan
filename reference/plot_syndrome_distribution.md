# Plot Syndrome Distribution

Produces a column chart showing the number of unique patients per
infectious syndrome. Supports three modes: pooled across all centres
(`"overall"`), faceted by centre (`"faceted"`), or a single specified
centre (`"single"`).

## Usage

``` r
plot_syndrome_distribution(
  data,
  mode = c("overall", "faceted", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  syndrome_col = "infectious_syndrome",
  center_col = "center_name",
  colour = "#4393C3",
  ncol = 2,
  base_size = 14,
  title = NULL
)
```

## Arguments

- data:

  Data frame. Must contain the columns named by `patient_col` and
  `syndrome_col`. In `"faceted"` and `"single"` modes, `center_col` is
  also required.

- mode:

  Character. One of `"overall"` (default), `"faceted"`, or `"single"`.

- center:

  Character or `NULL`. Centre name to display when `mode = "single"`.
  Ignored otherwise.

- patient_col:

  Character. Column containing patient identifiers. Default
  `"PatientInformation_id"`.

- syndrome_col:

  Character. Column containing infectious syndrome labels. Default
  `"infectious_syndrome"`.

- center_col:

  Character. Column containing hospital / centre names. Default
  `"center_name"`.

- colour:

  Character. Fill colour for the bars. Default `"#4393C3"`.

- ncol:

  Integer. Number of columns in faceted layout. Default `2`.

- base_size:

  Numeric. Base font size. Default `14`.

- title:

  Character. Custom plot title. Auto-generated if `NULL`.

## Value

A `ggplot` object.
