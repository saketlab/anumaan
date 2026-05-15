# Plot Unique Patient Count by Hospital

Produces a column chart showing the number of unique patients per
hospital. Bars are sorted in descending order (highest count on the
left) and each bar is labelled with its count.

## Usage

``` r
plot_patients_by_hospital(
  data,
  patient_col = "PatientInformation_id",
  center_col = "center_name",
  colour = "#4393C3",
  bar_width = 0.68,
  base_size = 14,
  title = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame. Must contain the columns named by `patient_col` and
  `center_col`.

- patient_col:

  Character. Column containing patient identifiers. Default
  `"PatientInformation_id"`.

- center_col:

  Character. Column containing hospital / centre names. Default
  `"center_name"`.

- colour:

  Character. Fill colour for the bars. Default `"#4393C3"`.

- bar_width:

  Numeric. Width of bars (0-1). Default `0.68`.

- base_size:

  Numeric. Base font size. Default `14`.

- title:

  Character. Custom plot title. Auto-generated if `NULL`.

- syndrome_col:

  Character or `NULL`. Syndrome filter column. Default `NULL`.

- syndrome_name:

  Character or `NULL`. Syndrome value to retain. Requires
  `syndrome_col`. Default `NULL`.

## Value

A `ggplot` object.
