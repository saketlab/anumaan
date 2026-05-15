# Boxplot of Length of Stay by Age Group

Displays the distribution of LOS (y-axis, log10 scale) across patient
age groups (x-axis) using boxplots. Supports `"faceted"` (one panel per
centre), `"overall"` (all centres pooled), and `"single"` (one centre)
modes.

## Usage

``` r
plot_los_by_agebin(
  data,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  agebin_col = "Age_bin",
  age_col = "age_years",
  age_breaks = NULL,
  age_labels = NULL,
  center_col = "center_name",
  admission_col = "admission_date",
  discharge_col = "final_outcome_date",
  min_los = 1,
  max_los = 365,
  age_levels = NULL,
  fill_colour = "#2C7FB8",
  fill_alpha = 0.7,
  outlier_alpha = 0.25,
  box_width = 0.6,
  ncol = 3,
  base_size = 14,
  title = NULL,
  subtitle = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  A data frame.

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`.

- center:

  Character. Centre name required when `mode = "single"`.

- patient_col:

  Character. Column of patient identifiers.

- agebin_col:

  Character. Pre-existing age group column. Used when `age_breaks` is
  `NULL`. Default `"Age_bin"`.

- age_col:

  Character. Continuous age column used when `age_breaks` is provided.
  Default `"age_years"`.

- age_breaks:

  Numeric vector or `NULL`. Breakpoints for deriving custom bins from
  `age_col` (e.g. `c(0, 5, 15, 30, 50, 70, Inf)`). `NULL` = use
  `agebin_col` directly.

- age_labels:

  Character vector or `NULL`. Labels for each bin created by
  `age_breaks`. Must have length `length(age_breaks) - 1`. `NULL` =
  auto-generate from breaks.

- center_col:

  Character. Column of centre identifiers.

- admission_col:

  Character. Admission date column. Default `"admission_date"`.

- discharge_col:

  Character. Discharge / final outcome date column. Default
  `"final_outcome_date"`.

- min_los:

  Numeric. Minimum LOS (days) to include. Default `1`.

- max_los:

  Numeric. Maximum LOS (days) to include. Default `365`.

- age_levels:

  Character vector or `NULL`. Custom bin order when using `agebin_col`.
  `NULL` = auto-sort.

- fill_colour:

  Character. Box fill colour. Default `"#2C7FB8"`.

- fill_alpha:

  Numeric. Box fill transparency (0-1). Default `0.7`.

- outlier_alpha:

  Numeric. Outlier point transparency (0-1). Default `0.25`.

- box_width:

  Numeric. Width of each box. Default `0.6`.

- ncol:

  Integer. Number of facet columns (`mode = "faceted"` only). Default
  `3`.

- base_size:

  Numeric. Base font size. Default `14`.

- title:

  Character or `NULL`. Plot title. Auto-generated if `NULL`.

- subtitle:

  Character or `NULL`. Plot subtitle. `NULL` = none.

- syndrome_col:

  Character or `NULL`. Syndrome column for pre-filtering. Default
  `NULL`.

- syndrome_name:

  Character or `NULL`. Syndrome value to retain. Requires
  `syndrome_col`. Default `NULL`.

## Value

A `ggplot` object.

## Details

LOS is computed as
`difftime(discharge_col, admission_col, units = "days")`. One row per
unique patient x centre x admission x discharge is retained. Rows with
`LOS < min_los` or `LOS > max_los` are dropped.

**Two ways to supply age groups:**

1.  **Pre-existing column** (default) – set `agebin_col` to an existing
    categorical column (e.g. `"Age_bin"`). Bins are sorted automatically
    (`<1`, `1-5`, ..., `85+`); override order with `age_levels`.

2.  **Custom breaks from continuous age** – set `age_breaks` to a
    numeric vector (e.g. `c(0, 5, 15, 30, 50, 70, Inf)`) and `age_col`
    to the continuous age column (default `"age_years"`). Bins are
    derived via [`cut()`](https://rdrr.io/r/base/cut.html) and labelled
    automatically. `agebin_col` is ignored when `age_breaks` is
    provided.
