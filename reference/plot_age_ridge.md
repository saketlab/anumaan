# Ridge / Density Plot of Patient Age

Visualises the distribution of patient age on the x-axis using either a
stacked ridge plot (`mode = "all"`) or a single-centre density curve
(`mode = "single"`).

## Usage

``` r
plot_age_ridge(
  data,
  mode = c("all", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  age_col = "age_years",
  center_col = "center_name",
  min_age = 0,
  max_age = 120,
  colours = NULL,
  scale = 1.2,
  alpha = 0.7,
  base_size = 14,
  title = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  A data frame.

- mode:

  Character. One of `"all"` (ridge plot, all centres stacked on y-axis
  ordered by median age ascending) or `"single"` (density curve for one
  centre specified by `center`).

- center:

  Character. Centre name required when `mode = "single"`.

- patient_col:

  Character. Column of patient identifiers.

- age_col:

  Character. Column of patient age (numeric, in years). Default
  `"age_years"`.

- center_col:

  Character. Column of centre identifiers.

- min_age:

  Numeric. Minimum age to include. Default `0`.

- max_age:

  Numeric. Maximum age to include. Default `120`.

- colours:

  Named character vector of colours for centres. If `NULL` (default),
  the RColorBrewer "Set2" palette is used.

- scale:

  Numeric. `scale` argument passed to
  [`ggridges::geom_density_ridges()`](https://wilkelab.org/ggridges/reference/geom_density_ridges.html).
  Controls overlap between ridges. Default `1.2`.

- alpha:

  Numeric. Fill transparency (0-1). Default `0.7`.

- base_size:

  Numeric. Base font size. Default `14`.

- title:

  Character or `NULL`. Plot title. Auto-generated if `NULL`.

- syndrome_col:

  Character or `NULL`. Column of syndrome labels for pre-filtering.
  Default `NULL`.

- syndrome_name:

  Character or `NULL`. Syndrome value to retain. Requires
  `syndrome_col`. Default `NULL`.

## Value

A `ggplot` object.

## Details

One row per unique patient x centre is retained before plotting. Rows
with age outside `[min_age, max_age]` are dropped.
