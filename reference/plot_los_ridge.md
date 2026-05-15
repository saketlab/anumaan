# Ridge / Density Plot of Length of Stay (LOS)

Visualises the distribution of length of stay (LOS) on a log10 x-axis
using either a stacked ridge plot (`mode = "all"`) or a single-centre
density curve (`mode = "single"`).

## Usage

``` r
plot_los_ridge(
  data,
  mode = c("all", "single"),
  center = NULL,
  patient_col = "patient_id",
  center_col = "center_name",
  admission_col = "admission_date",
  discharge_col = "final_outcome_date",
  filter_outcome = NULL,
  outcome_col = "final_outcome",
  min_los = 1,
  max_los = 365,
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

  Character. One of `"all"` (ridge plot, all centres stacked on y-axis,
  ordered by median LOS ascending) or `"single"` (density curve for one
  centre specified by `center`).

- center:

  Character. Centre name required when `mode = "single"`.

- patient_col:

  Character. Column of patient identifiers.

- center_col:

  Character. Column of centre identifiers.

- admission_col:

  Character. Column containing admission date/datetime. Default
  `"admission_date"`.

- discharge_col:

  Character. Column containing discharge / final outcome date. Default
  `"final_outcome_date"`.

- filter_outcome:

  Character or `NULL`. If not `NULL`, only rows whose outcome column
  (`outcome_col`) equals this value are used (e.g. `"Discharged"`). Set
  `NULL` to include all outcomes. Default `NULL`.

- outcome_col:

  Character. Column of patient outcomes, used only when `filter_outcome`
  is not `NULL`. Default `"final_outcome"`.

- min_los:

  Numeric. Minimum LOS (days) to include. Rows below this are removed
  (prevents `log(0)` issues). Default `1`.

- max_los:

  Numeric. Maximum LOS (days) to include. Default `365`.

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

LOS is computed as
`difftime(discharge_col, admission_col, units = "days")`. One row per
unique patient x centre x admission_date x discharge_date is retained
before plotting. Rows with `LOS < min_los` or `LOS > max_los` are
dropped.
