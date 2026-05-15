# Plot Distribution of Final Outcomes by Year

Produces a stacked bar chart where the x-axis shows the year derived
from `date_col`, the y-axis shows the number of unique patients, and
bars are filled by final outcome category. Each bar segment is labelled
with count and percentage; a total `n=` label sits above each bar.

## Usage

``` r
plot_outcome_by_year(
  data,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  outcome_col = "final_outcome",
  date_col = "final_outcome_date",
  center_col = "center_name",
  merge_referred = TRUE,
  palette = c(Death = "#E74C3C", Died = "#E74C3C", Discharged = "#2ECC71", LAMA =
    "#95A5A6", Referred = "#3498DB", `Transferred to other hospital` = "#3498DB"),
  ncol = 2,
  base_size = 14,
  title = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame. Long-format AMR dataset.

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`. Default
  `"faceted"`.

- center:

  Character. Required when `mode = "single"`. Exact centre name.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- outcome_col:

  Character. Final outcome column. Default `"final_outcome"`.

- date_col:

  Character. Date column used to extract the year (must be a `Date` or
  `POSIXct` column). Default `"final_outcome_date"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- merge_referred:

  Logical. Recode `"Transferred to other hospital"` to `"Referred"`.
  Default `TRUE`.

- palette:

  Named character vector mapping outcome values to colours. Unmatched
  outcomes get a default ggplot colour.

- ncol:

  Integer. Facet columns (`mode = "faceted"` only). Default 2.

- base_size:

  Numeric. Base font size. Default 14.

- title:

  Character. Custom plot title. Auto-generated if `NULL`. Default
  `NULL`.

- syndrome_col:

  Character or `NULL`. Syndrome filter column. Default `NULL`.

- syndrome_name:

  Character or `NULL`. Syndrome value to retain. Requires
  `syndrome_col`. Default `NULL`.

## Value

A `ggplot` object.

## Details

One outcome is taken per patient per centre (the first recorded value)
to avoid double-counting long-format rows.

Supports three display modes:

- **"faceted"** (default) – one panel per centre.

- **"overall"** – all centres pooled into a single chart.

- **"single"** – one specific centre; pass its name via `center`.
