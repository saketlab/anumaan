# Plot Top Pathogens by Unique Patients

Produces a horizontal bar chart of the top `n` organisms ranked by
unique patient count. Supports three display modes:

- **"faceted"** – one panel per centre, top *n* organisms within each
  centre independently. Facet strips show patient and organism counts
  for that centre.

- **"overall"** – all centres pooled; single chart of the top *n*
  organisms across the entire dataset.

- **"single"** – one specific centre; pass the centre name via `center`.

## Usage

``` r
plot_top_organisms(
  data,
  n = 5,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  organism_col = "organism_name",
  center_col = "center_name",
  fill_colour = "steelblue",
  ncol = 2,
  base_size = 14,
  title = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame. Long-format AMR dataset (one row per isolate/antibiotic
  record).

- n:

  Integer. Number of top organisms to show. Default 5.

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`. Default
  `"faceted"`.

- center:

  Character. Required when `mode = "single"`. Exact name of the centre
  to plot (must match a value in `center_col`).

- patient_col:

  Character. Column name for patient identifiers. Default
  `"PatientInformation_id"`.

- organism_col:

  Character. Column name for organism names. Default `"organism_name"`.

- center_col:

  Character. Column name for centre/facility names. Default
  `"center_name"`.

- fill_colour:

  Character. Bar fill colour (any R colour string or hex). Default
  `"steelblue"`.

- ncol:

  Integer. Number of columns in the facet grid (`mode = "faceted"`
  only). Default 2.

- base_size:

  Numeric. Base font size passed to
  [`eda_theme()`](https://saketlab.github.io/anumaan/reference/eda_theme.md).
  Default 14.

- title:

  Character. Custom plot title. If `NULL` a title is generated
  automatically from `mode` and `n`. Default `NULL`.

- syndrome_col:

  Character or `NULL`. Column name containing syndrome/infection
  category labels (e.g. `"infectious_syndrome"`). If `NULL` (default),
  no syndrome filtering is applied.

- syndrome_name:

  Character or `NULL`. The syndrome value to retain (e.g. `"BSI"`,
  `"VAP"`). Requires `syndrome_col` to be set. If `NULL` (default), all
  syndromes are included.

## Value

A `ggplot` object. Print it or pass to
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

Labels on each bar show absolute count and percentage of that centre's
(or overall) unique patient total.
