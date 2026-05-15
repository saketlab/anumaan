# Plot Antibiotic Resistance Heatmap

Produces a ggplot2 tile heatmap where:

- **x-axis** – antibiotics (ordered by total patients tested, most
  tested on the left)

- **y-axis** – centres

- **fill colour** – proportion of patients who were Resistant for that
  antibiotic at that centre (0 = fully susceptible, 1 = fully resistant)

## Usage

``` r
plot_abx_heatmap(
  data,
  n = NULL,
  mode = c("all", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  organism_col = "organism_name",
  center_col = "center_name",
  show_values = TRUE,
  midpoint = 0.5,
  low_colour = "#1a9850",
  mid_colour = "#fee08b",
  high_colour = "#d73027",
  base_size = 14,
  title = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame. Long-format AMR dataset.

- n:

  Integer or `NULL`. If an integer, only the top *n* antibiotics by
  total patients tested are shown. If `NULL`, all antibiotics are shown.
  Default `NULL`.

- mode:

  Character. `"all"` or `"single"`. Default `"all"`.

- center:

  Character. Required when `mode = "single"`. Exact centre name.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- antibiotic_col:

  Character. Antibiotic name column. Default `"antibiotic_name"`.

- value_col:

  Character. Susceptibility result column (R/S values). Default
  `"antibiotic_value"`.

- organism_col:

  Character. Organism name column. Default `"organism_name"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- show_values:

  Logical. Print proportion values inside tiles. Default `TRUE`.

- midpoint:

  Numeric. Midpoint of the colour gradient (0-1). Default `0.5`.

- low_colour:

  Character. Colour for proportion = 0 (fully susceptible). Default
  `"#1a9850"` (green).

- mid_colour:

  Character. Colour for the midpoint. Default `"#fee08b"` (yellow).

- high_colour:

  Character. Colour for proportion = 1 (fully resistant). Default
  `"#d73027"` (red).

- base_size:

  Numeric. Base font size. Default 14.

- title:

  Character. Custom plot title. Auto-generated if `NULL`. Default
  `NULL`.

- syndrome_col:

  Character or `NULL`. Column name containing syndrome/infection
  category labels (e.g. `"infectious_syndrome"`). If `NULL` (default),
  no syndrome filtering is applied.

- syndrome_name:

  Character or `NULL`. The syndrome value to retain (e.g. `"BSI"`,
  `"VAP"`). Requires `syndrome_col` to be set. If `NULL` (default), all
  syndromes are included.

## Value

A `ggplot` object.

## Details

Only R and S results are used. The same **worst-phenotype rule** as
[`plot_abx_susceptibility()`](https://saketlab.github.io/anumaan/reference/plot_abx_susceptibility.md)
is applied: per patient x organism x antibiotic, any single R marks the
episode as R.

Cells with no data (antibiotic not tested at that centre) are shown in
light grey.

Supports two display modes:

- **"all"** (default) – all centres on the y-axis.

- **"single"** – filter to one centre before plotting (produces a
  single-row heatmap useful for per-centre reports).
