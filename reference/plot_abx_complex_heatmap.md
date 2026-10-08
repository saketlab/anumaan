# Plot Antibiotic Resistance ComplexHeatmap

Produces a clustered
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
of the proportion of patients Resistant, with **centres on rows** and
**antibiotics on columns**. Unlike
[`plot_abx_heatmap`](https://saketlab.github.io/anumaan/reference/plot_abx_heatmap.md)
(a `ggplot2` tile heatmap with a fixed row/column order), this version
hierarchically clusters both rows and columns by default, which groups
centres with similar resistance profiles and antibiotics with similar
resistance patterns together.

## Usage

``` r
plot_abx_complex_heatmap(
  data,
  patient_col = "PatientInformation_id",
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  organism_col = "organism_name",
  center_col = "center_name",
  low_colour = "#1a9850",
  mid_colour = "#fee08b",
  high_colour = "#d73027",
  midpoint = 0.5,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  legend_name = "Resistance",
  row_title = "Centres",
  column_title = "Antibiotics",
  base_size = 14,
  column_names_max_height = grid::unit(10, "cm"),
  save_path = NULL,
  width = 15,
  height = 8,
  dpi = 300
)
```

## Arguments

- data:

  Data frame. Long-format AMR dataset.

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

- low_colour:

  Character. Colour for proportion = 0 (fully susceptible). Default
  `"#1a9850"` (green).

- mid_colour:

  Character. Colour for the midpoint. Default `"#fee08b"` (yellow).

- high_colour:

  Character. Colour for proportion = 1 (fully resistant). Default
  `"#d73027"` (red).

- midpoint:

  Numeric. Midpoint of the colour gradient (0-1). Default `0.5`.

- cluster_rows:

  Logical. Hierarchically cluster centres. Default `TRUE`.

- cluster_columns:

  Logical. Hierarchically cluster antibiotics. Default `TRUE`.

- legend_name:

  Character. Heatmap legend title. Default `"Resistance"`.

- row_title:

  Character. Row-axis title. Default `"Centres"`.

- column_title:

  Character. Column-axis title. Default `"Antibiotics"`.

- base_size:

  Numeric. Base font size for titles/labels. Default 14.

- column_names_max_height:

  [`grid::unit`](https://rdrr.io/r/grid/unit.html). Vertical space
  reserved for the (45-degree-rotated) column-name labels below the
  heatmap body. Increase this (and/or `height`) if long labels – e.g.
  antibiotic class names – are being clipped. Default
  `grid::unit(10, "cm")`.

- save_path:

  Character or `NULL`. If supplied, the heatmap is drawn to a PNG at
  this path (via
  [`grDevices::png()`](https://rdrr.io/r/grDevices/png.html) +
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html))
  in addition to being returned. Default `NULL` (not saved).

- width:

  Numeric. PNG width in inches. Used only when `save_path` is supplied.
  Default 15.

- height:

  Numeric. PNG height in inches. Used only when `save_path` is supplied.
  Default 8.

- dpi:

  Numeric. PNG resolution. Used only when `save_path` is supplied.
  Default 300.

## Value

A
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object (draw with `ComplexHeatmap::draw(ht)` if not saving to file).

## Details

Only R and S results are used. The same **worst-phenotype rule** as
[`plot_abx_heatmap()`](https://saketlab.github.io/anumaan/reference/plot_abx_heatmap.md)
is applied: per patient x organism x antibiotic, any single R marks the
episode as R. Centre/antibiotic combinations with no tests are shown as
0 (fully susceptible) in the underlying matrix – there is no "no data"
grey cell in this version, since `Heatmap()` requires a complete numeric
matrix.

Requires the Bioconductor package `ComplexHeatmap` and the CRAN package
`circlize`, neither of which is a hard dependency of `anumaan`. Install
with:
`if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("ComplexHeatmap")`
and `install.packages("circlize")`.
