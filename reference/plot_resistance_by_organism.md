# Resistance proportion for top N organisms

100 vs Susceptible patients for the top `n` organisms (ranked by total
distinct patient count). A patient x organism pair is classified as
Resistant if any antibiotic result for that pair is `"R"`
(worst-phenotype rule at patient x organism x antibiotic level).

## Usage

``` r
plot_resistance_by_organism(
  data,
  n = 5,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "patient_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_name",
  value_col = "result",
  center_col = "center_name",
  colours = NULL,
  base_size = 11,
  ncol = 3,
  title = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame (one row per patient x organism x antibiotic result).

- n:

  Number of top organisms to display.

- mode:

  One of `"faceted"`, `"overall"`, or `"single"`.

- center:

  Centre name used when `mode = "single"`.

- patient_col:

  Column for patient ID.

- organism_col:

  Column for organism name.

- antibiotic_col:

  Column for antibiotic name.

- value_col:

  Column for resistance result (`"R"` / `"S"`).

- center_col:

  Column for centre name.

- colours:

  Named character vector with `"R"` and `"S"` entries.

- base_size:

  Base font size for
  [`eda_theme()`](https://saketlab.github.io/anumaan/reference/eda_theme.md).

- ncol:

  Facet columns (faceted mode only).

- title:

  Optional title override.

- syndrome_col:

  Optional column used to pre-filter by syndrome.

- syndrome_name:

  Optional value of `syndrome_col` to keep.
