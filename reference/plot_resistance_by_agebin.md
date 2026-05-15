# Resistance proportion by age group

100 vs Susceptible patients within each age bin. A patient is classified
as Resistant if they have ANY (organism × antibiotic) result of `"R"`;
Susceptible only when every result is `"S"` (worst-phenotype rule
applied at patient × organism × antibiotic level, then collapsed to
patient level).

## Usage

``` r
plot_resistance_by_agebin(
  data,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "patient_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_name",
  value_col = "result",
  agebin_col = "Age_bin",
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

  Data frame (one row per patient × organism × antibiotic result).

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

- agebin_col:

  Column for age bin label.

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
