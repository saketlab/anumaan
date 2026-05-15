# Plot Antibiotic Susceptibility Pattern (Stacked R/S Bars)

Produces a horizontal stacked bar chart showing the number of admissions
tested as Resistant (R) or Susceptible (S) for the top `n` antibiotics.
Percentage labels are shown inside each bar segment.

## Usage

``` r
plot_abx_susceptibility(
  data,
  n = 5,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  organism_col = "organism_name",
  center_col = "center_name",
  colours = c(R = "#D73027", S = "#1A9850"),
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

- n:

  Integer. Number of top antibiotics to show (ranked by total patients
  tested). Default 5.

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`. Default
  `"faceted"`.

- center:

  Character. Required when `mode = "single"`. Exact name of the centre
  to plot.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- antibiotic_col:

  Character. Antibiotic name column. Default `"antibiotic_name"`.

- value_col:

  Character. Susceptibility result column containing "R" and "S" values.
  Default `"antibiotic_value"`.

- organism_col:

  Character. Organism name column (used in worst-phenotype
  deduplication). Default `"organism_name"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- colours:

  Named character vector. Fill colours for R and S. Default
  `c("R" = "#D73027", "S" = "#1A9850")`.

- ncol:

  Integer. Facet columns (`mode = "faceted"` only). Default 2.

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

**Worst-phenotype rule:** when a patient has multiple records for the
same organism-antibiotic combination, any single R result marks the
whole episode as R. This prevents double-counting repeat cultures from
the same infection episode.

**Counting unit:** bars show unique patients (not isolate rows), so a
patient tested against the same antibiotic twice is counted once.

Supports three display modes:

- **"faceted"** – top *n* antibiotics per centre independently, one
  panel per centre.

- **"overall"** – all centres pooled; top *n* antibiotics globally.

- **"single"** – one specific centre; pass the centre name via `center`.
