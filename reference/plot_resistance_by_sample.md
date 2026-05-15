# Plot Distribution of Antibiotic Resistance Across Sample Types

Produces a horizontal 100% stacked bar chart showing the proportion of
Resistant (R) and Susceptible (S) deduplicated AST results for the top
`n` sample types. Each bar sums to 100% and is labelled with the R and S
percentages.

## Usage

``` r
plot_resistance_by_sample(
  data,
  n = 10,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  sample_col = "sample_type",
  organism_col = "organism_name",
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
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

  Integer. Number of top sample types to show (ranked by deduplicated
  test count). Default 10.

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`. Default
  `"faceted"`.

- center:

  Character. Required when `mode = "single"`. Exact centre name.

- patient_col:

  Character. Patient ID column. Used for deduplication. Default
  `"PatientInformation_id"`.

- sample_col:

  Character. Sample/specimen type column. Default `"sample_type"`.

- organism_col:

  Character. Organism name column. Used for deduplication. Default
  `"organism_name"`.

- antibiotic_col:

  Character. Antibiotic name column. Default `"antibiotic_name"`.

- value_col:

  Character. Susceptibility result column (R/S values). Default
  `"antibiotic_value"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- colours:

  Named character vector for R and S bars. Default
  `c("R" = "#D73027", "S" = "#1A9850")`.

- ncol:

  Integer. Facet columns (`mode = "faceted"` only). Default 2.

- base_size:

  Numeric. Base font size. Default 14.

- title:

  Character. Custom title. Auto-generated if `NULL`. Default `NULL`.

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

**Deduplication (worst-phenotype):** before counting, rows are collapsed
by patient x centre x sample type x organism x antibiotic. If the same
antibiotic appears more than once for a patient-organism-sample
combination (data entry error or repeat test), any single R result marks
the episode as R. This prevents duplicate rows from inflating counts.

**Counting unit:** deduplicated AST test results – one record per
patient x organism x antibiotic x sample type combination.

**Top N ranking:** sample types are ranked by total deduplicated test
count, so the most-tested specimen types appear first.

**Default sample column:** `"sample_type"` – the raw column. Switch to
`"specimen_normalized"` for standardised values.

Supports three display modes:

- **"faceted"** (default) – top *n* sample types per centre
  independently, one panel per centre.

- **"overall"** – all centres pooled; top *n* globally.

- **"single"** – one specific centre.
