# Plot Death vs Discharged Counts for Top Pathogens

Produces a side-by-side (dodged) bar chart comparing Death and
Discharged patient counts for the top `n` organisms. LAMA and Referred
are excluded – only the two definitive outcomes are shown.

## Usage

``` r
plot_death_discharged(
  data,
  n = 5,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  resistance_filter = NULL,
  patient_col = "PatientInformation_id",
  organism_col = "organism_name",
  outcome_col = "final_outcome",
  center_col = "center_name",
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  death_label = "Death",
  discharged_label = "Discharged",
  colours = NULL,
  bar_width = 0.65,
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

  Integer. Number of top organisms to show. Default 5.

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`. Default
  `"faceted"`.

- center:

  Character. Required when `mode = "single"`. Exact centre name.

- resistance_filter:

  Character or `NULL`. When `"R"`, only patients who had at least one
  Resistant result for the organism are included (worst-phenotype rule:
  any R = R). When `"S"`, only patients with exclusively Susceptible
  results are included. `NULL` (default) includes all patients
  regardless of resistance status.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- organism_col:

  Character. Organism name column. Default `"organism_name"`.

- outcome_col:

  Character. Final outcome column. Default `"final_outcome"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- antibiotic_col:

  Character. Antibiotic name column. Only used when `resistance_filter`
  is not `NULL`. Default `"antibiotic_name"`.

- value_col:

  Character. Susceptibility result column (R/S). Only used when
  `resistance_filter` is not `NULL`. Default `"antibiotic_value"`.

- death_label:

  Character. Exact string used for death in the outcome column. Default
  `"Death"`.

- discharged_label:

  Character. Exact string used for discharged. Default `"Discharged"`.

- colours:

  Named character vector with keys matching `death_label` and
  `discharged_label`. Default
  `c("Death" = "#E74C3C", "Discharged" = "#2ECC71")`.

- bar_width:

  Numeric. Width of each individual bar (0-1). Default `0.65`.

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
  `"VAP"`). Requires `syndrome_col`. If `NULL` (default), all syndromes
  are included.

## Value

A `ggplot` object.

## Details

Each bar is labelled with the count and the percentage that outcome
represents out of all Death + Discharged patients for that organism
(e.g. `45 (32%)` on a Death bar means 45 patients died, which is 32% of
the organism's known outcomes). This makes the mortality rate
immediately readable without a separate calculation.

**Deduplication:** uses `distinct(centre, patient, organism, outcome)` –
one record per patient per organism per outcome – consistent with the
original EDA script.

**Top N organisms** are ranked by total unique patients with a known
outcome (Death or Discharged). For `mode = "faceted"` the ranking is
done independently per centre; for `"overall"` it is done globally.

Supports three display modes:

- **"faceted"** (default) – one panel per centre.

- **"overall"** – all centres pooled into one chart.

- **"single"** – one specific centre.
