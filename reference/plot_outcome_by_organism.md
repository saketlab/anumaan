# Plot Final Outcome Proportions for Resistant or Susceptible Patients

Produces a 100% horizontal stacked bar chart showing the distribution of
final outcomes (Death, Discharged, LAMA, Referred) for the top `n`
organisms, split by whether the patient had a Resistant or Susceptible
result.

## Usage

``` r
plot_outcome_by_organism(
  data,
  n = 5,
  resistance_filter = c("R", "S"),
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  organism_col = "organism_name",
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  outcome_col = "final_outcome",
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

- n:

  Integer. Number of top organisms to show. Default 5.

- resistance_filter:

  Character. `"R"` or `"S"`. Which patient group to display. Default
  `"R"`.

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`. Default
  `"faceted"`.

- center:

  Character. Required when `mode = "single"`. Exact centre name.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- organism_col:

  Character. Organism name column. Default `"organism_name"`.

- antibiotic_col:

  Character. Antibiotic name column. Default `"antibiotic_name"`.

- value_col:

  Character. Susceptibility result column (R/S). Default
  `"antibiotic_value"`.

- outcome_col:

  Character. Final outcome column. Default `"final_outcome"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- merge_referred:

  Logical. Recode `"Transferred to other hospital"` to `"Referred"`.
  Default `TRUE`.

- palette:

  Named character vector mapping outcome values to colours. Unmatched
  outcomes get grey automatically.

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

Each bar represents one organism. The fill shows outcome proportions. A
`n=` label is printed above each bar showing the total number of unique
patients in that organism x resistance group.

**Resistance grouping:** `resistance_filter = "R"` shows patients who
had *at least one* R result for the given organism. `"S"` shows patients
with only S results. Because these groups overlap (a patient resistant
to one antibiotic may be susceptible to another), calling the function
twice – once for "R" and once for "S" – gives the most complete picture.

**Top N organisms** are ranked by total unique patients per centre (or
globally for `mode = "overall"`), calculated before the R/S split so the
same organisms appear in both R and S plots.

Supports three display modes:

- **"faceted"** (default) – one panel per centre.

- **"overall"** – all centres pooled.

- **"single"** – one specific centre.
