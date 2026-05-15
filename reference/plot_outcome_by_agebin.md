# Plot Final Outcome Proportions by Age Bin

Produces a horizontal 100% stacked bar chart showing how final outcomes
(Death, Discharged, LAMA, Referred) are distributed across age groups.
Each bar represents one age bin and sums to 100%. An `n=` label above
each bar shows the total patients in that age group.

## Usage

``` r
plot_outcome_by_agebin(
  data,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
  agebin_col = "Age_bin",
  outcome_col = "final_outcome",
  center_col = "center_name",
  merge_referred = TRUE,
  age_levels = NULL,
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

- agebin_col:

  Character. Age bin column. Default `"Age_bin"`.

- outcome_col:

  Character. Final outcome column. Default `"final_outcome"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- merge_referred:

  Logical. Recode `"Transferred to other hospital"` to `"Referred"`.
  Default `TRUE`.

- age_levels:

  Character vector or `NULL`. Custom ordering of age bin labels. If
  `NULL` (default), bins are sorted automatically (`<1` first, numeric
  ranges in order, open-ended bins last).

- palette:

  Named character vector. Outcome -\> fill colour mapping. Unmatched
  outcomes receive grey. Default uses red/green/grey/blue preset.

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

**Deduplication:** one record per patient per centre – first recorded
outcome and age bin per patient per centre. Each patient contributes
exactly one bar segment per centre.

**Age bin ordering:** bins are ordered numerically (`<1`, `1-5`, `5-10`,
..., `85+`) automatically. Supply `age_levels` to override with a custom
order.

Supports three display modes:

- **"faceted"** (default) – one panel per centre.

- **"overall"** – all centres pooled; one chart.

- **"single"** – one specific centre.
