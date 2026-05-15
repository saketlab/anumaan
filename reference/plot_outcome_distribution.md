# Plot Distribution of Final Outcomes

Produces a bar chart of final outcome counts per centre, with each bar
labelled with the count and percentage of that centre's (or overall)
unique patients. Bars are coloured by outcome category.

## Usage

``` r
plot_outcome_distribution(
  data,
  mode = c("faceted", "overall", "single"),
  center = NULL,
  patient_col = "PatientInformation_id",
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

- mode:

  Character. One of `"faceted"`, `"overall"`, or `"single"`. Default
  `"faceted"`.

- center:

  Character. Required when `mode = "single"`. Exact centre name.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- outcome_col:

  Character. Final outcome column. Default `"final_outcome"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- merge_referred:

  Logical. Recode `"Transferred to other hospital"` to `"Referred"`.
  Default `TRUE`.

- palette:

  Named character vector mapping outcome values to colours. Any outcome
  not in the vector gets a default ggplot colour. Default uses a preset
  palette for common outcome categories.

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

**Deduplication:** a patient may appear multiple times in a long-format
dataset (one row per isolate or antibiotic). This function takes the
*first* recorded outcome per patient per centre so each patient is
counted once.

**Referred merging:** by default, `"Transferred to other hospital"` is
recoded to `"Referred"` since they represent the same clinical event
across different centres. Set `merge_referred = FALSE` to keep them
separate.

Supports three display modes:

- **"faceted"** (default) – one panel per centre.

- **"overall"** – all centres pooled into one chart.

- **"single"** – one specific centre; pass the name via `center`.
