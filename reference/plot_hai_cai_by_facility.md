# Plot HAI vs CAI Infection Distribution by Facility

Derives Healthcare-Associated Infection (HAI) vs Community-Associated
Infection (CAI) status from dates, then produces a stacked bar chart
showing unique patient counts per centre. Each bar carries count and
percentage labels inside each segment, and an `n=` total label above.

## Usage

``` r
plot_hai_cai_by_facility(
  data,
  patient_col = "PatientInformation_id",
  center_col = "center_name",
  admission_col = "admission_date",
  culture_col = "culture_date",
  hai_threshold_hours = 48,
  hai_label = "HAI",
  cai_label = "CAI",
  colours = NULL,
  base_size = 14,
  title = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame. Long-format AMR dataset.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- admission_col:

  Character. Admission date column. Default `"admission_date"`.

- culture_col:

  Character. Culture/specimen date column. Default `"culture_date"`.

- hai_threshold_hours:

  Numeric. Hours from admission to first culture at or above which a
  case is classified as HAI. Default `48`.

- hai_label:

  Character. Label for HAI bars. Default `"HAI"`.

- cai_label:

  Character. Label for CAI bars. Default `"CAI"`.

- colours:

  Named character vector. Fill colours keyed by `hai_label` and
  `cai_label`. Default red/green.

- base_size:

  Numeric. Base font size. Default `14`.

- title:

  Character. Custom title. Auto-generated if `NULL`.

- syndrome_col:

  Character or `NULL`. Column name containing syndrome/infection
  category labels. `NULL` = no filter.

- syndrome_name:

  Character or `NULL`. Syndrome value to retain. Requires
  `syndrome_col`.

## Value

A `ggplot` object.

## Details

**Classification rule:** for each patient x centre, the earliest culture
date is used. If the time from admission to that culture is \\\geq\\
`hai_threshold_hours` (default 48 h), the patient is classified as HAI;
otherwise CAI. The existing `type_of_infection` column in the dataset is
*ignored* – classification is derived entirely from dates.

**Deduplication:** one classification per patient per centre (earliest
culture date).
