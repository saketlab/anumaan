# Plot Mono vs Polymicrobial Infections by Facility

Produces a dodged bar chart comparing the number of unique patients with
Monomicrobial vs Polymicrobial infections at each centre. Hospitals
appear on the x-axis; the two bars per hospital are coloured by
infection type. Count labels are printed above each bar.

## Usage

``` r
plot_mono_poly_by_facility(
  data,
  patient_col = "PatientInformation_id",
  poly_col = "is_polymicrobial",
  center_col = "center_name",
  mono_label = "Monomicrobial",
  poly_label = "Polymicrobial",
  colours = NULL,
  bar_width = 0.68,
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

- poly_col:

  Character. Column with binary infection type (`0` = mono, `1` = poly).
  Default `"is_polymicrobial"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- mono_label:

  Character. Label for monomicrobial bars. Default `"Monomicrobial"`.

- poly_label:

  Character. Label for polymicrobial bars. Default `"Polymicrobial"`.

- colours:

  Named character vector. Fill colours keyed by `mono_label` and
  `poly_label`. Default blue/red.

- bar_width:

  Numeric. Width of each individual bar. Default `0.68`.

- base_size:

  Numeric. Base font size. Default `14`.

- title:

  Character. Custom title. Auto-generated if `NULL`.

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

**Deduplication:** `distinct(patient, centre, infection type)` – each
unique patient-centre-infection type combination is counted once. A
patient who has both mono and polymicrobial episodes at the same centre
contributes to both bars, matching the original EDA script.

**Infection type:** the `poly_col` column is expected to contain binary
values (`0` = Monomicrobial, `1` = Polymicrobial). The labels shown in
the legend and plot are controlled by `mono_label` and `poly_label`.
