# Plot Patient Distribution by Location Type Across Facilities

Produces a dodged bar chart showing the number of unique patients in
each location type (ICU, Ward, Other) per facility. Raw location strings
are normalised automatically via case-insensitive pattern matching, so
centre-specific ICU names (e.g. `"Respiratory ICU (RICU)"`) all map to
`"ICU"`.

## Usage

``` r
plot_location_by_facility(
  data,
  patient_col = "PatientInformation_id",
  location_col = "location",
  center_col = "center_name",
  icu_pattern = "icu",
  ward_pattern = "ward",
  icu_label = "ICU",
  ward_label = "Ward",
  other_label = "Other",
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

- location_col:

  Character. Location/unit column. Default `"location"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- icu_pattern:

  Character. Case-insensitive regex pattern that identifies ICU values.
  Default `"icu"`.

- ward_pattern:

  Character. Case-insensitive regex pattern that identifies Ward values.
  Default `"ward"`.

- icu_label:

  Character. Display label for ICU bars. Default `"ICU"`.

- ward_label:

  Character. Display label for Ward bars. Default `"Ward"`.

- other_label:

  Character. Display label for all other locations. Default `"Other"`.

- colours:

  Named character vector. Fill colours keyed by `icu_label`,
  `ward_label`, and `other_label`. Default pink/green/grey.

- bar_width:

  Numeric. Width of each individual bar. Default `0.68`.

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

**Normalisation:** a location value is classified as:

- **ICU** – if the lower-cased string contains `icu_pattern` (default
  `"icu"`)

- **Ward** – else if it contains `ward_pattern` (default `"ward"`)

- **Other** – everything else

**Deduplication:** `distinct(patient, centre, location type)` – each
unique patient-centre-location combination is counted once.
