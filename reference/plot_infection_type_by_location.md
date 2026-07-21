# Plot Infection Type (HAI/CAI) by Location (ICU/Ward), Faceted by Hospital

100 types (e.g. HAI vs CAI) among patients in each location (e.g. ICU vs
Ward). Each bar is normalised *within* its location, so the two bars in
a panel directly compare \\P(\text{infection type} \mid
\text{location})\\ – i.e. whether infection type is associated with
location, per hospital.

## Usage

``` r
plot_infection_type_by_location(
  data,
  patient_col = "PatientInformation_id",
  center_col = "center_name",
  location_col = "location",
  infection_col = "infection_type",
  location_levels = c("ICU", "Ward"),
  infection_levels = c("HAI", "CAI"),
  style = c("bar", "heatmap"),
  colours = NULL,
  ncol = 5,
  base_size = 12,
  show_counts = TRUE,
  title = NULL
)
```

## Arguments

- data:

  Data frame. Long-format AMR / cleaned stewardship data.

- patient_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- center_col:

  Character. Centre/facility column. Default `"center_name"`.

- location_col:

  Character. Location column. Default `"location"`.

- infection_col:

  Character. Infection-type column. Default `"infection_type"`.

- location_levels:

  Character vector. Location values to keep (become the x-axis bars).
  Default `c("ICU", "Ward")`.

- infection_levels:

  Character vector. Infection-type values to keep (become the stacked
  fills). Default `c("HAI", "CAI")`.

- style:

  Character. `"bar"` (default) draws a 100 bar per hospital (x =
  location, fill = infection type, normalised within location).
  `"heatmap"` draws a tile grid per hospital (x = location, y =
  infection type, fill = within-hospital joint proportion).

- colours:

  Named character vector of fill colours keyed by `infection_levels`.
  Auto-generated if `NULL`.

- ncol:

  Integer. Number of facet columns. Default `5`.

- base_size:

  Numeric. Base font size. Default `12`.

- show_counts:

  Logical. Print unique-patient counts inside each segment. Default
  `TRUE`.

- title:

  Character. Custom title. Auto-generated if `NULL`.

## Value

A `ggplot` object (one panel per centre).

## Details

**Deduplication:** unique patients per centre x location x infection
type. Only the requested `location_levels` and `infection_levels` are
kept (other/blank/compound values are dropped).
