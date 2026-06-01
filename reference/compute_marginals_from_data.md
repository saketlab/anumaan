# Compute Marginal Resistance Rates from Preprocessed Wide Data

Takes the wide-format tibble produced by
[`preprocess_for_profiles()`](https://saketlab.github.io/anumaan/reference/preprocess_for_profiles.md)
and computes marginal resistance rates per pathogen x antibiotic class,
optionally stratified by geography, year, and/or patient outcome. The
counting unit is always `isolate_id`: a patient with multiple isolates
contributes one count per isolate.

## Usage

``` r
compute_marginals_from_data(
  data_wide,
  col_map,
  panel_map,
  stratify_by = NULL,
  outcome_col = NULL,
  min_n_tested = 30L,
  external_marginals = NULL,
  ext_col_map = list(pathogen_col = "pathogen", class_col = "antibiotic_class",
    geography_col = "geography", year_col = "year", rate_col = "resistance_prevalence")
)
```

## Arguments

- data_wide:

  Tibble. Output of `preprocess_for_profiles()$data_wide`. One row per
  isolate; antibiotic-class columns contain `"S"`, `"R"`, or `NA` (not
  tested).

- col_map:

  Named list. `col_map_resolved` from
  [`preprocess_for_profiles()`](https://saketlab.github.io/anumaan/reference/preprocess_for_profiles.md).

- panel_map:

  Named list. Same panel used in
  [`preprocess_for_profiles()`](https://saketlab.github.io/anumaan/reference/preprocess_for_profiles.md).

- stratify_by:

  Character vector or `NULL`. Dimensions to stratify over. Valid:
  `"geography"`, `"year"`. Default `NULL`.

- outcome_col:

  Character or `NULL`. Column name for patient outcome. When supplied,
  separate marginals are computed for each outcome value (e.g. Died vs
  Discharged). Default `NULL`.

- min_n_tested:

  Integer. Minimum tested-isolate count per (stratum x pathogen x class)
  cell. Cells below this threshold are dropped with a logged reason.
  Default `30L`.

- external_marginals:

  Data frame or `NULL`. Pre-modelled marginal resistance rates (e.g. GBD
  ST-GPR estimates). When provided, locally computed rates are replaced
  where a match exists on pathogen x class \[x geography x year\]. Must
  have columns matching `ext_col_map`. Default `NULL`.

- ext_col_map:

  Named list. Column names in `external_marginals`. Default
  `list(pathogen_col="pathogen", class_col="antibiotic_class", geography_col="geography", year_col="year", rate_col="resistance_prevalence")`.

## Value

Tibble with columns: `pathogen`, `antibiotic_class`, \[stratum
columns\], `n_tested`, `n_resistant`, `marginal_resistance`,
`marginal_source` (`"computed"` or `"external"`).

## Details

Optionally, externally modelled marginals (e.g. from GBD ST-GPR) can be
supplied to override the locally computed values while retaining the
pairwise correlation structure from the local data.
