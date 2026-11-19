# Exploratory Data Analysis in anumaan

``` r
library(anumaan)
library(dplyr)
```

## Overview

This vignette covers the **exploratory data analysis (EDA) layer**. It
sits between the two other pipelines: it is run on data that has already
been through
[Preprocessing](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md),
and it is what a stewardship team typically reviews *before* committing
to [DALY burden
estimation](https://saketlab.github.io/anumaan/articles/daly-pipeline.md).
It does not transform data (like Pipeline 1) or estimate burden (like
Pipeline 3) – it only visualises what is already there, to catch
data-quality problems and characterise the cohort before further
modelling.

The functions in this layer are named `plot_*` and share the styling
helper
[`eda_theme()`](https://saketlab.github.io/anumaan/reference/eda_theme.md).

Most `plot_*` functions here share a common `mode` argument:

- **`"faceted"`** – one panel per centre (the default for most
  functions)
- **`"overall"`** – all centres pooled into a single chart
- **`"single"`** – one specific centre, passed via `center = "..."`
- **`"all"`** – every centre in one chart, used by the heatmap and ridge
  plots in place of `"faceted"` and `"overall"`

and a common `syndrome_col`/`syndrome_name` pair for restricting the
plot to one infectious syndrome (e.g. `syndrome_name = "BSI"`) when that
column is available.

### Connecting from Pipeline 1

The `plot_*` functions’ *default* column names (`organism_name`,
`antibiotic_value`, `PatientInformation_id`, …) are convenience defaults
for raw-ish data – they do **not** automatically match the standardised
column names
[Preprocessing](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md)
actually produces. When passing real preprocessed output into these
functions, map the arguments to Pipeline 1’s real output columns
explicitly:

| EDA argument | Pipeline 1’s actual output column |
|----|----|
| `patient_col` | `patient_id` |
| `organism_col` | `organism_normalized` (from [`prep_standardize_organisms()`](https://saketlab.github.io/anumaan/reference/prep_standardize_organisms.md)) |
| `antibiotic_col` | `antibiotic_normalized` (from [`prep_standardize_antibiotics()`](https://saketlab.github.io/anumaan/reference/prep_standardize_antibiotics.md)) |
| `value_col` | `ast_value_harmonized` (from [`prep_harmonize_ast()`](https://saketlab.github.io/anumaan/reference/prep_harmonize_ast.md)) |
| `sample_col` | `specimen_normalized` (from [`prep_standardize_specimens()`](https://saketlab.github.io/anumaan/reference/prep_standardize_specimens.md)) |
| `agebin_col` | `Age_bin` (from [`prep_assign_age_bins()`](https://saketlab.github.io/anumaan/reference/prep_assign_age_bins.md)) |
| `age_col` | `Age` (from [`prep_fill_age()`](https://saketlab.github.io/anumaan/reference/prep_fill_age.md)) |
| `admission_col` | `date_of_admission` |
| `culture_col` | `date_of_culture` |
| `discharge_col` | `date_of_final_outcome` |
| `center_col` | `center_name` (a user-supplied raw column – Pipeline 1 does not invent hospital identity, it just carries it through) |
| `outcome_col` | `final_outcome` (the raw outcome text persists unchanged; [`prep_standardize_final_outcome()`](https://saketlab.github.io/anumaan/reference/prep_standardize_final_outcome.md) *adds* a separate `outcome_std` column alongside it rather than replacing it) |

Every example below uses these Pipeline-1-aligned names, with the
arguments passed explicitly rather than relying on defaults.

## A Small Synthetic Dataset

The examples below use a small synthetic long-format dataset (one row
per isolate x antibiotic) with the same column names Pipeline 1 would
hand off – as if it were the output of the
[Preprocessing](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md)
vignette’s `prepped` object, scaled up to 60 patients in two hospitals
so the `mode = "faceted"` examples have something to facet over.

``` r
set.seed(42)
n_pts <- 60
centers <- c("General Hospital", "City Medical Center")
organisms <- c(
  "Escherichia coli", "Klebsiella pneumoniae",
  "Staphylococcus aureus", "Pseudomonas aeruginosa"
)
antibiotics <- c("Amikacin", "Ciprofloxacin", "Meropenem", "Ceftriaxone")

# One row per patient: every patient-level value is drawn once, and derived
# columns (age bin, infection type, outcome date) follow from the others.
patients <- data.frame(
  patient_id = paste0("pt_", sprintf("%03d", 1:n_pts)),
  center_name = sample(centers, n_pts, replace = TRUE),
  organism_normalized = sample(organisms, n_pts, replace = TRUE),
  specimen_normalized = sample(c("Blood", "Urine", "Sputum"), n_pts, replace = TRUE),
  final_outcome = sample(
    c("Death", "Discharged", "LAMA", "Transferred to other hospital"),
    n_pts,
    replace = TRUE, prob = c(0.15, 0.65, 0.1, 0.1)
  ),
  is_polymicrobial = sample(c(0, 1), n_pts, replace = TRUE, prob = c(0.75, 0.25)),
  location = sample(c("ICU", "Ward"), n_pts, replace = TRUE),
  Age = sample(0:85, n_pts, replace = TRUE),
  infectious_syndrome = sample(c("BSI", "UTI", "RTI"), n_pts, replace = TRUE),
  date_of_admission = as.Date("2023-01-01") + sample(0:700, n_pts, replace = TRUE),
  stringsAsFactors = FALSE
)

days_to_culture <- sample(0:5, n_pts, replace = TRUE)
patients$date_of_culture <- patients$date_of_admission + days_to_culture
patients$date_of_final_outcome <- patients$date_of_culture + sample(2:25, n_pts, replace = TRUE)
patients$infection_type <- ifelse(days_to_culture >= 2, "HAI", "CAI")
patients$Age_bin <- as.character(cut(
  patients$Age,
  breaks = c(0, 18, 45, 65, Inf),
  labels = c("0-17", "18-44", "45-64", "65+"),
  right = FALSE
))

# Long format: three antibiotic results per patient
eda_data <- patients[rep(seq_len(n_pts), each = 3), ]
rownames(eda_data) <- NULL
eda_data$antibiotic_normalized <- rep(antibiotics, length.out = n_pts * 3)
eda_data$ast_value_harmonized <- sample(c("R", "S"), n_pts * 3, replace = TRUE, prob = c(0.4, 0.6))

nrow(eda_data)
#> [1] 180
```

## Theme

[`eda_theme()`](https://saketlab.github.io/anumaan/reference/eda_theme.md)
is the shared visual style behind every plot in this vignette – based on
[`ggpubr::theme_pubr()`](https://rpkgs.datanovia.com/ggpubr/reference/theme_pubr.html)
when `ggpubr` is installed, falling back to
[`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
otherwise. You do not normally call it directly; it is applied
internally by every `plot_*` function.

``` r
eda_theme(base_size = 14)
```

## Organism and Susceptibility Overview

### Top organisms

``` r
plot_top_organisms(
  eda_data,
  n = 3, mode = "faceted",
  patient_col = "patient_id", organism_col = "organism_normalized", center_col = "center_name"
) +
  ggplot2::theme(plot.margin = ggplot2::margin(10, 40, 10, 10))
```

![](eda-workflow_files/figure-html/top-organisms-1.png)

Pool across all centres with `mode = "overall"`, or restrict to one
centre with `mode = "single", center = "General Hospital"`.

### Susceptibility pattern (stacked R/S bars)

[`plot_abx_susceptibility()`](https://saketlab.github.io/anumaan/reference/plot_abx_susceptibility.md)
applies a **worst-phenotype rule**: if a patient has multiple records
for the same organism-antibiotic combination, a single R result marks
the whole episode as R, so repeat cultures from one infection episode
are not double-counted.

``` r
plot_abx_susceptibility(
  eda_data,
  n = 3, mode = "overall",
  patient_col = "patient_id", antibiotic_col = "antibiotic_normalized",
  value_col = "ast_value_harmonized", organism_col = "organism_normalized",
  center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/abx-susceptibility-1.png)

### Resistance heatmap

``` r
plot_abx_heatmap(
  eda_data,
  mode = "all",
  patient_col = "patient_id", antibiotic_col = "antibiotic_normalized",
  value_col = "ast_value_harmonized", organism_col = "organism_normalized",
  center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/abx-heatmap-1.png)

### Resistance by specimen type

``` r
plot_resistance_by_sample(
  eda_data,
  n = 3, mode = "overall",
  patient_col = "patient_id", sample_col = "specimen_normalized",
  organism_col = "organism_normalized", antibiotic_col = "antibiotic_normalized",
  value_col = "ast_value_harmonized", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/resistance-by-sample-1.png)

### Resistance by age bin and by organism

``` r
plot_resistance_by_agebin(
  eda_data,
  mode = "overall",
  patient_col = "patient_id", organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized", value_col = "ast_value_harmonized",
  agebin_col = "Age_bin", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/resistance-by-agebin-1.png)

``` r
plot_resistance_by_organism(
  eda_data,
  n = 3, mode = "overall",
  patient_col = "patient_id", organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized", value_col = "ast_value_harmonized",
  center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/resistance-by-organism-1.png)

## Outcomes

### Overall distribution

By default `merge_referred = TRUE` recodes
`"Transferred to other hospital"` to `"Referred"`, since both represent
the same clinical event across centres.

``` r
plot_outcome_distribution(
  eda_data,
  mode = "overall",
  patient_col = "patient_id", outcome_col = "final_outcome", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/outcome-distribution-1.png)

### Outcome by organism, split by resistance status

`resistance_filter = "R"` shows patients with *at least one* R result
for that organism; `"S"` shows patients with only S results. Because
these groups overlap, calling the function once for each gives the
fuller picture.

``` r
plot_outcome_by_organism(
  eda_data,
  n = 3, resistance_filter = "R", mode = "overall",
  patient_col = "patient_id", organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized", value_col = "ast_value_harmonized",
  outcome_col = "final_outcome", center_col = "center_name",
  title = "Final outcomes among resistant patients (top 3 organisms)"
)
```

![](eda-workflow_files/figure-html/outcome-by-organism-1.png)

### Death vs. discharged

``` r
plot_death_discharged(
  eda_data,
  n = 3, mode = "overall",
  patient_col = "patient_id", organism_col = "organism_normalized",
  outcome_col = "final_outcome", center_col = "center_name",
  antibiotic_col = "antibiotic_normalized", value_col = "ast_value_harmonized"
)
```

![](eda-workflow_files/figure-html/death-discharged-1.png)

### Outcome by age bin and by year

``` r
plot_outcome_by_agebin(
  eda_data,
  mode = "overall",
  patient_col = "patient_id", agebin_col = "Age_bin",
  outcome_col = "final_outcome", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/outcome-by-agebin-1.png)

``` r
plot_outcome_by_year(
  eda_data,
  mode = "overall",
  patient_col = "patient_id", outcome_col = "final_outcome",
  date_col = "date_of_final_outcome", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/outcome-by-year-1.png)

## Case-Mix and Facility Patterns

### Patients per hospital

``` r
plot_patients_by_hospital(eda_data, patient_col = "patient_id", center_col = "center_name")
```

![](eda-workflow_files/figure-html/patients-by-hospital-1.png)

### Monomicrobial vs. polymicrobial, by facility

``` r
plot_mono_poly_by_facility(
  eda_data,
  patient_col = "patient_id", poly_col = "is_polymicrobial", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/mono-poly-1.png)

### HAI vs. CAI, by facility

Computed from the gap between `admission_col` and each patient’s
earliest `culture_col`, using `hai_threshold_hours` (default 48h). It
does not read an existing infection-type column.

``` r
plot_hai_cai_by_facility(
  eda_data,
  patient_col = "patient_id", center_col = "center_name",
  admission_col = "date_of_admission", culture_col = "date_of_culture"
)
```

![](eda-workflow_files/figure-html/hai-cai-1.png)

### Location (ICU vs. Ward), by facility

``` r
plot_location_by_facility(
  eda_data,
  patient_col = "patient_id", location_col = "location", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/location-facility-1.png)

### Infection type by location

``` r
plot_infection_type_by_location(
  eda_data,
  style = "bar", patient_col = "patient_id", center_col = "center_name",
  location_col = "location", infection_col = "infection_type"
)
```

![](eda-workflow_files/figure-html/infection-type-location-1.png)

### Syndrome distribution

``` r
plot_syndrome_distribution(
  eda_data,
  mode = "overall", patient_col = "patient_id",
  syndrome_col = "infectious_syndrome", center_col = "center_name"
)
```

![](eda-workflow_files/figure-html/syndrome-distribution-1.png)

## Length of Stay and Age Distributions

### LOS ridge plot

``` r
plot_los_ridge(
  eda_data,
  mode = "all", patient_col = "patient_id", center_col = "center_name",
  admission_col = "date_of_admission", discharge_col = "date_of_final_outcome"
)
```

![](eda-workflow_files/figure-html/los-ridge-1.png)

### Age ridge plot

``` r
plot_age_ridge(eda_data, mode = "all", patient_col = "patient_id", age_col = "Age", center_col = "center_name")
```

![](eda-workflow_files/figure-html/age-ridge-1.png)

### LOS by age bin

``` r
plot_los_by_agebin(
  eda_data,
  mode = "overall", patient_col = "patient_id", agebin_col = "Age_bin",
  age_col = "Age", center_col = "center_name",
  admission_col = "date_of_admission", discharge_col = "date_of_final_outcome"
)
```

![](eda-workflow_files/figure-html/los-by-agebin-1.png)

## Summary

The EDA layer is meant to be run interactively, plot by plot, while
deciding whether the preprocessed cohort is ready for burden estimation
– not wrapped into a single automated entry point the way
[`run_preprocess()`](https://saketlab.github.io/anumaan/reference/run_preprocess.md)
is for Pipeline 1. A typical review pass:

1.  Check organism and susceptibility patterns for anything unexpected
    ([`plot_top_organisms()`](https://saketlab.github.io/anumaan/reference/plot_top_organisms.md),
    [`plot_abx_heatmap()`](https://saketlab.github.io/anumaan/reference/plot_abx_heatmap.md))
2.  Check outcome distributions and case-mix balance across centres
3.  Check LOS/age distributions for implausible values before they reach
    burden estimation
4.  Only then move on to the [DALY Burden Estimation
    vignette](https://saketlab.github.io/anumaan/articles/daly-pipeline.md)
