# Getting Started with anumaan

## What is anumaan?

**anumaan** is an R package for preprocessing antimicrobial resistance
(AMR) surveillance data and estimating disease burden using GBD
methodology. It is structured around three sequential pipelines that can
be used together or separately.

------------------------------------------------------------------------

## Pipeline 1 — Preprocessing

**Vignette:** [Preprocessing
Workflow](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md)

The preprocessing pipeline standardises raw hospital or surveillance
data into an analysis-ready format. It covers:

- Column name standardisation and schema alignment
- Date parsing and chronological validation
- Organism, specimen, antibiotic, and AST value cleaning, including the
  0.1.0.9021 specimen reference additions and collapsed sterile-site
  classes
- Demographic derivations (age, HAI/CAI, LOS)
- Event creation, deduplication, and readmission classification
- Contaminant detection and polymicrobial weighting
- Attrition tracking and analysis-readiness filtering

All functions start with `prep_`. The pipeline entry point is
[`run_preprocess()`](https://saketlab.github.io/anumaan/reference/run_preprocess.md)
for automated runs or the modular `prep_*` functions for staged,
inspectable workflows.

``` r
library(anumaan)

# Full automated pipeline
result <- run_preprocess(
  data   = your_ast_data,
  config = amr_config(hai_cutoff = 3, event_gap_days = 14)
)

# Or modular
data <- prep_standardize_organisms(data, organism_col = "organism_name")
data <- prep_harmonize_ast(data, ast_col = "antibiotic_value")
data <- prep_create_event_ids(data, patient_col = "patient_id",
                               date_col = "date_of_culture")
```

See the full step-by-step walkthrough in the [Preprocessing Workflow
vignette](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md).

------------------------------------------------------------------------

## Pipeline 2 — Exploratory Data Analysis

**Vignette:** [Exploratory Data
Analysis](https://saketlab.github.io/anumaan/articles/eda-workflow.md)

The EDA pipeline visualises preprocessed data before it moves on to
burden estimation. It does not transform data – it is the review step a
stewardship team runs on Pipeline 1’s output to catch data-quality
issues and characterise the cohort. It covers:

- Organism and antibiotic-susceptibility overviews (top organisms, R/S
  patterns, resistance heatmaps)
- Outcome distributions, including outcome by organism/resistance
  status, by age bin, and by year
- Case-mix and facility patterns (HAI/CAI, ICU vs. ward, polymicrobial
  share, syndrome distribution)
- Length-of-stay and age distributions

All functions start with `plot_`, plus the shared
[`eda_theme()`](https://saketlab.github.io/anumaan/reference/eda_theme.md)
styling helper – deliberately distinct from the `prep_` (Pipeline 1) and
`daly_` (Pipeline 3) naming.

``` r
library(anumaan)

plot_top_organisms(data, n = 10, mode = "faceted")
plot_abx_heatmap(data, mode = "all")
plot_outcome_distribution(data, mode = "overall")
plot_los_ridge(data, mode = "all")
```

See the full walkthrough in the [Exploratory Data Analysis
vignette](https://saketlab.github.io/anumaan/articles/eda-workflow.md).

------------------------------------------------------------------------

## Pipeline 3 — DALY Burden Estimation

**Vignette:** [DALY Burden
Estimation](https://saketlab.github.io/anumaan/articles/daly-pipeline.md)

The DALY pipeline takes analysis-ready preprocessed data and produces
resistance-profile probability distributions and GBD-style burden
estimates. It covers:

- Resistance profile estimation via convex optimisation (Pathway 1) —
  works from either facility line-list data or pre-computed aggregate
  marginals (GBD ST-GPR, GLASS, national surveillance networks)
- Marginal and pairwise co-resistance computation with Pearson
  back-calculation
- Profile probability estimation solving a simplex-constrained QP
- Bootstrap uncertainty intervals for profile probabilities
- Years of Life Lost (YLL) — associated and attributable
- Years Lived with Disability (YLD) — associated and attributable
- Total DALY burden per pathogen, hospital, and organism group

All estimation functions start with `daly_`. Profile-specific functions
are `compute_*`, `estimate_*`, `enumerate_*`, `build_*`, `validate_*`,
`check_*`, and `bootstrap_*`.

``` r
# From aggregate marginals (GBD / GLASS / national surveillance)
profiles <- estimate_profiles_convex(
  marginals = your_marginals_table,
  panel_map = list(
    "Klebsiella pneumoniae" = c("Carbapenems", "3GC", "Fluoroquinolones")
  )
)

# Assign LOS relative risk and compute fatal resistance prevalence
profiles_rr <- daly_assign_rr_to_profiles(profiles, rr_table)
R_k         <- daly_calc_resistance_prevalence_fatal(profiles_rr)
```

See the full walkthrough in the [DALY Burden Estimation
vignette](https://saketlab.github.io/anumaan/articles/daly-pipeline.md).

### Bayesian multivariate probit backend selection

For Pathway 2 Bayesian resistance-profile estimation,
[`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)
now accepts a `compute` argument so the same statistical model can be
run on either a standard CPU backend or an OpenCL-enabled backend.

``` r
fit <- fit_bayesian_multivariate_probit(
  event_class_data = event_level_ast,
  class_cols = c("3GC", "Carbapenems", "Fluoroquinolones"),
  fixed_effects = c("Age_normalised", "gender", "location"),
  random_effects = list(
    list(name = "hospital", group_col = "center_name", terms = "intercept")
  ),
  pathogen = "klebsiella pneumoniae",
  residual_structure = "correlated",
  compute = list(
    backend = "cpu"
  )
)
```

``` r
fit_opencl <- fit_bayesian_multivariate_probit(
  event_class_data = event_level_ast,
  class_cols = c("3GC", "Carbapenems", "Fluoroquinolones"),
  fixed_effects = c("Age_normalised", "gender", "location"),
  random_effects = list(
    list(name = "hospital", group_col = "center_name", terms = "intercept")
  ),
  pathogen = "klebsiella pneumoniae",
  residual_structure = "correlated",
  compute = list(
    backend = "opencl",
    opencl_platform_id = 0L,
    opencl_device_id = 0L
  )
)
```

This backend choice affects compilation and sampling only. Priors,
likelihood, diagnostics, posterior predictive checks, resistance
profiles, and downstream DALY quantities remain unchanged.

------------------------------------------------------------------------

## Installation

``` r
# From GitHub
remotes::install_github("saketlab/anumaan")

# Build vignettes locally when you want the full article set
remotes::install_github("saketlab/anumaan", build_vignettes = TRUE)
```

### Optional dependencies

| Feature                        | Packages                                                                                                                |
|--------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| Convex QP solver (recommended) | `osqp`, `Matrix`                                                                                                        |
| Convex QP solver (fallback)    | `quadprog`                                                                                                              |
| Mixed-effects LOS modelling    | `lme4`, `glmmTMB`                                                                                                       |
| Spatial analysis               | `sf`, `spdep`, `leaflet`                                                                                                |
| Python ICD-10 embedding        | `reticulate` + `alethia`                                                                                                |
| pkgdown site build             | `pkgdown`                                                                                                               |
| Enhanced EDA theme             | `ggpubr` (plots fall back to [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html) when absent) |

------------------------------------------------------------------------

## Data Format

All three pipelines expect **long format** input — one row per isolate ×
antibiotic combination:

| Column             | Role                                                           |
|--------------------|----------------------------------------------------------------|
| `patient_id`       | Unique patient identifier                                      |
| `isolate_id`       | Unique isolate identifier (counting unit for resistance stats) |
| `organism_name`    | Organism name (raw)                                            |
| `antibiotic_name`  | Antibiotic name (raw)                                          |
| `antibiotic_value` | Susceptibility result: R, I, or S                              |
| `date_of_culture`  | Culture collection date                                        |
| `specimen_type`    | Specimen source                                                |

Optional but used when present: `date_of_admission`,
`date_of_final_outcome`, `DOB`, `Age`, `gender`, `final_outcome`,
`infection_type`, `state`.

[`prep_standardize_specimens()`](https://saketlab.github.io/anumaan/reference/prep_standardize_specimens.md)
now returns `sterile_classification` collapsed to `Sterile`,
`Non-Sterile`, or `Others/Ambiguous`. The 0.1.0.9021 reference also adds
stewardship specimen labels such as `Brain abscess`, `Instrument`,
`Lung aspirate`, `Lymph node`, `Hair`, and `Superficial Biopsy`.

------------------------------------------------------------------------

## Package Reference

Full function documentation is at the [Reference
page](https://saketlab.github.io/anumaan/reference/index.md).

``` r
# Search all functions
help(package = "anumaan")

# Common entry points
?run_preprocess
?amr_config
?estimate_profiles_convex
?compute_resistance_profiles
```

Report bugs at <https://github.com/saketlab/anumaan/issues>.
