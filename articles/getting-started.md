# Getting Started with anumaan

## What is anumaan?

**anumaan** is an R package for preprocessing antimicrobial resistance
(AMR) surveillance data and estimating disease burden using GBD
methodology. It is structured around two sequential, independent
pipelines that can be used together or separately.

------------------------------------------------------------------------

## Pipeline 1 — Preprocessing

**Vignette:** [Preprocessing
Workflow](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md)

The preprocessing pipeline standardises raw hospital or surveillance
data into an analysis-ready format. It covers:

- Column name standardisation and schema alignment
- Date parsing and chronological validation
- Organism, specimen, antibiotic, and AST value cleaning
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

## Pipeline 2 — DALY Burden Estimation

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

------------------------------------------------------------------------

## Installation

``` r
# From GitHub
remotes::install_github("saketlab/anumaan")
```

### Optional dependencies

| Feature                        | Packages                 |
|--------------------------------|--------------------------|
| Convex QP solver (recommended) | `osqp`, `Matrix`         |
| Convex QP solver (fallback)    | `quadprog`               |
| Mixed-effects LOS modelling    | `lme4`, `glmmTMB`        |
| Spatial analysis               | `sf`, `spdep`, `leaflet` |
| Python ICD-10 embedding        | `reticulate` + `alethia` |

------------------------------------------------------------------------

## Data Format

Both pipelines expect **long format** input — one row per isolate ×
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
