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
- Organism, specimen, antibiotic, and AST value cleaning
- Age filling and age bins
- Clinical derivations (HAI/CAI, length of stay, mortality)
- Event creation, deduplication, and readmission classification
- Contaminant detection and polymicrobial weighting
- MDR/XDR classification

Most modular functions start with `prep_`.
[`run_preprocess()`](https://saketlab.github.io/anumaan/reference/run_preprocess.md)
runs the core standardisation and derivation steps in one call,
configured through
[`amr_config()`](https://saketlab.github.io/anumaan/reference/amr_config.md).
The modular `prep_*` functions give a staged, inspectable workflow and
cover every step listed above.

``` r
library(anumaan)

# One-call pipeline
result <- run_preprocess(
  data   = your_ast_data,
  config = amr_config(hai_cutoff = 3, event_gap_days = 14)
)
clean_data <- result$data

# Or modular
data <- prep_standardize_organisms(your_ast_data,
                                   organism_col = "organism_name")
data <- prep_harmonize_ast(data,
                           antibiotic_col = "antibiotic_name",
                           ast_col        = "antibiotic_value")
data <- prep_create_event_ids(data,
                              patient_col = "patient_id",
                              date_col    = "date_of_culture")
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
  patterns, resistance heatmaps, resistance by specimen, age bin, and
  organism)
- Outcome distributions, including outcome by organism/resistance
  status, by age bin, and by year
- Case-mix and facility patterns (patients per hospital, HAI/CAI, ICU
  vs. ward, polymicrobial share, syndrome distribution)
- Length-of-stay and age distributions

The EDA functions start with `plot_` and share the
[`eda_theme()`](https://saketlab.github.io/anumaan/reference/eda_theme.md)
styling helper.

``` r
library(anumaan)

plot_top_organisms(data, n = 10, mode = "faceted",
                   patient_col = "patient_id")
plot_abx_heatmap(data, mode = "all",
                 patient_col = "patient_id")
plot_outcome_distribution(data, mode = "overall",
                          patient_col = "patient_id")
plot_los_ridge(data, mode = "all",
               admission_col = "date_of_admission",
               discharge_col = "date_of_final_outcome")
```

See the full walkthrough in the [Exploratory Data Analysis
vignette](https://saketlab.github.io/anumaan/articles/eda-workflow.md).

------------------------------------------------------------------------

## Pipeline 3 — DALY Burden Estimation

**Vignette:** [DALY Burden
Estimation](https://saketlab.github.io/anumaan/articles/daly-pipeline.md)

The DALY pipeline takes analysis-ready preprocessed data and produces
resistance-profile probabilities and GBD-style burden estimates. It runs
in six steps: resistance profiles, length of stay (LOS), relative risk
of death, population attributable fractions (PAF), years of life lost
(YLL), and years lived with disability (YLD). Burden functions start
with `daly_`.

### Resistance profiles

A resistance profile is one combination of resistant/susceptible results
across a panel of antibiotic classes. The package estimates the
probability of every profile for each pathogen through one of two
pathways.

**Pathway 1 — convex optimisation.** Works from aggregate data: marginal
resistance per class and pairwise co-resistance. It finds the profile
distribution that best matches those targets.

**Pathway 2 — Bayesian multivariate probit.** Fits a model directly to
event-level AST data, with covariates and optional random effects, and
completes the classes that were not tested.

``` r
fit <- fit_bayesian_multivariate_probit(
  event_class_data = event_level_ast,
  class_cols       = c("Carbapenems", "3GC", "Fluoroquinolones"),
  fixed_effects    = c("Age_normalised", "gender"),
  random_effects   = c("center_name"),
  pathogen         = "Klebsiella pneumoniae"
)
profiles <- compute_event_profile_probabilities(fit)
agg      <- aggregate_profiles_for_daly(profiles, hospital_col = "center_name")
```

Calibration checks for Pathway 2 are
[`validate_marginal_calibration()`](https://saketlab.github.io/anumaan/reference/validate_marginal_calibration.md),
[`validate_pairwise_calibration()`](https://saketlab.github.io/anumaan/reference/validate_pairwise_calibration.md),
[`validate_complete_profile_calibration()`](https://saketlab.github.io/anumaan/reference/validate_complete_profile_calibration.md)
and
[`mask_and_validate_ast()`](https://saketlab.github.io/anumaan/reference/mask_and_validate_ast.md).

#### Bayesian multivariate probit backend selection

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

### Length of stay

LOS is the time from admission to discharge among discharged admissions.
[`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md)
fits a Gamma model of LOS for each pathogen and antibiotic class,
adjusted for age, sex, infection type, ICU and hospital. It returns how
much longer resistant infections stay than susceptible ones (relative
LOS), and for each hospital the adjusted LOS of resistant admissions
(`LOSR`) and of the same admissions had they been susceptible (`LOSS`),
with bootstrap confidence intervals.
[`daly_fit_los_rr_distribution()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr_distribution.md)
and
[`daly_fit_los_distn_by_profile()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_distn_by_profile.md)
are distribution-based alternatives.

### Relative risk of death

[`daly_fit_mortality_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_mortality_rr.md)
fits a logistic model of death on resistance for each pathogen and
antibiotic class, adjusted for age, sex and infection type. It returns
an odds ratio (`OR_death`) and a relative risk (`RR_death`).

### Population attributable fractions

The per-class estimates are mapped onto profiles: each profile takes the
largest relative risk among its resistant classes
([`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)).
The PAF is then the share of LOS, or of deaths, that would not have
occurred had every infection been susceptible
([`daly_calc_paf_los()`](https://saketlab.github.io/anumaan/reference/daly_calc_paf_los.md),
[`daly_calc_paf_mortality()`](https://saketlab.github.io/anumaan/reference/daly_calc_paf_mortality.md)).

### Years of life lost

[`daly_calc_yll_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yll_associated.md)
sums the remaining life expectancy, by age bin and sex, of patients who
died with the infection.
[`daly_calc_yll_attributable()`](https://saketlab.github.io/anumaan/reference/daly_calc_yll_attributable.md)
multiplies that by the mortality PAF to give the share due to
resistance.

### Years lived with disability

[`daly_calc_yld_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_associated.md)
combines the number of non-fatal cases, the pathogen’s share of them,
the profile probabilities, LOS and a disability weight.
[`daly_calc_yld_attributable()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_attributable.md)
uses the excess LOS of resistant over susceptible infections instead.

See the full walkthrough in the [DALY Burden Estimation
vignette](https://saketlab.github.io/anumaan/articles/daly-pipeline.md).

------------------------------------------------------------------------

## Installation

``` r
# From GitHub
remotes::install_github("saketlab/anumaan")

# Build vignettes locally when you want the full article set
remotes::install_github("saketlab/anumaan", build_vignettes = TRUE)
```

### Optional dependencies

| Feature | Packages |
|----|----|
| Bayesian multivariate probit fitting ([`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)) | `cmdstanr`, `posterior` |
| Joint Bernoulli resistance-profile fitting ([`profile_fit_joint_bernoulli()`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)) | `cmdstanr`, `posterior` |
| Sampler diagnostic plots ([`profile_check_sampler()`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md) and others) | `bayesplot` |
| Mixed-effects LOS modelling | `lme4`, `glmmTMB` |
| Spatial analysis | `sf`, `spdep`, `leaflet` |
| Python ICD-10 embedding | `reticulate` + `alethia` |
| pkgdown site build | `pkgdown` |
| Enhanced EDA theme | `ggpubr` (plots fall back to [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html) when absent) |

------------------------------------------------------------------------

## Data Format

All three pipelines expect **long format** input — one row per isolate ×
antibiotic combination:

| Column | Role |
|----|----|
| `patient_id` | Unique patient identifier |
| `isolate_id` | Unique isolate identifier (counting unit for resistance stats) |
| `organism_name` | Organism name (raw) |
| `antibiotic_name` | Antibiotic name (raw) |
| `antibiotic_value` | Susceptibility result: R, I, or S |
| `date_of_culture` | Culture collection date |
| `specimen_type` | Specimen source |

Optional but used when present: `date_of_admission`,
`date_of_final_outcome`, `DOB`, `Age`, `gender`, `final_outcome`,
`infection_type`, `state`.

[`prep_standardize_specimens()`](https://saketlab.github.io/anumaan/reference/prep_standardize_specimens.md)
returns `sterile_classification` collapsed to `Sterile`, `Non-Sterile`,
or `Others/Ambiguous`. The specimen reference also includes stewardship
specimen labels such as `Brain abscess`, `Instrument`, `Lung aspirate`,
`Lymph node`, `Hair`, and `Superficial Biopsy`.

------------------------------------------------------------------------

## Package Reference

Full function documentation is at the [Reference
page](https://saketlab.github.io/anumaan/reference/index.html).

``` r
# Search all functions
help(package = "anumaan")

# Common entry points
?run_preprocess
?amr_config
?fit_bayesian_multivariate_probit
?aggregate_profiles_for_daly
```

Report bugs at <https://github.com/saketlab/anumaan/issues>.
