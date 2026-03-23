# Extract LOS Vectors by Resistance Status

For a given organism (and optionally center), computes per-patient
resistance status and returns separate LOS vectors for resistant and
susceptible exposures.

## Usage

``` r
get_los_by_resistance(
  abx_data,
  los_data,
  organism,
  center = NULL,
  patient_id_col = "PatientInformation_id",
  organism_col = "organism_clean",
  center_col = "center_name",
  antibiotic_col = "antibiotic_name",
  sus_col = "sus"
)
```

## Arguments

- abx_data:

  Data frame of antibiotic susceptibility results with columns for
  patient ID, organism, antibiotic name, and susceptibility value
  (already cleaned to "R"/"S").

- los_data:

  Data frame from
  [`prepare_los_data()`](https://saketlab.github.io/anumaan/reference/prepare_los_data.md).

- organism:

  Character. Organism name (lowercase, trimmed) to filter on.

- center:

  Optional character. Center name to filter on. Default `NULL` (all
  centers).

- patient_id_col:

  Character. Patient ID column. Default `"PatientInformation_id"`.

- organism_col:

  Character. Organism column. Default `"organism_clean"`.

- center_col:

  Character. Center column. Default `"center_name"`.

- antibiotic_col:

  Character. Antibiotic name column. Default `"antibiotic_name"`.

- sus_col:

  Character. Susceptibility column (values "R"/"S"). Default `"sus"`.

## Value

List with elements `R` (LOS vector for resistant exposures) and `S` (LOS
vector for susceptible exposures).
