# Collapse drug-level data to antibiotic-class binary wide matrix

Standardises antibiotic values (Intermediate -\> R/S per GBD rules),
collapses to class level (class = 1 if ANY drug in class = R), and
pivots wide: one row per patient, one column per antibiotic class (0/1).
Column names are sanitised with make.names(); the original-to-safe name
mapping is stored as the "class_name_map" attribute.

Standardises antibiotic values (Intermediate -\> R/S per GBD rules),
collapses to class level (class = 1 if ANY drug in class = R), and
pivots wide: one row per patient, one column per antibiotic class (0/1).
Column names are sanitised with make.names(); the original-to-safe name
mapping is stored as the "class_name_map" attribute.

## Usage

``` r
.build_class_resistance_wide(
  data,
  patient_id_col = "PatientInformation_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  untested_fill = 0L
)

.build_class_resistance_wide(
  data,
  patient_id_col = "PatientInformation_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  untested_fill = 0L
)
```
