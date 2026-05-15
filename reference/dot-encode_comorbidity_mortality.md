# Encode Comorbidity Column for Mortality Model

Standardises a free-text or numeric comorbidity column to a consistent
coding for use as a covariate in `daly_fit_mortality_or()`.

## Usage

``` r
.encode_comorbidity_mortality(
  data,
  comorbidity_col = "comorbidities",
  patient_id_col = "PatientInformation_id"
)
```

## Arguments

- data:

  Patient-level data frame.

- comorbidity_col:

  Character. Default `"comorbidities"`.

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

## Value

`data` with column `comorbidity_encoded` added. Attribute
`"comorbidity_encoding"` records the strategy used: `"numeric"`,
`"binary"`, `"ordinal"`, or `"absent"`.

## Details

Three encoding strategies are applied in order:

1.  **Numeric** (Charlson / Elixhauser index already present): used
    as-is.

2.  **Binary text** (`"present"` / `"none"`, etc.): recoded to integer 0
    / 1.

3.  **Ordinal text** (`"none"` / `"mild"` / `"moderate"` / `"severe"`):
    recoded to an ordered factor.

Missing / unknown values are set to `NA` in all cases.
