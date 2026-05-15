# Standardize Final Outcome Column

Normalizes outcome variants to a two-level controlled vocabulary:
"Survived", "Died", or NA.

## Usage

``` r
prep_standardize_final_outcome(data, col = "final_outcome")
```

## Arguments

- data:

  Data frame.

- col:

  Character. Outcome column. Default "final_outcome".

## Value

Data frame with `outcome_std` column added.

## Details

- Survived:

  Survived / Alive / Discharge / Discharged / Recovered

- Died:

  Died / Death / Expired / Deceased / Dead

- NA:

  Unknown / Missing / Absconded / LAMA / DAMA
