# Classify MDR and XDR

Classifies isolates as MDR (Multidrug Resistant) and XDR (Extensively
Drug Resistant) in a single pass using Magiorakos 2012 criteria.

## Usage

``` r
prep_classify_mdr_xdr(
  data,
  definition = "Magiorakos",
  organism_group_col = "org_group"
)
```

## Arguments

- data:

  Data frame with class-level resistance data.

- definition:

  Character. Classification criteria. Default "Magiorakos".

- organism_group_col:

  Character. Organism group column for pathogen-specific thresholds.
  Default "org_group".

## Value

Data frame with mdr, mdr_confidence, mdr_method, n_resistant_categories,
resistant_categories, xdr, xdr_confidence, and xdr_method columns added.

## Details

MDR: resistant to at least one agent in three or more antimicrobial
categories. XDR: susceptible to at most two antimicrobial categories.

Requires `class_result_event` column (run
[`prep_collapse_class_level()`](https://saketlab.github.io/anumaan/reference/prep_collapse_class_level.md)
first, then rename `class_resistance` to `class_result_event`).

## References

Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
