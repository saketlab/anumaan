# Map Diagnosis Text to ICD Candidates

Maps free-text diagnosis strings to ICD-10 candidate descriptions using
one of three methods:

## Usage

``` r
prep_map_diagnosis_to_icd(
  data,
  text_col = "diagnosis_text",
  reference = NULL,
  method = c("exact", "fuzzy", "python_embedding"),
  icd_desc_col = "description3",
  icd_code_col = "icd_code_who_eq",
  top_k = 5L,
  threshold = 0,
  model = "FremyCompany/BioLORD-2023",
  id_col = NULL
)
```

## Arguments

- data:

  Data frame.

- text_col:

  Character. Column containing prepared diagnosis text (output of
  [`prep_diagnosis_text()`](https://saketlab.github.io/anumaan/reference/prep_diagnosis_text.md)).
  Default `"diagnosis_text"`.

- reference:

  Data frame or NULL. ICD-10 reference table. If NULL, loads from
  `inst/extdata/icd10_who.csv`.

- method:

  Character. One of `"exact"`, `"fuzzy"`, `"python_embedding"`. Default
  `"exact"`.

- icd_desc_col:

  Character. Column in `reference` to match against. Default
  `"description3"` (most concise ICD labels).

- icd_code_col:

  Character. Column in `reference` holding ICD codes. Default
  `"icd_code_who_eq"`.

- top_k:

  Integer. Maximum ICD candidates to return per input string. Default 5.
  Ignored for `"exact"` (returns all exact matches).

- threshold:

  Numeric. Minimum similarity score (0-1) to retain a candidate. Default
  0.0 (keep all). For `"fuzzy"`: similarity is
  `1 - normalised_distance`.

- model:

  Character. Sentence-transformers model name. Used only for
  `"python_embedding"`. Default `"FremyCompany/BioLORD-2023"`.

- id_col:

  Character or NULL. Identifier column to carry through into the output.
  Default NULL (output contains only match columns).

## Value

Long data frame with columns: `diagnosis_text`, `icd_prediction`,
`icd_code`, `icd_score`, `icd_rank`, `icd_method`. If `id_col` is
supplied, it is included as the first column.

## Details

- `"exact"`:

  Case-insensitive exact string match against ICD descriptions. Fast,
  high precision, low recall.

- `"fuzzy"`:

  String distance matching via stringdist. Handles typos and minor
  variations. Requires stringdist.

- `"python_embedding"`:

  Semantic embedding similarity using the Python `alethia` package via
  reticulate. Highest recall. Requires Python, `alethia`, and a
  sentence-transformers model.
