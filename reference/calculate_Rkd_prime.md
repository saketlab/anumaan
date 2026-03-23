# Calculate non-fatal prevalence of resistance (R'\_{k,d})

Computes isolate-wise prevalence of resistance for each pathogen-drug
combination, with optional facility stratification. Antibiotics are
collapsed to classes, retaining the maximum resistance prevalence within
each class.

## Usage

``` r
calculate_Rkd_prime(
  ast_data,
  isolate_col = "isolate_id",
  pathogen_col = "pathogen",
  antibiotic_col = "antibiotic",
  ast_result_col = "ast_result",
  drug_class_col = NULL,
  antibiotic_class_map = NULL,
  facility_col = NULL,
  facility_name = NULL,
  pathogen_name = NULL
)
```

## Arguments

- ast_data:

  Data frame containing AST results.

- isolate_col:

  Character. Unique isolate ID column.

- pathogen_col:

  Character. Pathogen column (K).

- antibiotic_col:

  Character. Antibiotic name column (d).

- ast_result_col:

  Character. AST result column ("R", "S", "I").

- drug_class_col:

  Character or NULL. Pre-existing drug class column in `ast_data`. Used
  instead of `antibiotic_class_map` when supplied.

- antibiotic_class_map:

  Data frame with columns: antibiotic, drug_class.

- facility_col:

  Character or NULL. Facility column if present.

- facility_name:

  Character or NULL. If provided, filters data to the specified facility
  before computation.

- pathogen_name:

  Character vector or NULL. If provided, filters data to the specified
  pathogen(s).

## Value

Data frame with columns: pathogen, drug_class, R_kd_prime, N_tested,
N_resistant (+ facility if provided)
