# Build QP Constraint Matrix and Target Vector

Constructs the constraint matrix **M** and target vector **v** used in
the simplex-constrained weighted least-squares QP:

## Usage

``` r
build_constraint_matrix(profiles_enum, r_marg, co_mat = NULL)
```

## Arguments

- profiles_enum:

  Tibble from
  [`enumerate_binary_profiles()`](https://saketlab.github.io/anumaan/reference/enumerate_binary_profiles.md).
  The `profile_delta` column is ignored; only the binary class columns
  are used.

- r_marg:

  Named numeric vector. Marginal resistance rates, one per class. Names
  must match column names in `profiles_enum`.

- co_mat:

  Numeric matrix or `NULL`. Square symmetric pairwise co-resistance
  prevalence matrix with row/column names matching `names(r_marg)`. `NA`
  cells trigger independence fallback. Default `NULL`.

## Value

Named list:

- `M`:

  Numeric matrix \\(n + n(n-1)/2) \times 2^n\\.

- `v`:

  Numeric vector of constraint targets.

- `constraint_names`:

  Character vector labelling each row of M.

- `fallback_pairs`:

  Character vector of class pairs that used the independence fallback.

- `capped_pairs`:

  Named numeric vector of pairs whose pairwise value was capped, showing
  original and capped values.

## Details

\$\$\min_p \\Mp - v\\^2_W \quad \text{s.t.} \quad p \ge 0,\\ \sum p =
1\$\$

Marginal rows: \\M\_{d,\delta} = 1\\ iff class \\d\\ is resistant in
profile \\\delta\\. Pairwise rows: \\M\_{d_1 d_2,\delta} = 1\\ iff both
classes are resistant. When a pairwise value is unavailable (too few
co-tested isolates), the product of marginals (independence assumption)
is substituted. Pairwise values that exceed \\\min(P(A), P(B))\\ are
capped.
