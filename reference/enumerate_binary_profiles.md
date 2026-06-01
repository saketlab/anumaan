# Enumerate All Binary Resistance Profiles for a Set of Antibiotic Classes

Generates all \\2^n\\ binary resistance profiles for `n` antibiotic
classes using vectorised bit manipulation. Each row is one profile; each
class column contains 1 (resistant) or 0 (susceptible).

## Usage

``` r
enumerate_binary_profiles(classes)
```

## Arguments

- classes:

  Character vector. Ordered antibiotic class names. The order determines
  which bit position each class occupies and must be consistent across
  all calls within a single analysis.

## Value

Tibble with `2^n` rows and `n + 1` columns: `profile_delta` (character
label, e.g. `"RSS"`) and one integer column per class (0 / 1).
