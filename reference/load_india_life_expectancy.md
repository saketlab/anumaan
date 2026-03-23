# Load and Parse India Life Expectancy Lookup Table

Reads the statewise life expectancy Excel file and returns a tidy data
frame with columns `sex`, `age_bin`, and `life_expectancy` for the India
column only. Three sex sections are parsed: `"Combined"`, `"Male"`, and
`"Female"`.

## Usage

``` r
load_india_life_expectancy(le_path)
```

## Arguments

- le_path:

  Character. Full path to the life expectancy xlsx file.

## Value

Data frame with columns `sex` (`"Combined"`, `"Male"`, `"Female"`),
`age_bin`, `life_expectancy`.
