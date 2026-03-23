# Find Package Data File

Helper function to locate data files in inst/extdata. Works both when
the package is installed and during development.

## Usage

``` r
find_extdata_file(filename)
```

## Arguments

- filename:

  Character. Name of the file to find (e.g., "organisms.csv").

## Value

Character. Full path to the file, or empty string if not found.
