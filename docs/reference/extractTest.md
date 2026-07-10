# Extract results from the differential TE-bin usage test

Extract results from the differential TE-bin usage test

## Usage

``` r
extractTest(object, filter = TRUE)
```

## Arguments

- object:

  An object from **ChimericDrivenTest**.

- filter:

  Logical. If `TRUE` (default), filter out features with `NA` adjusted
  p-values and features flagged as structural zeros (`allZero`).
