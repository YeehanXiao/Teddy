# Extract the result from the differentially expressed TE-chimeric exon test

Extract the result from the differentially expressed TE-chimeric exon
test

## Usage

``` r
extractTest(object, filter = TRUE)
```

## Arguments

- object:

  An object of DE TE-chimeric exon test, out from
  **ChimericDrivenTest**.

- filter:

  Logical. If `TRUE` (default), filter out features with `NA` adjusted
  p-values and features flagged as structural zeros (`allZero`).
