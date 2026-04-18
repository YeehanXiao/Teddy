# Calculating the fold change

TE-chimeric exon usage fold changes are calculated based on the
coefficients of the GLM fit

## Usage

``` r
calculateFoldchange(object, genes, crossVar = "condition")
```

## Arguments

- object:

  An object of test

- genes:

  Specify the genes for the calculation of foldchanges

- corssVar:

  Default as "condition"
