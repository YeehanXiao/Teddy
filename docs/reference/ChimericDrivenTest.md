# Testing for differential TE-chimeric exon usage

Detect to what extent TE-chimeric exon affect the expression of the
corresponding transcript

## Usage

``` r
ChimericDrivenTest(
  SEobject,
  condition = NULL,
  annotation = NULL,
  maxit = 100,
  niter = 10,
  quiet = FALSE,
  warning = FALSE,
  filterAllZero = FALSE
)
```

## Arguments

- SEobject:

  An object of RangedSummarizedExperiment class, out from **countAnno**.

- condition:

  Vector, indicating the experimental condition of the samples.

- annotation:

  The genome annotation object, out from **prepareAnno**.

- maxit:

  Control parameter: maximum number of iterations to allow for
  convergence when calculating dispersion.

- niter:

  Number of times to iterate between estimation of means and estimation
  of dispersion.

- quiet:

  Whether to suppress messages at each step.

- warning:

  Whether to print warnings at each step.

- filterAllZero:

  Logical, whether to remove all-zero rows before fitting the DESeq2
  model.
