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
  warning = FALSE
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

  control parameter: maximum number of iterations to allow for
  convergence when calculating dispersion.

- niter:

  whether to print messages at each step

- quiet:

  number of times to iterate between estimation of means and estimation
  of dispersion

- warning:

  whether to print warning at each step
