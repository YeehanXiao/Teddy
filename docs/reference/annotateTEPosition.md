# Annotate TE positions within transcript structures

Annotate TE positions within transcript structures

## Usage

``` r
annotateTEPosition(
  chi_GTF,
  te,
  gtf = NULL,
  combineSE = NULL,
  minoverlap = 5L,
  boundary_tolerance = NULL
)
```

## Arguments

- chi_GTF:

  GRanges object returned by processGTF().

- te:

  GRanges object containing TE annotations.

- gtf:

  GRanges object or path to the complete transcript annotation.

- combineSE:

  SummarizedExperiment containing the complete GTF in metadata.

- minoverlap:

  Minimum overlap between exon and TE.

- boundary_tolerance:

  Optional distance from an exon boundary used for positional
  refinement. NULL preserves the original TEDDY classification only.

## Value

A GRanges object containing exact TE-exon pairs and positional
annotations.
