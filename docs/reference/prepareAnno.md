# Preparing the genome annotation object

Flatten exon appearing multiple times among different transcripts in GTF
file

Flatten exon appearing multiple times among different transcripts in GTF
file

## Usage

``` r
prepareAnno(
  gtffile,
  singleGens = TRUE,
  transposon = NULL,
  minoverlap = 10,
  cores = 4
)
```

## Arguments

- gtffile:

  GTF file.

- singleGens:

  Whether to allocate the exon overlapping with two genes to a single
  gene. Default is TRUE.

- transposon:

  A GRanges object with transposon data

- minoverlap:

  Minimum overlap for
  [`findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).
  Default is 10.

- cores:

  Number of cores used for parallel processing. Default is 1.
