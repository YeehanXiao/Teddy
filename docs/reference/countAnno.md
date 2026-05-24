# Counting reads on exon bins

Count reads overlapping exon bins defined in the flattened GTF/GFF
object.

## Usage

``` r
countAnno(
  annotation,
  bamfile,
  isPairedEnd = TRUE,
  strandSpecific = 0,
  allowMultiOverlap = TRUE,
  isLongRead = FALSE,
  nthreads = 1,
  ...
)
```

## Arguments

- annotation:

  A GRanges object, typically from prepareAnno().

- bamfile:

  Character vector of BAM file(s).

- isPairedEnd:

  Logical scalar or vector, whether the library is paired-end.

- strandSpecific:

  Integer: 0 (unstranded), 1 (forward), or 2 (reverse).

- allowMultiOverlap:

  Logical, whether to allow reads to overlap multiple features.

- isLongRead:

  Logical, whether this is long-read data (e.g. ONT, PacBio).

- nthreads:

  Number of BAM files processed in parallel.

- ...:

  Additional arguments passed to Rsubread::featureCounts().

## Value

A RangedSummarizedExperiment object containing a counts assay.
