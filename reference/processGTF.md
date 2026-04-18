# Process GTF and SummarizedExperiment for TE overlap analysis

This function processes a GTF file or a SummarizedExperiment object to
identify overlaps with transposable elements (TE). It can handle both
external (local) GTF files and in-memory SummarizedExperiment objects.
The function extracts exon and transcript information, finds overlaps
with TE annotations, and filters out transcripts with zero expression.

## Usage

``` r
processGTF(te, gtf_path = NULL, combineSE = NULL, minoverlap = 0, threads = 9)
```

## Arguments

- te:

  GRanges object. Transposable element annotations, required for overlap
  analysis. The \`te\` object must contain at least three metadata
  columns: \`names\`, \`family\`, and \`class\`.

- combineSE:

  SummarizedExperiment. Pre-processed RNA-seq quantification object
  (optional, required if gtfFile is NULL).

- minoverlap:

  Integer. Minimum required overlap (in base pairs) between GTF exons
  and TE annotations.

- threads:

  Integer. Number of threads to use for parallel processing. Default is
  4.

- gtfFile:

  Character. Path to the GTF file (optional, required if combineSE is
  NULL).

## Value

A GRanges object with transcripts overlapping TEs, annotated with TE
information.
