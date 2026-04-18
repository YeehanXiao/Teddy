# Extract GTF information

This function filters and extracts GTF information from a
\`SummarizedExperiment\` object based on the specified element type
(exon, transcript, or both) and an FPKM threshold.

## Usage

``` r
extractGTF(combineSE, filter = 1, type = c("exon", "transcript", "both"))
```

## Arguments

- combineSE:

  A `SummarizedExperiment` object, typically the output from the
  `combineSE` function. It should contain GTF metadata in its
  `@metadata$gtf` slot and FPKM values in one of its assays.

- filter:

  A numeric value specifying the minimum FPKM threshold for transcripts
  to be included in the output. Default is 1.

- type:

  A character string specifying the type of genetic elements to extract.
  Valid options are "exon", "transcript", or "both". Default is "exon".

## Value

A filtered subset of the GTF metadata from the \`combineSE\` object,
including only the specified types of genetic elements that meet the
FPKM threshold.
