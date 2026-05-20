# Build an isoform support matrix across multiple GTF files

Compare multiple query GTF files against a reference transcriptome and
return a presence/absence matrix based on gffcompare tracking files.

## Usage

``` r
buildIsoformSupport(
  reference,
  gtffiles,
  conditions,
  out_dir = ".",
  out_prefix = "cmp",
  params = "",
  FPKM = NULL,
  threshold = 0.1,
  cores = 4,
  fill_missing_annot = FALSE
)
```

## Arguments

- reference:

  Character. Path to the reference GTF.

- gtffiles:

  Character vector. Paths to query GTF files.

- conditions:

  Character vector. Sample, platform, or depth names corresponding to
  \`gtffiles\`.

- out_dir:

  Character. Output directory for gffcompare results.

- out_prefix:

  Character. Prefix for output files.

- params:

  Additional parameters passed to gffcompare.

- FPKM:

  Optional numeric matrix. Rows correspond to isoforms and columns
  correspond to \`conditions\`.

- threshold:

  Numeric. FPKM cutoff used when refining presence.

- cores:

  Integer. Number of cores for parallel execution.

- fill_missing_annot:

  Logical. Whether to fill missing annotation using
  \`lookup_annot_tx()\`. Default is FALSE.

## Value

A tibble with one row per reference transcript and one logical column
per condition.
