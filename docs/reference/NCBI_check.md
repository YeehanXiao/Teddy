# Standardize transposon annotations for TEDDY

This helper prepares a transposon annotation object for TEDDY by
converting chromosome naming to NCBI or UCSC style and ensuring that a
TE-name metadata column named \`names\` is available. If \`names\` is
absent, common TE-name columns such as \`name\`, \`repName\`,
\`repname\`, or \`TE_name\` are used to create it.

## Usage

``` r
NCBI_check(transposon, replace_name = FALSE, ncbi_style = TRUE)
```

## Arguments

- transposon:

  A [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object containing transposon annotations.

- replace_name:

  Logical. If \`TRUE\` and the object contains a \`name\` column but not
  a \`names\` column, the \`name\` column is renamed to \`names\`. If
  \`FALSE\`, a new \`names\` column is created while preserving the
  original \`name\` column. Default is \`FALSE\`.

- ncbi_style:

  Logical. If \`TRUE\` (default), chromosome names are converted to NCBI
  style (e.g., '1', '2'). If \`FALSE\`, they are converted to UCSC style
  (e.g., 'chr1', 'chr2').

## Value

A [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object with standardized chromosome names and a character metadata
column named \`names\`.
