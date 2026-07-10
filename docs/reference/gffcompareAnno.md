# R wrapper to run gffcompare

The function to compare one or more GTF/GFF files (the query files) with
a reference annotation using gffcompare. The final annotated GTF is
renamed to the user-specified \`outfile\`. Tracking, loci and summary
files are retained, while temporary \`.tmap\` and \`.refmap\` files are
removed by default.

## Usage

``` r
gffcompareAnno(
  reference,
  gtffile,
  outfile,
  params = "",
  cleanup = TRUE,
  overwrite = FALSE
)
```

## Arguments

- reference:

  Use a reference annotation file to guide compare assembly process.

- gtffile:

  GTF/GFF file(s) with assembled transcripts.

- outfile:

  The name of the final output annotated GTF.

- params:

  Other parameters passed to gffcompare.

- cleanup:

  Logical. Whether to remove \`.tmap\` and \`.refmap\` files generated
  by gffcompare. Default is TRUE.

- overwrite:

  Logical. Whether to overwrite existing output files. Default is FALSE.
