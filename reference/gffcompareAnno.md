# R wrapper to Run gffcompare

The function to compare and merge accuracy of one or more GFF files (the
“query” files), when compared with a reference annotation (also provided
as GFF).

## Usage

``` r
gffcompareAnno(reference, gtffile, outfile, params = "")
```

## Arguments

- reference:

  Use a reference annotation file to guide compare assembly process.

- gtffile:

  GTF files with gffcompare annotation transcripts.

- outfile:

  The name of the output annotated merged GTF.

- params:

  Other parameters
