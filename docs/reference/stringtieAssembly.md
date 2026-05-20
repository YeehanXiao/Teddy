# R wrapper to Run stringtie

Run stringtie for transcriptome assembly

## Usage

``` r
stringtieAssembly(
  bam,
  reference = NULL,
  outfile,
  params = "",
  longRead = FALSE
)
```

## Arguments

- bam:

  Character, BAM file (sorted by genomic coordinates)

- reference:

  Character, reference annotation GTF (optional)

- outfile:

  Character, output GTF file name

- params:

  Character, extra parameters (default "")

- longRead:

  Logical, whether the data is long-read (PacBio/ONT). Default FALSE
  (short-read).
