# Consolidation of information on transcript assembly of multiple samples

Transcript-quantification is prerequisite for many downstream
investigations. Several metrics have been proposed for measuring
abundance in transcript level based on RNA-seq data.

## Usage

``` r
stringtieCombine(
  reference = NULL,
  bamFiles = NULL,
  gtfFiles = NULL,
  params = "",
  longRead = FALSE,
  cores = 1
)
```

## Arguments

- reference:

  Compared with a reference annotation file

- bamFiles:

  Character vector. BAM files (sorted by coordinate).

- gtfFiles:

  Character vector. Output GTF files (will be overwritten after
  re-quantification).

- params:

  Other parameters for StringTie

- longRead:

  Logical scalar or vector: TRUE for ONT/PacBio

- cores:

  Number of cores to use in parallel (default: 1)
