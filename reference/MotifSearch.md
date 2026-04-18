# Search the motif binding sites through TE-chimeric transcripts

Search the motif binding sites through TE-chimeric transcripts

## Usage

``` r
MotifSearch(
  object,
  te,
  pwm,
  minoverlap = 15,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  min.score = 0.9,
  filter = 0
)
```

## Arguments

- object:

  A Granges object of target TE-chimeric transcripts.

- te:

  A Granges object of TE location annotations, typically derived from
  annotated TE reference.

- pwm:

  The pwm matrix for the target motif.

- minoverlap:

  A non-negative integer specifying the minimum number for ranges to be
  considered overlapping. Default is 15.

- genome:

  The genome to use, from the BSgenome package (e.g.,
  BSgenome.Mmusculus.UCSC.mm10).

- min.score:

  A numeric value specifying the minimum score for a PWM match to be
  considered significant. Default is 0.9.

- filter:

  A non-negative integer specifying the threshold of matches present for
  a TE-chimeric transcripts. This allows for filtering out transcripts
  with a low number of significant motif matches, focusing the analysis
  on more likely candidates. Default is 0.
