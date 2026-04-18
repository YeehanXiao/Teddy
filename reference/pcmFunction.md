# Convert Probability Matrix to Position Weight Matrix

Convert Probability Matrix to Position Weight Matrix

## Usage

``` r
pcmFunction(pcm)
```

## Arguments

- pcm:

  A numeric matrix with 4 rows representing nucleotide probabilities for
  'A', 'C', 'G', and 'T'. Each column represents a position in the
  motif. The values should be proportions and each column should sum up
  to 1.

## Value

A PWM object that can be used in further analysis with motif matching
functions.
