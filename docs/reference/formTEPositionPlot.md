# Plot TE position within a transcript structure

Plot TE position within a transcript structure

## Usage

``` r
formTEPositionPlot(
  GTF,
  TEposition,
  txid,
  geneName = NULL,
  showTElabel = TRUE,
  clusterGap = 0.015,
  exonFill = "#ECD1D6",
  TEfill = "#3A3B4F",
  borderColor = "black",
  TEborderColor = "white"
)
```

## Arguments

- GTF:

  A GRanges object containing the complete transcript annotation.

- TEposition:

  A GRanges object returned by
  [`annotateTEPosition`](https://yeehanxiao.github.io/Teddy/reference/annotateTEPosition.md).

- txid:

  Character. Transcript ID to plot.

- geneName:

  Character. Gene name. If NULL, inferred from GTF.

- showTElabel:

  Logical. Whether to show TE labels.

- clusterGap:

  Numeric. Maximum plotting distance between adjacent TE segments to
  group them into the same visual cluster.

- exonFill:

  Character. Fill color for host-derived exon regions.

- TEfill:

  Character. Fill color for TE-derived regions.

- borderColor:

  Character. Color for exon borders, introns, and label lines.

- TEborderColor:

  Character. Border color separating adjacent TE segments.

## Value

A ggplot object.
