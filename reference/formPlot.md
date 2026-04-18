# Plot the formation of a specific transcript

Plot the formation of a specific transcript

## Usage

``` r
formPlot(GTF, txid = NULL, rank = 1, geneName = "Gene", TEname = "MT2B2")
```

## Arguments

- GTF:

  A Granges object, extracted from the output of
  stringtieCombine[`stringtieCombine`](https://yeehanxiao.github.io/Teddy/reference/stringtieCombine.md),
  filtered to rows with the type "exon". Use the `loadData` function to
  view the specific format.

- txid:

  The transcript ID of the TE-chimeric transcript to be plotted.

- rank:

  Integer. The index of the exon in which the TE element is located in
  the chimeric event. This information can be obtained from the
  `annotation` or the output of
  [`stringtieCombine`](https://yeehanxiao.github.io/Teddy/reference/stringtieCombine.md).

- geneName:

  The gene name to which the transcript belongs.

- TEname:

  The TE element involved in the fusion.
