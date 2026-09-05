# Annotate significant chimeric transcript test results

Extract significant chimeric transcripts, compute fold change, match
gene symbol and transcript annotation including exon rank.

## Usage

``` r
annotateChiTestResult(chi_test, gtf, se, cutoff = 0.2)
```

## Arguments

- chi_test:

  A DESeqDataSet object returned from ChimericDrivenTest().

- gtf:

  A GRanges object returned from extractGTF().

- se:

  A SummarizedExperiment object.

- cutoff:

  Adjusted p-value cutoff to extract the significant results.

## Value

A data.frame with extracted test results, fold change, gene symbol, and
exon annotation.
