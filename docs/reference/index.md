# Package index

## Transcriptome reconstruction and quantification

- [`stringtieAssembly()`](https://yeehanxiao.github.io/Teddy/reference/stringtieAssembly.md)
  : R wrapper to Run stringtie
- [`stringtieMerge()`](https://yeehanxiao.github.io/Teddy/reference/stringtieMerge.md)
  : R wrapper to Run stringtie tool
- [`gffcompareAnno()`](https://yeehanxiao.github.io/Teddy/reference/gffcompareAnno.md)
  : R wrapper to run gffcompare
- [`stringtieCombine()`](https://yeehanxiao.github.io/Teddy/reference/stringtieCombine.md)
  : Consolidation of information on transcript assembly of multiple
  samples
- [`buildIsoformSupport()`](https://yeehanxiao.github.io/Teddy/reference/buildIsoformSupport.md)
  : Build an isoform support matrix across multiple GTF files

## TE-chimeric feature annotation

- [`prepareAnno()`](https://yeehanxiao.github.io/Teddy/reference/prepareAnno.md)
  : Preparing the genome annotation object
- [`countAnno()`](https://yeehanxiao.github.io/Teddy/reference/countAnno.md)
  : Counting reads on exon bins
- [`processGTF()`](https://yeehanxiao.github.io/Teddy/reference/processGTF.md)
  : Process GTF and SummarizedExperiment for TE overlap analysis

## Differential TE-chimeric exon usage

- [`ChimericDrivenTest()`](https://yeehanxiao.github.io/Teddy/reference/ChimericDrivenTest.md)
  : Testing for differential TE-chimeric exon usage
- [`extractTest()`](https://yeehanxiao.github.io/Teddy/reference/extractTest.md)
  : Extract the result from the differentially expressed TE-chimeric
  exon test
- [`calculateFoldchange()`](https://yeehanxiao.github.io/Teddy/reference/calculateFoldchange.md)
  : Calculating the fold change
- [`pick_canonical_tx()`](https://yeehanxiao.github.io/Teddy/reference/pick_canonical_tx.md)
  : Pick the most expressed reference transcript for a gene
- [`pick_chimeric_tx()`](https://yeehanxiao.github.io/Teddy/reference/pick_chimeric_tx.md)
  : Pick TE-chimeric transcript(s) for a gene

## Visualization

- [`formPlot()`](https://yeehanxiao.github.io/Teddy/reference/formPlot.md)
  : Plot the formation of a specific transcript
- [`diffBinPlot()`](https://yeehanxiao.github.io/Teddy/reference/diffBinPlot.md)
  : Plot the gene model in the control and reference groups
- [`plotIsoformDynamics()`](https://yeehanxiao.github.io/Teddy/reference/plotIsoformDynamics.md)
  : Plot Isoform Expression Dynamics

## Sequence and motif analysis

- [`export_transcript()`](https://yeehanxiao.github.io/Teddy/reference/export_transcript.md)
  : Export transcript sequences
- [`MotifSearch()`](https://yeehanxiao.github.io/Teddy/reference/MotifSearch.md)
  : Search the motif binding sites through TE-chimeric transcripts
- [`pcmFunction()`](https://yeehanxiao.github.io/Teddy/reference/pcmFunction.md)
  : Convert Probability Matrix to Position Weight Matrix

## Utilities

- [`extractGTF()`](https://yeehanxiao.github.io/Teddy/reference/extractGTF.md)
  : Extract GTF Based on Type and FPKM Filter
- [`clean_sampleNames()`](https://yeehanxiao.github.io/Teddy/reference/clean_sampleNames.md)
  : Clean sample names Remove replicate suffixes from sample names.
- [`NCBI_check()`](https://yeehanxiao.github.io/Teddy/reference/NCBI_check.md)
  : Standardize transposon annotations for TEDDY
