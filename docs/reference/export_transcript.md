# Export transcript sequences

Export spliced transcript sequences by retrieving exon sequences from an
exon-level GTF \`GRanges\` object. Exons can either be concatenated into
full transcript sequences or exported separately.

## Usage

``` r
export_transcript(
  txid,
  GTF,
  genome,
  outdir = "transcript_fasta",
  write = TRUE,
  collapse = TRUE
)
```

## Arguments

- txid:

  Character vector. Transcript IDs to export.

- GTF:

  A \`GRanges\` object containing exon-level transcript annotation, such
  as the output of \`extractGTF(type = "exon")\`.

- genome:

  Genome sequence object accepted by \`Biostrings::getSeq()\`, such as a
  \`BSgenome\` object or \`FaFile\`.

- outdir:

  Character. Output directory for FASTA files. Default is
  \`"transcript_fasta"\`.

- write:

  Logical. If \`TRUE\`, write sequences to FASTA files and invisibly
  return output file paths. If \`FALSE\`, return a \`DNAStringSet\`.

- collapse:

  Logical. If \`TRUE\`, concatenate exon sequences into one spliced
  transcript sequence. If \`FALSE\`, return or write exon sequences
  separately.

## Value

If \`write = TRUE\`, invisibly returns output FASTA file paths. If
\`write = FALSE\`, returns a \`DNAStringSet\`.
