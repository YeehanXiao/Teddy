Short-read RNA-seq preprocessing WDL.

The workflow supports STAR-based alignment for Illumina paired-end RNA-seq data.

The WDL command blocks define the STAR alignment parameters and downstream samtools processing steps, including BAM sorting and indexing.
The template JSON uses placeholder paths, sample names, STAR genome index paths, and paired-end FASTQ inputs; users should replace them with local inputs.

The output sorted BAM files can be used directly as input for downstream TEDDY analysis.