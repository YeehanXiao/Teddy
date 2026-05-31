Long-read RNA-seq preprocessing WDL.

The workflow supports minimap2-based alignment for:
- ont_drna
- ont_cdna
- pacbio_hifi

The WDL command blocks define the minimap2 presets and downstream samtools processing steps.
The template JSON uses placeholder paths and sample names; users should replace them with local inputs.
