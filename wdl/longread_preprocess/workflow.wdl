version 1.0
import "map.wdl" as map

workflow runMinimapLR_full {
  input {
    Array[String] samples
    String        platform

    String fastqdir
    String readSuffix = ".fastq.gz"
    String bamSuffix  = ".bam"

    File   ref_fa
    File?  index

    Boolean do_map = true
    Boolean do_use = false
    Boolean is_paired = false

    Int    threads = 16
    String outputs
    String out_bamdir
    String rmpath
    String index_dir
  }

  # Build mmi index if not provided
  if (!defined(index)) {
    call map.Minimap2Index as build_index {
      input:
        ref_fa    = ref_fa,
        threads   = threads,
        index_dir = index_dir
    }
  }
  File mmi = select_first([index, build_index.mmi])

  scatter (s in samples) {
    if (do_map) {
      call map.Minimap2Map as minimap2_map {
        input:
          mmi        = mmi,
          fastqdir   = fastqdir,
          sampleName = s,
          tail       = readSuffix,
          platform   = platform,
          outputs    = out_bamdir,
          threads    = threads
      }
    }

    if (do_use) {
      call map.UseExistingBam as use_bam {
        input:
          bam_file = sub(out_bamdir, "/$", "") + "/" + s + bamSuffix,
          sample   = s
      }
    }

    File in_bam = select_first([minimap2_map.bamout, use_bam.bamout])

    call map.rmdup {
      input:
        sampleName = s,
        inbam      = in_bam,
        rmpath     = rmpath,
        threads    = threads,
        is_paired  = is_paired
    }
  }

  output {
    Array[File] rmdup_bams = rmdup.rmdupbam
    Array[File] rmdup_bais = rmdup.bamindex
    File used_mmi          = mmi
  }
}