version 1.0

task Minimap2Index {
  input {
    File ref_fa
    Int threads = 8
    String index_dir   
  }

  String mmi_name = basename(ref_fa, ".fa") + ".mmi"
  String mmi_path = index_dir + "/" + mmi_name

  command <<<
    set -euo pipefail
    mkdir -p "~{index_dir}"
    minimap2 -t ~{threads} -d "~{mmi_path}" "~{ref_fa}"
  >>>

  output {
    File mmi = mmi_path
  }

  parameter_meta {
    ref_fa: { description: "Reference FASTA" }
    index_dir: { description: "Directory to store minimap2 index" }
  }
}

task Minimap2Map {
  input {
    Boolean do_map = true
    File    mmi
    String  fastqdir        # 
    String  sampleName
    String  tail = ".fq.gz" # ".fastq.gz" / ".fq.gz"
    String  platform        # "ont_drna" | "ont_cdna" | "pacbio_hifi"
    String  outputs
    Int     threads = 16
  }

  String outdir = sub(outputs, "//$", "")
  String fastq  = fastqdir + "/" + sampleName + tail
  String outbam = outdir + "/" + sampleName + ".bam"

  command <<<
    set -euo pipefail
    if [[ "~{do_map}" == "true" ]]; then
      mkdir -p "~{outdir}"
      preset=""
      case "~{platform}" in
        ont_drna)    preset="-ax splice -uf -k14" ;;
        ont_cdna)    preset="-ax splice -k14"     ;;
        pacbio_hifi) preset="-ax splice:hq"       ;;
        *) echo "Unsupported platform: ~{platform}" >&2; exit 2 ;;
      esac
         minimap2 -t ~{threads} $preset "~{mmi}" "~{fastq}" \
         | samtools view -@ ~{threads} -bS - > "~{outbam}"
    else
      :
    fi
  >>>


  output {
    File? bamout = select_first([outbam])
  }


  parameter_meta {
    platform:  {description: "ont_drna | ont_cdna | pacbio_hifi"}
    fastqdir:  {description: "Directory containing single-end long-read FASTQ (.fq.gz)"}
    sampleName:{description: "Sample name, will auto append .fq.gz"}
  }
}
task UseExistingBam {
  input {
    Boolean do_use = false
    File    bam_file
    String  sample
  }
  command <<<
    set -euo pipefail
    if [[ "~{do_use}" == "true" ]]; then
      test -s "~{bam_file}"
    fi
  >>>
  output {
    File? bamout = bam_file
  }
}

task rmdup {
  input {
    String sampleName
    File   inbam                 # 
    String rmpath                # 
    Int    threads = 8
    Boolean is_paired = false    # ONT/PacBio false；Illumina PE true
  }

  String outdir = sub(rmpath, "//$", "")
  String tmpdir = "tmp_" + sampleName
  String outbam = outdir + "/" + sampleName + ".bam"

  command <<<
    set -euo pipefail
    mkdir -p "~{outdir}" "~{tmpdir}"

    if [[ "~{is_paired}" == "true" ]]; then
      # name sort → fixmate → coordinate sort → markdup -r
     samtools sort -@ ~{threads} -n -o "~{tmpdir}/~{sampleName}.qname.bam" "~{inbam}"
     samtools fixmate -m "~{tmpdir}/~{sampleName}.qname.bam" "~{tmpdir}/~{sampleName}.fixmate.bam"
     samtools sort -@ ~{threads} -o "~{tmpdir}/~{sampleName}.coord.bam" "~{tmpdir}/~{sampleName}.fixmate.bam"
     samtools markdup -@ ~{threads} -r "~{tmpdir}/~{sampleName}.coord.bam" "~{outbam}"
    else
      # → markdup -r
      samtools sort -@ ~{threads} -o "~{tmpdir}/~{sampleName}.coord.bam" "~{inbam}"
      samtools markdup -@ ~{threads} -r "~{tmpdir}/~{sampleName}.coord.bam" "~{outbam}"
    fi

    samtools index -@ ~{threads} "~{outbam}"
    rm -rf "~{tmpdir}"
  >>>

  output {
    File rmdupbam = outbam
    File bamindex = outbam + ".bai"
  }
}

