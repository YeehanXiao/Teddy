version 1.0
task STARbam {
   input {
      Int threads
      String index
      String outputs
      String inputs
      String sampleName
    }

   String outpath = sub(outputs, "//$", "") 
   File Read1 = "~{inputs}/~{sampleName}.R1.fastq.gz"
   File Read2 = "~{inputs}/~{sampleName}.R2.fastq.gz"
    
  command <<<
     mkdir -p ~{outpath}
     STAR --runThreadN ~{threads} --genomeDir ~{index} --readFilesCommand zcat --readFilesIn ~{Read1} ~{Read2} --outFileNamePrefix ~{outpath}/~{sampleName} --outSAMtype BAM SortedByCoordinate 
  >>>


  output {
      File bamout = "~{outpath}/~{sampleName}Aligned.sortedByCoord.out.bam"
      String bampath = outpath
    }

  parameter_meta {
        threads:{description: "Threads used", category: "comment"}
        index: {description: "Index of genome reference", category: "required"}
        inputs: {description: "input file path", category: "required"}
        outputs: {description: "output file path", category: "required"}
        sampleName: {description: "Name of the sample", category: "required"}
    }
}

task featureCount {
   input{
       Int threads 
       String GTF
       String bampath
       String countpath
       String countName
   }

   String bamFile = "~{bampath}/~{countName}.bam"

   command <<<
    mkdir -p ~{countpath}
    cd ~{bampath}
    featureCounts -T ~{threads} -a ~{GTF} -o ~{countpath}/~{countName}.txt -p -B -C -f -t gene -g gene_id  ~{bamFile}
    >>>

   output {
      File count = "~{countpath}/~{countName}.txt"
   }

   parameter_meta {
        threads:{description: "Threads used", category: "comment"}
        GTF: {description: "GTF file for annotation", category: "required"}
        countpath: {description: "output file path", category: "required"}
        bampath: {description:"Dir of bam files",category:"required"}
        countName: {description: "name of count", category: "required"}
   }
}


task rmdup {
  input {
    String sampleName
    String bampath
    String rmpath
    String picardpath
  }

  command <<<
    mkdir -p ~{rmpath}
    java -jar ~{picardpath} MarkDuplicates \
      --REMOVE_DUPLICATES TRUE \
      -I ~{bampath}/~{sampleName}Aligned.sortedByCoord.out.bam \
      -O ~{rmpath}/~{sampleName}.bam \
      -M ~{rmpath}/~{sampleName}.txt \
      --ASSUME_SORT_ORDER coordinate
  >>>

  output {
    File rmdupbam = "~{rmpath}/~{sampleName}.bam"
    String newbampath = rmpath
  }
}

task bamsort{
    input {
      String sampleName
      String picardpath
      String bampath
      String sortpath
  }

  command <<<
    mkdir -p ~{sortpath}
    java -jar  ~{picardpath} SortSam -I ~{bampath}/~{sampleName}.bam -O ~{sortpath}/~{sampleName}.bam --SORT_ORDER coordinate
   >>>

  output {
    String outpath = sortpath
    File sortbam = "~{sortpath}/~{sampleName}.bam"
 } 

}

task buildindex{
  input {
     String sampleName
     String picardpath
     String bampath
  }
  
  command <<<
      mkdir -p ~{bampath}
      java -jar  ~{picardpath} BuildBamIndex -I ~{bampath}/~{sampleName}.bam
   >>>
   output {
     String bamindexpath = bampath
   }
}


task bam2bw{
  input {
    String sampleName
    String bampath
    String bwpath
  }

  command <<<
      mkdir -p ~{bwpath}
      bamCoverage -b ~{bampath}/~{sampleName}.bam -p 5 -o ~{bwpath}/~{sampleName}_forward.bw --filterRNAstrand forward --normalizeUsing RPKM
      bamCoverage -b ~{bampath}/~{sampleName}.bam -p 5 -o ~{bwpath}/~{sampleName}_reverse.bw --filterRNAstrand reverse --normalizeUsing RPKM
  >>>
}
