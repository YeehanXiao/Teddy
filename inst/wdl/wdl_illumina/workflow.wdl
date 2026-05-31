version 1.0
import "map.wdl" as map
workflow runSTAR{
  input {
    Array[String] samples
    String inputs
    String index
    Int threads
    String GTF
    String bamdir
    String countpath
    String bwpath
    String rmpath
    String sortpath
    String picardpath
  }
    scatter (sample in samples) { 
      call map.STARbam {input: sampleName = sample, outputs = bamdir, index = index,inputs = inputs, threads = threads}
      call map.rmdup {input: sampleName = sample, bampath = STARbam.bampath, rmpath = rmpath, picardpath = picardpath}
      call map.bamsort {input: sampleName = sample, bampath = rmdup.newbampath, sortpath = sortpath, picardpath = picardpath}
      call map.buildindex {input: sampleName = sample, bampath = bamsort.outpath, picardpath = picardpath} 
      call map.featureCount {input: GTF = GTF, bampath =  bamsort.outpath, countpath = countpath, countName = sample, threads = threads}
      call map.bam2bw {input: sampleName = sample, bampath = buildindex.bamindexpath, bwpath = bwpath}
     }
    
  }
  