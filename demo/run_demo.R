suppressPackageStartupMessages({
  library(Teddy)
  library(SummarizedExperiment)
  library(parallel)
})

demo_dir <- normalizePath("demo/data")
work_dir <- normalizePath("demo/output", mustWork = FALSE)

gtf_dir  <- file.path(work_dir, "gtf")
GTF_dir  <- file.path(work_dir, "GTF")
meta_dir <- file.path(work_dir, "meta")

dir.create(gtf_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(GTF_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

reference <- file.path(demo_dir, "chr8_demo.gtf")
te <- readRDS(file.path(demo_dir, "chr8_demo_TE.rds"))

bamfiles <- c(
  biopos = file.path(demo_dir, "chr8_demo_biopos.bam"),
  bioneg = file.path(demo_dir, "chr8_demo_bioneg.bam")
)

stopifnot(
  file.exists(reference),
  file.exists(bamfiles),
  file.exists(paste0(bamfiles, ".bai"))
)

## 1. Per-sample transcript reconstruction
gtffiles <- mapply(
  FUN = function(bam, sample) {
    outfile <- file.path(gtf_dir, paste0(sample, ".gtf"))

    Teddy::stringtieAssembly(
      bam = bam,
      reference = reference,
      outfile = outfile,
      params = "-p 2"
    )

    outfile
  },
  bam = unname(bamfiles),
  sample = names(bamfiles),
  USE.NAMES = FALSE
)

## 2. Merge reconstructed transcriptomes
merged_gtf <- file.path(GTF_dir, "merged.gtf")

Teddy::stringtieMerge(
  reference = reference,
  gtfFiles = gtffiles,
  outfile = merged_gtf,
  params = "-p 2 -F 0 -T 0"
)

## 3. Compare and annotate against the reference
annotated_prefix <- file.path(GTF_dir, "annotated.gtf")

Teddy::gffcompareAnno(
  reference = reference,
  gtffile = merged_gtf,
  outfile = annotated_prefix,
  overwrite = TRUE
)

annotated_candidates <- c(
  annotated_prefix,
  file.path(GTF_dir, "Teddy.annotated.gtf")
)

annotated_gtf <- annotated_candidates[file.exists(annotated_candidates)][1]

if (is.na(annotated_gtf)) {
  stop("Annotated GTF was not generated.")
}

## 4. Generate TE-annotated exonic bins
anno <- Teddy::prepareAnno(
  gtffile = annotated_gtf,
  transposon = te,
  cores = 2
)

## 5. Transcript-level quantification across both samples
combineSE <- Teddy::stringtieCombine(
  reference = annotated_gtf,
  params = "-p 2",
  gtfFiles = gtffiles,
  bamFiles = unname(bamfiles),
  longRead = FALSE,
  cores = 2
)

## 6. Bin-level quantification
se <- Teddy::countAnno(
  annotation = anno,
  bamfile = unname(bamfiles),
  nthreads = 2
)

## 7. Extract TE-chimeric transcript features
chi_GTF <- Teddy::processGTF(
  te = te,
  combineSE = combineSE,
  minoverlap = 5,
  threads = 2
)

saveRDS(anno,      file.path(meta_dir, "anno.rds"))
saveRDS(combineSE, file.path(meta_dir, "combineSE.rds"))
saveRDS(se,        file.path(meta_dir, "se.rds"))
saveRDS(chi_GTF,   file.path(meta_dir, "chi_GTF.rds"))

summary_table <- data.frame(
  output = c(
    "TE-annotated exonic bins",
    "transcript-level features",
    "bin-level features",
    "TE-chimeric transcript features"
  ),
  n = c(
    length(anno),
    nrow(combineSE),
    nrow(se),
    length(chi_GTF)
  )
)

write.table(
  summary_table,
  file.path(work_dir, "demo_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(summary_table)
cat("\nTEDDY demo completed successfully.\n")
