# =========================================================================
# Script: benchmark_ChimeraTEdata.R
# Purpose: Benchmark TEDDY against ChimeraTE using their official example data.
# Author: Yihan Xiao (Teddy Developer)
#
# Note: 
# The benchmark dataset used in this script is the official example data 
# from the ChimeraTE repository. You can download it directly from:
# https://github.com/GrzeB/ChimeraTE.git

# Please replace the absolute paths (e.g., "/mnt/datadisk/xiaoyihan/...") 
# with your local directory paths where the ChimeraTE example data is stored.
# =========================================================================

start_time <- Sys.time()

# -------------------------------------------------------------------------
# Part 1: Environment Setup & Data Loading
# ---------------------------------------------------------------------
library(Rstringtie)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
bamfiles <-list.files("/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE/example_data/mode1/", pattern = "*.bam$", full.names = TRUE)

# -------------------------------------------------------------------------
# Part 2: Assembly with Loose Parameters (Fair Comparison)
# -------------------------------------------------------------------------

params <- "-p 8 -f 0.02 -c 0.5 -j 1 -g 20"

mclapply(bamfiles, FUN = function(x, reference, params) {
  bam <- basename(x)
  outfile <- file.path(basename(x), gsub("\\.bam$", ".loose.gtf", bam))
  Teddy::stringtieAssembly(bam = x,
                           reference = reference,
                           outfile = outfile,
                           params = params)
}, reference = "/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE/example_data/mode1/test.gtf", params = params, mc.cores = 8)

gtffiles <- list.files("/mnt/datadisk/xiaoyihan/manuals/compare/", pattern = "loose",  full.names = TRUE)
stringtieMerge(reference = "/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE/example_data/mode1/test.gtf",
               gtfFiles = gtffiles,
               outfile = "../gtf/fly.gtf",
               params = "-p 8"
)

# -------------------------------------------------------------------------
# Part 3: TE Processing and Annotation using TEDDY Modules
# -------------------------------------------------------------------------
te <- read.delim("/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE/example_data/mode1/dmel_TEs_sample.gtf", header = FALSE, sep = "\t")

te_G <- GRanges(seqnames = te[, 1],
                IRanges(start = te[, 4],
                        end = te[, 5]),
                strand = te[, 7],
                names = te[, 9],class = te[,9])
processGTF(te = te_G, combineSE = combineSE,minoverlap = 0)

gffcompareAnno(reference="/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE/example_data/mode1/test.gtf", gtffile="../gtf/fly.gtf", outfile="../gtf/fly_loose.gtf")
gffcompareAnno(reference="/mnt/datadisk/xiaoyihan/manuals/compare/ChimeraTE/example_data/mode1/test.gtf", gtffile=gtffiles[2], outfile="../gtf/fly_loose_CR2.gtf")

anno <- prepareAnno(gtffile = "../gtf/fi_fly.annotated.gtf", transposon = te_G)


se <- countAnno(annotation = anno, bamfile = bamfiles)

combineSE <- stringtieCombine(reference ="../gtf/fi_fly.annotated.gtf", bamFiles = bamfiles ,
                              params = "-p 70", 
                              gtfFiles = gtffiles)
anno[anno$transposon!="none"]



# -------------------------------------------------------------------------
# Part 4: Load ChimeraTE Baseline and Evaluate
# -------------------------------------------------------------------------
# =========================================================================
# Note: 
# Load ChimeraTE baseline results for Venn diagram comparison.
# These expected gene lists were derived directly from the official 
# ChimeraTE example data README to ensure a fair and strict benchmark.
# =========================================================================
library(GenomicRanges)
expected_data <- readRDS("/Users/xiaoyihan/Core/MyTools/Develop_Final/Teddy/benchmark_reproducibility/ChiemraTEData_benchmark/intermediate/ChimeraTE_results.rds")
chimera_genes <- expected_data$chimera_genes
init_gr <- expected_data$init_gr
exon_gr <- expected_data$exon_gr
term_gr <- expected_data$term_gr

chimera_genes <- unique(c(init_gr$gene$gene_id,
                          exon_gr$gene$gene_id,
                          term_gr$gene$gene_id))

CHI_GTF <- processGTF(te = te_G, combineSE = combineSE,minoverlap = 0)

exon_count_df <- as.data.frame(mcols(GTF))%>%
  dplyr::group_by(transcript_id, gene_id, gene_name) %>%
  dplyr::summarise(
    exon_count = length(unique(exon_number)),
    max_exon_number = max(exon_number, na.rm = TRUE),
    .groups = "drop"
  )
mcols(CHI_GTF)$structure <- {tx_rank <- as.integer(as.character(mcols(CHI_GTF)$tx_exon_rank)); ec <- setNames(as.integer(exon_count_df$exon_count), exon_count_df$transcript_id); mx <- setNames(as.integer(exon_count_df$max_exon_number), exon_count_df$transcript_id); ecv <- ec[as.character(mcols(CHI_GTF)$transcript_id)]; ecv[is.na(ecv)] <- mx[as.character(mcols(CHI_GTF)$transcript_id)][is.na(ecv)]; ifelse(tx_rank==1, "TE-initiated", ifelse(!is.na(tx_rank) & !is.na(ecv) & tx_rank==ecv, "TE-terminated", "TE-exonized")) }


setdiff(unique(CHI_GTF$gene_name),unique(c(init_gr$gene$gene_id,exon_gr$gene$gene_id,term_gr$gene$gene_id)))
#[1] "FBgn0262731"
setdiff(unique(c(init_gr$gene$gene_id,exon_gr$gene$gene_id,term_gr$gene$gene_id)),CHI_GTF$gene_name)
#[1] "FBgn0031188" 

end_time <- Sys.time()
elapsed_time <- end_time - start_time
CHI_GTF <- processGTF(te = te_G, combineSE = combineSE,minoverlap = 0)
print(paste("Total script execution time: ", elapsed_time, " seconds"))
# =========================================================================
# Performance Note: 
# In our local testing environment, the core analytical pipeline 
# (including StringTie assembly, GTF merging, and TEDDY's TE-classification) 
# completed in approximately 5.35 seconds for this dataset, demonstrating 
# high computational efficiency.
# =========================================================================
# -------------------------------------------------------------------------
# Part 5: Visualization (Gviz Tracks)
# -------------------------------------------------------------------------
#---- Target Gene 1: FBgn0031188 ----
target_gene <- "FBgn0031188"
pad <- 5000


gene_exons <- GTF[mcols(GTF)$gene_name == target_gene]
gene_span  <- range(gene_exons)
plot_from  <- start(gene_span) - pad
plot_to    <- end(gene_span) + pad
plot_chr   <- as.character(seqnames(gene_span))

which_gr <- GRanges(plot_chr, IRanges(plot_from, plot_to))
bw_file <- "Benchmark/fly/CR1_reverse.bw"
which_gr <- GRanges(plot_chr, IRanges(plot_from, plot_to))

bw_gr <- import.bw(bw_file, which = which_gr)
bw_gr_pos <- import.bw("Benchmark/fly/CR1_forward.bw", which = which_gr)

bw_gr_orig <- import.bw("Benchmark/fly/CR1_Aligned.sortedByCoord.out.bw", which = which_gr)

data_dt <- DataTrack(range = bw_gr, genome = "", chromosome = plot_chr,
                     name = "CR_reverse_RNAseq",  type = "polygon", window = 38,
                     aggregation = "mean",
                     col = "black",fill= "black")
displayPars(data_dt) <- list(
  col = "black",        
  fill = FALSE,      
  lwd = 0.8,            
  alpha = 1,            
  showAxis = TRUE
)
data_dt_pos <- DataTrack(range = bw_gr_pos, genome = "", chromosome = plot_chr,
                         name = "CR1_coverage_PosStrand",  type = "polygon", window = 50,
                         aggregation = "mean",
                         col = "black")
data_dt_orig <- DataTrack(range = bw_gr_orig, genome = "", chromosome = plot_chr,
                          name = "CR1_coverage_PosStrand",  type = "polygon", window = 50,
                          aggregation = "mean",
                          col = "black")

Teddy_tx <- AnnotationTrack(
  GTF[GTF$transcript_id %in% c("MSTRG.82.1","FBtr0335047")], 
  name = "Tx",
  transcriptAnnotation = "",
  col = NA, fill = "#918579",shape = "box"
) 
Gene_track <- GeneRegionTrack( GTF[GTF$transcript_id %in% c("MSTRG.82.1")], transcript = c("MSTRG.82.1"),
                               fill = "darkblue", col = NA, col.line = "darkblue")
Gene_track_2 <- GeneRegionTrack( GTF[GTF$transcript_id %in% c("FBtr0335047")], transcript = c("FBtr0335047"),
                                 fill = "darkblue", col = NA, col.line = "darkblue")

S2_highlight <- HighlightTrack(trackList = list(data_dt),
                               range = GTF[GTF$transcript_id %in% c("MSTRG.82.1")][4,],
                               col = NA)

genomeAxis <- GenomeAxisTrack(name = "axis", col = "black")

annotation_te <- te_G[grep(pattern = "S2", te_G$names)]


te_annotationTrack <- AnnotationTrack(annotation_te,
                                      feature = annotation_te$names,
                                      col = NA, fill = "darkblue",
                                      shape = "arrow",
                                      showFeatureId = FALSE,name = "TE" )
displayPars(te_annotationTrack) <- list(
  col = "black",
  fontcolor.feature = "darkblue",
  cex.feature = 0.85,
  just.group = "left",
  lwd = 0.35,                   
  arrowHeadWidth = 1,         
  arrowHeadLength = 0.3,col = NA,fill = "darkblue"
)
plot_chr <- as.character(plot_chr)

bw_gr <- bw_gr[seqnames(bw_gr) == plot_chr]


plotTracks(list(genomeAxis,te_annotationTrack,Gene_track,Gene_track_2,data_dt), 
           from = plot_from, to = plot_to, 
           chromosome =plot_chr, scale = 0.2,
           geneSymbols = TRUE,
           transcriptAnnotation = "transcript")

#---- Target Gene 2: FBgn0262731 ----
target_gene <- "FBgn0262731"
pad <- 5000

gene_exons <- GTF[mcols(GTF)$gene_name == target_gene]
gene_span  <- range(gene_exons)
plot_from  <- start(gene_span) - pad
plot_to    <- end(gene_span) + pad
plot_chr   <- as.character(seqnames(gene_span))

which_gr <- GRanges(plot_chr, IRanges(plot_from, plot_to))
bw_gr <- import.bw(bw_file, which = which_gr)
data_dt <- DataTrack(range = bw_gr, genome = NA, chromosome = plot_chr,
                     name = "CR1_coverage", type = "histogram", window = "auto")


gene_exons <- GTF[mcols(GTF)$gene_name == target_gene]
gene_span  <- range(gene_exons)
plot_from  <- start(gene_span) - pad
plot_to    <- end(gene_span) + pad
plot_chr   <- as.character(seqnames(gene_span))

which_gr <- GRanges(plot_chr, IRanges(plot_from, plot_to))

bw_file <- "Benchmark/fly/CR2_reverse.bw"
which_gr <- GRanges(plot_chr, IRanges(plot_from, plot_to))

bw_gr <- import.bw(bw_file, which = which_gr)
bw_gr_pos <- import.bw("Benchmark/fly/CR2_forward.bw", which = which_gr)


data_dt <- DataTrack(range = bw_gr, genome = "", chromosome = plot_chr,
                     name = "CR2_coverage_NegStrand",  type = "polygon", window = 50,
                     aggregation = "mean",
                     col = "black",fill= "black")

data_dt_pos <- DataTrack(range = bw_gr_pos, genome = "", chromosome = plot_chr,
                         name = "CR1_coverage_PosStrand",  type = "polygon", window = 50,
                         aggregation = "mean",
                         col = "black")


Teddy_tx <- AnnotationTrack(
  GTF[GTF$transcript_id %in% c("FBtr0089254","FBtr0089253","MSTRG.12.3")], 
  name = "Tx",
  transcriptAnnotation = "",
  col = NA, fill = "#918579",shape = "box"
) 
Gene_track <- GeneRegionTrack( GTF[GTF$transcript_id %in% c("FBtr0089254")], transcript = c("FBtr0089254"),
                               fill = "darkblue", col = NA, col.line = "darkblue")

FB_highlight <- HighlightTrack(trackList = list(data_dt,data_dt_pos,data_dt_orig),
                               range = te_G[grep(pattern = "FB", te_G$names)][2],
                               col = NA)

genomeAxis <- GenomeAxisTrack(name = "axis", col = "black")

annotation_te <- te_G[grep(pattern = "FB", te_G$names)]
annotation_te <- annotation_te[2]

te_annotationTrack <- AnnotationTrack(annotation_te,
                                      feature = annotation_te$names,
                                      col = NA, fill = "darkblue",
                                      shape = "arrow",
                                      showFeatureId = TRUE,name = "FB" )





plotTracks(list(genomeAxis,te_annotationTrack,Teddy_tx,Gene_track,FB_highlight), 
           from = plot_from, to = plot_to, 
           chromosome =plot_chr, scale = 0.2,
           geneSymbols = TRUE,
           transcriptAnnotation = "transcript")
# -------------------------------------------------------------------------
# Part 6: Venn Diagram Output & Execution Time
# -------------------------------------------------------------------------


Teddy_genes <- unique(CHI_GTF$gene_name)


library(eulerr)
sets <- list(
  TEDDY = unique(Teddy_genes),
  ChimeraTE = unique(chimera_genes)
)
fit <- euler(sets)
plot(fit,
     fills = list(fill = c( "#EBC9C7","#918579"), alpha = 0.55),
     edges = TRUE,
     quantities = list(fontsize = 12, type = "counts"), 
     labels = NULL,
     main = "Benchmark on ChimeraTE fly sample ")

