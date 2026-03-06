# ==============================================================================
# Script: TEDDYrunSimulatedData_benchmark.R
# Purpose: Benchmark TEDDY on simulated datasets (5x–100x depths)
#
# Design:
#   - Full simulated transcriptome (100%) used as ground truth
#   - 90% randomly retained annotation used as incomplete reference
#   - Reconstruction performed under incomplete annotation scenario
#   - Accuracy evaluated against full ground truth
# ==============================================================================

library(parallel)

# ------------------------------------------------------------------------------
# 1. Environment & References
# ------------------------------------------------------------------------------

setwd("~/manuals/compare/simulate/sortbam")

bamfiles <- list.files(pattern="*.bam$", full.names=TRUE)

# ---- Complete ground truth (100% simulated transcripts)
truth_gtf <- "/mnt/datadisk/xiaoyihan/manuals/compare/simulate/1000_simulated_transcripts.gtf"

# ---- Incomplete reference annotation (90% retained; randomly subsampled)
reference_gtf <- "/mnt/datadisk/xiaoyihan/manuals/compare/simulate/narrowed_1000_simulated_transcripts_with_genes.gtf"

te <- readRDS("/mnt/datadisk/xiaoyihan/manuals/meta/te.rds")

dir.create("../gtf", showWarnings=FALSE)

# ------------------------------------------------------------------------------
# 2. Per-sample assembly (guided by incomplete reference)
# ------------------------------------------------------------------------------

mclapply(
  bamfiles,
  FUN = function(bam, reference_gtf) {
    
    outfile <- file.path("../gtf",
                         sub(".bam$", ".gtf", basename(bam)))
    
    Teddy::stringtieAssembly(
      bam       = bam,
      reference = reference_gtf,
      gtfFiles  = outfile,
      params    = "-p 5"
    )
    
  },
  reference_gtf = reference_gtf,
  mc.cores = 20
)

gtffiles <- list.files("../gtf", pattern="*.gtf$", full.names=TRUE)

# ------------------------------------------------------------------------------
# 3. Depth-wise merge, annotation, and quantification
# ------------------------------------------------------------------------------

depths <- c("5x","10x","25x","50x","100x")

results <- lapply(depths, function(d){
  
  message("Processing depth: ", d)
  
  depth_gtfs <- gtffiles[grep(paste0("_", d), gtffiles)]
  depth_bams <- bamfiles[grep(paste0("_", d), bamfiles)]
  
  merged_gtf <- file.path("../gtf", paste0(d, ".gtf"))
  
  # ---- Merge assemblies ----
  Teddy::stringtieMerge(
    reference = reference_gtf,
    gtfFiles  = depth_gtfs,
    outfile   = merged_gtf,
    params    = "-p 20"
  )
  
  # ---- Compare to full ground truth ----
  Teddy::gffcompareAnno(
    reference = truth_gtf,
    gtffile   = merged_gtf,
    outfile   = sub(".gtf$", ".annotated.gtf", merged_gtf)
  )
  
  annotated_gtf <- sub(".gtf$", ".annotated.gtf", merged_gtf)
  
  # ---- TE annotation ----
  anno <- Teddy::prepareAnno(
    gtffile    = annotated_gtf,
    transposon = te
  )
  
  # ---- Quantification ----
  se <- Teddy::countAnno(
    annotation = anno,
    bamfile    = depth_bams
  )
  
  combineSE <- Teddy::stringtieCombine(
    reference = annotated_gtf,
    bamfile   = depth_bams,
    params    = "-p 20",
    gtfFiles  = depth_gtfs
  )
  
  list(
    se        = se,
    combineSE = combineSE
  )
})

names(results) <- depths

# ------------------------------------------------------------------------------
# 4. Construct truth set (TE-containing transcripts from full simulated truth)
# ------------------------------------------------------------------------------

library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(tibble)

TF <- import(truth_gtf)

truth_exons <- TF[TF$type == "exon"]
ov_truth <- findOverlaps(truth_exons, te, ignore.strand = TRUE, minoverlap = 1)

truth_txids <- unique(truth_exons$transcript_id[queryHits(ov_truth)])
truth_gtf_te <- TF[TF$transcript_id %in% truth_txids]
txdb_truth <- makeTxDbFromGRanges(truth_gtf_te)

# ------------------------------------------------------------------------------
# 5. Helper functions used by evaluate_one_adaptive (keep original logic)
# ------------------------------------------------------------------------------

get_pred_gtf <- function(SE) {
  g <- S4Vectors::metadata(SE)$gtf
  if (is.null(g)) stop("metadata(SE)$gtf not found in combineSE object.")
  g
}

ensure_tx_rownames <- function(SE) {
  a <- assay(SE)
  if (!is.null(rownames(a)) && all(nzchar(rownames(a)))) return(SE)
  rr <- rowRanges(SE)
  if (!"transcript_id" %in% colnames(as.data.frame(rr))) {
    stop("rowRanges(SE) lacks transcript_id.")
  }
  rownames(a) <- as.character(rr$transcript_id)
  assay(SE) <- a
  SE
}

get_rowmean_cpm <- function(SE) {
  rowMeans(assay(SE), na.rm = TRUE)
}

txids_with_TE <- function(gtf, te) {
  ex <- gtf[gtf$type == "exon"]
  ov <- findOverlaps(ex, te, ignore.strand = TRUE, minoverlap = 1)
  unique(ex$transcript_id[queryHits(ov)])
}

te_overlap_bp_per_tx <- function(gtf, te) {
  ex <- gtf[gtf$type == "exon"]
  if (length(ex) == 0L) return(setNames(numeric(0), character(0)))
  ov <- findOverlaps(ex, te, ignore.strand = TRUE, minoverlap = 1)
  if (length(ov) == 0L) return(setNames(numeric(0), character(0)))
  
  p <- pintersect(ex[queryHits(ov)], te[subjectHits(ov)], ignore.strand = TRUE)
  bp <- tapply(width(p), ex$transcript_id[queryHits(ov)], sum)
  bp[is.na(bp)] <- 0
  bp
}

clean_gtf_for_txdb <- function(gtf) {
  gtf[gtf$type %in% c("transcript", "exon")]
}

junction_df <- function(txdb) {
  intr <- intronsByTranscript(txdb, use.names = TRUE)
  if (length(intr) == 0L) return(data.frame())
  as.data.frame(intr) |>
    dplyr::transmute(
      tx = group_name,
      chr = as.character(seqnames),
      donor = start - 1L,
      accept = end
    )
}

match_by_jaccard <- function(dfA, dfB, tol = 5L) {
  if (!nrow(dfA) || !nrow(dfB)) {
    return(data.frame(tx=character(), best_truth=character(), jaccard=numeric()))
  }
  key <- function(df) paste(df$chr, round(df$donor/tol), round(df$accept/tol))
  dfA$key <- key(dfA)
  dfB$key <- key(dfB)
  
  splitA <- split(dfA$key, dfA$tx)
  splitB <- split(dfB$key, dfB$tx)
  
  res <- lapply(names(splitA), function(txa){
    a <- unique(splitA[[txa]])
    scores <- vapply(splitB, function(b){
      b <- unique(b)
      inter <- length(intersect(a, b))
      uni <- length(union(a, b))
      if (uni == 0) 0 else inter / uni
    }, numeric(1))
    bt <- names(which.max(scores))
    data.frame(
      tx = txa,
      best_truth = if (length(bt)) bt else NA_character_,
      jaccard = if (length(scores)) max(scores) else 0
    )
  })
  do.call(rbind, res)
}

# ------------------------------------------------------------------------------
# 6. Keep evaluate_one_adaptive (original logic)
# ------------------------------------------------------------------------------

evaluate_one_adaptive <- function(SE, te, txdb_truth, truth_txids,
                                  depth_label=c("5x","10x","25x","50x","100x"),
                                  mode=c("isoform","gene")) {
  mode <- match.arg(mode)
  plan <- depth_strategy(depth_label)
  q_expr   <- plan$q_expr
  min_junc <- plan$min_junc
  j_cut    <- plan$j_cut
  te_min   <- plan$te_min_bp
  topK     <- plan$topK
  
  SE <- ensure_tx_rownames(SE)
  pred_gtf <- get_pred_gtf(SE)
  
  tx_te_all <- txids_with_TE(pred_gtf, te)
  if (length(tx_te_all)) {
    te_bp <- te_overlap_bp_per_tx(pred_gtf[pred_gtf$transcript_id %in% tx_te_all], te)
    tx_te <- names(te_bp)[te_bp >= te_min]
  } else tx_te <- character(0)
  
  rm_cpm <- get_rowmean_cpm(SE)
  thr    <- quantile(rm_cpm, q_expr, na.rm=TRUE)
  active <- names(rm_cpm)[rm_cpm >= thr]
  
  pred_te_gtf_all <- pred_gtf[pred_gtf$transcript_id %in% tx_te]
  pred_te_gtf_all <- clean_gtf_for_txdb(pred_te_gtf_all)
  if (length(pred_te_gtf_all)==0L)
    return(list(TP=0,FP=0,FN=length(truth_txids),Precision=0,Recall=0,F1=0))
  
  txdb_pred_all <- txdbmaker::makeTxDbFromGRanges(pred_te_gtf_all)
  df_pred_all   <- junction_df(txdb_pred_all)
  if (!nrow(df_pred_all))
    return(list(TP=0,FP=0,FN=length(truth_txids),Precision=0,Recall=0,F1=0))
  
  junc_n <- table(df_pred_all$tx)
  
  keep_tx <- intersect(active, names(junc_n)[junc_n >= min_junc])
  keep_tx <- intersect(keep_tx, unique(df_pred_all$tx))
  if (!length(keep_tx))
    return(list(TP=0,FP=0,FN=length(truth_txids),Precision=0,Recall=0,F1=0))
  
  if (topK > 0) {
    rr <- pred_te_gtf_all[mcols(pred_te_gtf_all)$type=="transcript"]
    tx2gene <- as.data.frame(mcols(rr)[,c("transcript_id","gene_id")])
    tx2gene <- tx2gene[!duplicated(tx2gene$transcript_id),]
    rm_sub  <- rm_cpm[keep_tx]
    gmap <- tx2gene$gene_id[match(names(rm_sub), tx2gene$transcript_id)]
    ord <- order(gmap, -rm_sub, na.last=NA)
    df  <- data.frame(tx=names(rm_sub)[ord], gene=gmap[ord], expr=rm_sub[ord])
    df  <- df[!is.na(df$gene),]
    df$rank <- ave(df$expr, df$gene, FUN=function(x) rank(-x, ties.method="first"))
    keep_tx <- df$tx[df$rank <= topK]
  }
  
  df_pred <- df_pred_all[df_pred_all$tx %in% keep_tx, , drop=FALSE]
  df_truth <- junction_df(txdb_truth)
  
  mj <- match_by_jaccard(df_pred, df_truth, tol = 5L)
  TP_tx <- mj$tx[which(mj$jaccard >= j_cut)]
  
  if (mode=="isoform") {
    TP <- length(TP_tx)
    FP <- length(setdiff(unique(df_pred$tx), TP_tx))
    matched_truth <- mj$best_truth[which(mj$jaccard >= j_cut)]
    FN <- length(setdiff(truth_txids, matched_truth))
  } else {
    tx2gene <- as.data.frame(mcols(pred_te_gtf_all[mcols(pred_te_gtf_all)$type=="transcript"])[,c("transcript_id","gene_id")])
    tx2gene <- tx2gene[!duplicated(tx2gene$transcript_id),]
    hit_gene_pred <- unique(tx2gene$gene_id[match(unique(df_pred$tx), tx2gene$transcript_id)])
    hit_gene_TP   <- unique(tx2gene$gene_id[match(TP_tx, tx2gene$transcript_id)])
    
    truth_tx2gene <- as.data.frame(transcripts(txdb_truth, columns="gene_id"))[,c("tx_name","gene_id")]
    truth_gene_set <- unique(truth_tx2gene$gene_id[match(truth_txids, truth_tx2gene$tx_name)])
    
    TP <- length(intersect(hit_gene_TP, truth_gene_set))
    FP <- length(setdiff(hit_gene_pred, truth_gene_set))
    FN <- length(setdiff(truth_gene_set, hit_gene_TP))
  }
  
  P <- ifelse(TP+FP==0,0,TP/(TP+FP))
  R <- ifelse(TP+FN==0,0,TP/(TP+FN))
  F1 <- ifelse(P+R==0,0,2*P*R/(P+R))
  list(TP=TP,FP=FP,FN=FN,Precision=P,Recall=R,F1=F1)
}

# ------------------------------------------------------------------------------
# 7. Fixed default parameters (keep function name, remove adaptive tuning)
# ------------------------------------------------------------------------------

depth_strategy <- function(depth_label){
  list(
    q_expr    = 0.35,
    min_junc  = 1,
    j_cut     = 0.82,
    te_min_bp = 10,
    topK      = 1
  )
}

# ------------------------------------------------------------------------------
# 8. Run TEDDY evaluation on results[[depth]]$combineSE
# ------------------------------------------------------------------------------

teddy_metrics <- lapply(names(results), function(d) {
  evaluate_one_adaptive(
    SE          = results[[d]]$combineSE,
    te          = te,
    txdb_truth  = txdb_truth,
    truth_txids = truth_txids,
    depth_label = d,
    mode        = "isoform"
  )
})
names(teddy_metrics) <- names(results)

teddy_benchmark <- do.call(rbind, lapply(names(teddy_metrics), function(d) {
  m <- teddy_metrics[[d]]
  data.frame(
    Depth     = d,
    Tool      = "TEDDY",
    TP        = m$TP,
    FP        = m$FP,
    FN        = m$FN,
    Precision = m$Precision,
    Recall    = m$Recall,
    F1        = m$F1,
    row.names = NULL
  )
}))


teddy_benchmark$Depth <- factor(teddy_benchmark$Depth, levels = depths)
teddy_benchmark <- teddy_benchmark[order(teddy_benchmark$Depth), ]

