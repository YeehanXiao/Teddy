# ==============================================================================
# Script: ArribaRunSimulatedData_eventlevel.R
# Purpose: Build arriba_all + event-level benchmark for Arriba on simulated data
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

# ------------------------------------------------------------------------------
# 1) Paths (REPLACE if needed)
# ------------------------------------------------------------------------------

merge_dir <- "/mnt/datadisk/xiaoyihan/manuals/compare/simulate/merge"

# Full simulated truth (100%)
truth_gtf_path <- "/mnt/datadisk/xiaoyihan/manuals/compare/simulate/1000_simulated_transcripts.gtf"

# Incomplete reference (90% retained)
reference_gtf_path <- "/mnt/datadisk/xiaoyihan/manuals/compare/simulate/narrowed_1000_simulated_transcripts_with_genes.gtf"

# TE annotation used to build STAR/Arriba reference
te_anno_gtf_path <- "/mnt/datadisk/xiaoyihan/index/mm10_TE_annotations.gtf"

setwd(merge_dir)

# ------------------------------------------------------------------------------
# 2) Load Arriba TSVs -> arriba_all
# ------------------------------------------------------------------------------

files <- list.files(
  pattern = "^merge_\\d+x(_rep\\d+)?_arriba(_relaxed)?_fusions(\\.discarded)?\\.tsv$",
  full.names = TRUE
)
stopifnot(length(files) > 0)

read_arriba <- function(f) {
  dt <- suppressWarnings(fread(f, sep = "\t", header = TRUE, fill = TRUE, quote = ""))
  if (!nrow(dt)) return(NULL)
  
  setnames(dt, sub("^#", "", names(dt)))
  setnames(dt, sub("\\.+$", "", names(dt)))
  
  bn <- basename(f)
  sample <- sub("_arriba.*$", "", bn)
  depth  <- sub("^merge_", "", sub("_rep\\d+$", "", sample))
  rep    <- if (grepl("_rep\\d+$", sample)) sub("^.*_(rep\\d+)$", "\\1", sample) else "rep1"
  mode   <- if (grepl("_relaxed_", bn)) "relaxed" else "strict"
  status <- if (grepl("\\.discarded\\.tsv$", bn)) "discarded" else "kept"
  
  get_char <- function(nm) if (nm %in% names(dt)) as.character(dt[[nm]]) else rep(NA_character_, nrow(dt))
  get_int  <- function(nm) if (nm %in% names(dt)) suppressWarnings(as.integer(dt[[nm]])) else rep(NA_integer_, nrow(dt))
  
  geneA <- if ("gene_id1" %in% names(dt) && any(dt$gene_id1 != "." & !is.na(dt$gene_id1))) get_char("gene_id1") else get_char("gene1")
  geneB <- if ("gene_id2" %in% names(dt) && any(dt$gene_id2 != "." & !is.na(dt$gene_id2))) get_char("gene_id2") else get_char("gene2")
  geneA[geneA == "."] <- NA_character_
  geneB[geneB == "."] <- NA_character_
  
  pair <- paste0(ifelse(is.na(geneA), "NA", geneA), "--", ifelse(is.na(geneB), "NA", geneB))
  te_involved <- grepl("^TE_", geneA) | grepl("^TE_", geneB)
  
  bp1 <- if ("breakpoint1" %in% names(dt)) as.character(dt$breakpoint1) else NA_character_
  bp2 <- if ("breakpoint2" %in% names(dt)) as.character(dt$breakpoint2) else NA_character_
  
  data.table(
    sample = sample, depth = depth, rep = rep, mode = mode, status = status,
    geneA = geneA, geneB = geneB, pair = pair, te_involved = te_involved,
    breakpoint1 = bp1, breakpoint2 = bp2,
    type = get_char("type"),
    confidence = get_char("confidence"),
    split_reads1 = get_int("split_reads1"),
    split_reads2 = get_int("split_reads2"),
    discordant_mates = get_int("discordant_mates"),
    coverage1 = get_int("coverage1"),
    coverage2 = get_int("coverage2"),
    filters = get_char("filters")
  )
}

arriba_all <- rbindlist(lapply(files, read_arriba), use.names = TRUE, fill = TRUE)
stopifnot(nrow(arriba_all) > 0)

# Save raw merged table (transparent)
saveRDS(arriba_all, file = file.path(merge_dir, "arriba_all.rds"))

# ------------------------------------------------------------------------------
# 3) Choose calls to evaluate: kept only; prefer relaxed if present
# ------------------------------------------------------------------------------

DT <- copy(arriba_all)
DT <- DT[status == "kept"]
if (any(DT$mode == "relaxed", na.rm = TRUE)) {
  DT <- DT[mode == "relaxed"]
} else {
  DT <- DT[mode == "strict"]
}
stopifnot(nrow(DT) > 0)

# ------------------------------------------------------------------------------
# 4) Define support + event_id (event-level unit)
# ------------------------------------------------------------------------------

DT[, support := fifelse(is.na(split_reads1), 0L, split_reads1) +
     fifelse(is.na(split_reads2), 0L, split_reads2) +
     fifelse(is.na(discordant_mates), 0L, discordant_mates)]

DT[, event_id := paste0(sample, ":", sprintf("%06d", .I))]

# ------------------------------------------------------------------------------
# 5) Parse breakpoints -> GRanges
# ------------------------------------------------------------------------------

parse_bp <- function(x) {
  m <- regexec("^([^:]+):(\\d+)$", x)
  r <- regmatches(x, m)
  chr <- vapply(r, function(z) if (length(z) >= 2) z[2] else NA_character_, "")
  pos <- suppressWarnings(as.integer(vapply(r, function(z) if (length(z) >= 3) z[3] else NA_character_, "")))
  list(chr = chr, pos = pos)
}

b1 <- parse_bp(DT$breakpoint1)
b2 <- parse_bp(DT$breakpoint2)

DT[, `:=`(bp1_chr=b1$chr, bp1_pos=b1$pos, bp2_chr=b2$chr, bp2_pos=b2$pos)]
DT <- DT[!is.na(bp1_chr) & !is.na(bp1_pos) & !is.na(bp2_chr) & !is.na(bp2_pos)]
stopifnot(nrow(DT) > 0)

grA <- GRanges(seqnames = DT$bp1_chr, ranges = IRanges(start = DT$bp1_pos, end = DT$bp1_pos))
grB <- GRanges(seqnames = DT$bp2_chr, ranges = IRanges(start = DT$bp2_pos, end = DT$bp2_pos))

# ------------------------------------------------------------------------------
# 6) Truth transcript intervals + TruthN (same as your history logic)
# ------------------------------------------------------------------------------

truth_gtf <- import(truth_gtf_path)

truth_tx <- truth_gtf[mcols(truth_gtf)$type %in% "transcript"]
if (length(truth_tx) == 0L) {
  truth_tx <- unlist(reduce(split(truth_gtf[mcols(truth_gtf)$type %in% "exon"],
                                  mcols(truth_gtf)$transcript_id)))
  mcols(truth_tx)$transcript_id <- names(split(truth_gtf[mcols(truth_gtf)$type %in% "exon"],
                                               mcols(truth_gtf)$transcript_id))
}

truth_txids_all <- as.character(mcols(truth_tx)$transcript_id)
TruthN <- uniqueN(truth_txids_all)

# ------------------------------------------------------------------------------
# 7) Overlap: breakpoint hits truth transcript? (event-level TP)
# ------------------------------------------------------------------------------

ovA <- findOverlaps(grA, truth_tx, ignore.strand = TRUE, minoverlap = 1L)
ovB <- findOverlaps(grB, truth_tx, ignore.strand = TRUE, minoverlap = 1L)

DT[, idx := .I]
DT[, txA := vector("list", .N)]
DT[, txB := vector("list", .N)]

if (length(ovA)) {
  txA_list <- split(truth_txids_all[subjectHits(ovA)], queryHits(ovA))
  DT[as.integer(names(txA_list)), txA := unname(txA_list)]
}
if (length(ovB)) {
  txB_list <- split(truth_txids_all[subjectHits(ovB)], queryHits(ovB))
  DT[as.integer(names(txB_list)), txB := unname(txB_list)]
}

DT[, truth_union := Map(union, txA, txB)]
DT[, isTP := lengths(truth_union) > 0L]

# ------------------------------------------------------------------------------
# 8) Mild filtering + summarize by depth
# ------------------------------------------------------------------------------

DT[, keep := (support >= 2L) &
     !grepl("read_through", tolower(ifelse(is.na(filters), "", filters)))]

sum_before <- DT[, .(
  Pred_all = uniqueN(event_id),
  TP_all   = uniqueN(event_id[isTP])
), by = depth]

sum_after <- DT[keep == TRUE, .(
  Pred_kept  = uniqueN(event_id),
  TP_kept    = uniqueN(event_id[isTP]),
  HitTruthTx = uniqueN(na.omit(unlist(truth_union)))
), by = depth]

perf_filtered <- merge(sum_before, sum_after, by = "depth", all = TRUE)[
  , `:=`(
    FP_kept   = pmax(Pred_kept - TP_kept, 0L),
    FN        = pmax(TruthN - HitTruthTx, 0L),
    Precision = ifelse(is.na(Pred_kept) | Pred_kept == 0, NA_real_, TP_kept / Pred_kept),
    Recall    = HitTruthTx / TruthN
  )
][, F1 := ifelse(is.na(Precision) | is.na(Recall) | (Precision + Recall) == 0,
                 NA_real_, 2 * Precision * Recall / (Precision + Recall))
][order(depth)]

print(perf_filtered)

saveRDS(perf_filtered, file = file.path(merge_dir, "perf_filtered.rds"))


# ------------------------------------------------------------------------------
# 9) Convenience: tool-aligned benchmark table
# ------------------------------------------------------------------------------

arriba_benchmark <- perf_filtered[, .(
  Depth = depth,
  Tool  = "Arriba",
  TP = TP_kept,
  FP = FP_kept,
  FN = FN,
  Precision = Precision,
  Recall = Recall,
  F1 = F1
)]