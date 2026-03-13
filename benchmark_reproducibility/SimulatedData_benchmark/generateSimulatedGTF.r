suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(parallel)
})

# ==============================================================================
# Script: generate_simulated_gtf.R
# Purpose: Generate a ground-truth GTF by stochastically sampling TE-associated
#          transcript models from a reconstructed transcript annotation.
# Notes:
#   - The baseline annotation should be a reconstructed/curated transcriptome
#     background derived from the biological system of interest.
#   - Transcript simulation starts from exon structures already observed in the
#     baseline annotation.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. User parameters
# ------------------------------------------------------------------------------
baseline_gtf_path <- "input/reconstructed_early_embryo_transcriptome.gtf"
te_annotation_rds <- "input/mouse_te_annotation.rds"
output_gtf_path   <- "output/simulated_truth_1000_transcripts.gtf"

n_genes_to_sample <- 1000L
min_te_overlap    <- 15L
min_downstream_if_sparse <- 1L
max_downstream_if_sparse <- 2L
min_downstream_if_rich   <- 3L
max_downstream_if_rich   <- 5L

set.seed(1)

# ------------------------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------------------------
message("Loading baseline transcript annotation and TE annotation...")
gtf_all <- import(baseline_gtf_path)
te      <- readRDS(te_annotation_rds)

# ------------------------------------------------------------------------------
# 3. Keep exon features and assign exon ranks
# ------------------------------------------------------------------------------
message("Preparing exon-level annotation...")
gtf_exon <- gtf_all[gtf_all$type == "exon"]
gtf_exon <- gtf_exon[strand(gtf_exon) != "*"]

gtf_df <- as.data.frame(gtf_exon)
gtf_df$row_id <- seq_len(nrow(gtf_df))

gtf_df <- gtf_df %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  mutate(exon_number = row_number()) %>%
  ungroup()

gtf_df <- gtf_df %>%
  group_by(transcript_id) %>%
  mutate(tx_exon_rank = if (unique(strand) == "-") rev(exon_number) else exon_number) %>%
  ungroup() %>%
  arrange(row_id)

gtf_exon$exon_number <- gtf_df$exon_number
gtf_exon$tx_exon_rank <- gtf_df$tx_exon_rank

# ------------------------------------------------------------------------------
# 4. Annotate TE-overlapping exons
# ------------------------------------------------------------------------------
message("Annotating TE-overlapping exons...")
hits <- findOverlaps(gtf_exon, te, minoverlap = min_te_overlap)

gtf_exon$TE_name  <- "none"
gtf_exon$TE_class <- "none"

if (length(hits) > 0) {
  te_name_list <- split(te$name[subjectHits(hits)], queryHits(hits))
  te_name_list <- lapply(te_name_list, function(x) paste(unique(x), collapse = ","))
  
  te_class_list <- split(te$class[subjectHits(hits)], queryHits(hits))
  te_class_list <- lapply(te_class_list, function(x) paste(unique(x), collapse = ","))
  
  gtf_exon$TE_name[as.integer(names(te_name_list))]  <- unlist(te_name_list)
  gtf_exon$TE_class[as.integer(names(te_class_list))] <- unlist(te_class_list)
}

# ------------------------------------------------------------------------------
# 5. Define exon pools
# ------------------------------------------------------------------------------
first_exons <- gtf_exon[gtf_exon$tx_exon_rank == 1]
other_exons <- gtf_exon[gtf_exon$tx_exon_rank != 1]

# Optional: restrict to genes that have at least one TE-related exon anywhere
te_gene_ids <- unique(gtf_exon$gene_id[gtf_exon$TE_name != "none"])
te_gene_ids <- te_gene_ids[!is.na(te_gene_ids)]

candidate_gene_ids <- intersect(unique(first_exons$gene_id), te_gene_ids)
candidate_gene_ids <- candidate_gene_ids[!is.na(candidate_gene_ids)]

if (length(candidate_gene_ids) < n_genes_to_sample) {
  stop("Not enough candidate genes to sample the requested number of genes.")
}

# ------------------------------------------------------------------------------
# 6. Simulation helper
# ------------------------------------------------------------------------------
simulate_one_gene <- function(gene_id, first_exons, other_exons, transcript_id_prefix) {
  gene_first_exons <- first_exons[first_exons$gene_id == gene_id]
  if (length(gene_first_exons) == 0) return(NULL)
  
  initiated_exon <- gene_first_exons[sample(seq_along(gene_first_exons), 1)]
  gene_strand <- as.character(strand(initiated_exon))
  
  gene_other_exons <- other_exons[other_exons$gene_id == gene_id]
  if (length(gene_other_exons) == 0) return(NULL)
  
  if (gene_strand == "+") {
    valid_other_exons <- gene_other_exons[gene_other_exons$exon_number > initiated_exon$exon_number]
  } else {
    valid_other_exons <- gene_other_exons[gene_other_exons$exon_number < initiated_exon$exon_number]
  }
  
  if (length(valid_other_exons) == 0) {
    selected_exons <- initiated_exon
  } else {
    n_downstream <- if (length(valid_other_exons) < 5) {
      sample(min_downstream_if_sparse:max_downstream_if_sparse, 1)
    } else {
      sample(min_downstream_if_rich:max_downstream_if_rich, 1)
    }
    
    n_downstream <- min(n_downstream, length(valid_other_exons))
    sampled_idx <- sample(seq_len(length(valid_other_exons)), n_downstream)
    sampled_other_exons <- valid_other_exons[sampled_idx]
    
    exon_order <- order(start(c(initiated_exon, sampled_other_exons)))
    selected_exons <- c(initiated_exon, sampled_other_exons)[exon_order]
  }
  
  new_tx_id <- paste0(transcript_id_prefix, gene_id)
  
  exon_out <- selected_exons
  exon_out$type <- "exon"
  exon_out$transcript_id <- new_tx_id
  exon_out$gene_id <- gene_id
  
  tx_range <- range(exon_out)
  tx_out <- GRanges(
    seqnames = seqnames(tx_range),
    ranges = ranges(tx_range),
    strand = strand(tx_range)
  )
  tx_out$type <- "transcript"
  tx_out$gene_id <- gene_id
  tx_out$gene_name <- exon_out$gene_name[1]
  tx_out$transcript_id <- new_tx_id
  
  c(tx_out, exon_out)
}

# ------------------------------------------------------------------------------
# 7. Sample genes and generate truth GTF
# ------------------------------------------------------------------------------
message("Sampling genes and simulating transcript models...")
sampled_gene_ids <- sample(candidate_gene_ids, n_genes_to_sample, replace = FALSE)

simulated_list <- lapply(
  sampled_gene_ids,
  simulate_one_gene,
  first_exons = first_exons,
  other_exons = other_exons,
  transcript_id_prefix = "SIMTX_"
)

simulated_list <- Filter(Negate(is.null), simulated_list)

if (length(simulated_list) == 0) {
  stop("No simulated transcripts were generated.")
}

simulated_gtf <- sort(do.call(c, simulated_list))

# ------------------------------------------------------------------------------
# 8. Export
# ------------------------------------------------------------------------------
dir.create(dirname(output_gtf_path), recursive = TRUE, showWarnings = FALSE)
export(simulated_gtf, output_gtf_path)
message("Done: ", output_gtf_path)