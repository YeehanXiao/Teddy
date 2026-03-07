# ==============================================================================
# Script: generateSimulatedGTF.R
# Purpose: Generate the ground truth GTF by stochastically sampling 1,000 
#          TE-associated transcripts from a baseline reference annotation.
# ==============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 1. Path Configuration & Data Loading
# ------------------------------------------------------------------------------
# Define paths (Adjust input_gtf_path to your actual baseline reference annotation)
input_gtf_path <- "/mnt/datadisk/xiaoyihan/index/mm10_reference.gtf" 
te_path        <- "/mnt/datadisk/xiaoyihan/manuals/compare/simulate/meta/mouse_te.rds"
output_gtf     <- "1000_simulated_transcripts.gtf"

message("Loading reference GTF and TE annotations...")
GTF <- import(input_gtf_path)
te  <- readRDS(te_path)

# ------------------------------------------------------------------------------
# 2. GTF Pre-processing & Exon Ranking
# ------------------------------------------------------------------------------
message("Processing exons and calculating strand-aware exon ranks...")

# Filter for valid exons and remove unstranded entries
GTF <- GTF[GTF$type == "exon"]
GTF <- GTF[strand(GTF) != "*"]

# Convert to data frame for dplyr operations
GTF_df <- as.data.frame(GTF)
GTF_df$Number <- seq_len(nrow(GTF_df))

# Rank exons by start position (5' to 3' absolute direction)
GTF_df <- GTF_df %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  mutate(exon_number = row_number()) %>%
  ungroup() %>%
  arrange(Number)

GTF$exon_number <- GTF_df$exon_number

# Calculate actual transcriptional rank based on strand direction
GTF_df <- GTF_df %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  mutate(tx_exon_rank = case_when(
    unique(strand) == "-" ~ rev(exon_number),
    TRUE ~ exon_number
  )) %>%
  ungroup() %>%
  arrange(Number)

GTF$tx_exon_rank <- GTF_df$tx_exon_rank

# ------------------------------------------------------------------------------
# 3. TE Overlap Annotation
# ------------------------------------------------------------------------------
message("Annotating exons with Transposable Element (TE) overlaps...")

overlaps <- findOverlaps(GTF, te, minoverlap = 15)

# Extract class and repeat names for overlapping features
class_list <- split(te$class[subjectHits(overlaps)], queryHits(overlaps))
class_list <- lapply(class_list, function(x) paste(unique(x), collapse = ","))

repeat_list <- split(te$name[subjectHits(overlaps)], queryHits(overlaps))
repeat_list <- lapply(repeat_list, function(x) paste(unique(x), collapse = ","))

# Initialize and assign TE info
GTF$TE_name  <- "none"
GTF$TE_class <- "none"

GTF$TE_name[as.numeric(names(repeat_list))]  <- unlist(repeat_list)
GTF$TE_class[as.numeric(names(class_list))] <- unlist(class_list)

# Split exons by their relative position within the transcript
first_exons <- GTF[GTF$tx_exon_rank == 1]
other_exons <- GTF[GTF$tx_exon_rank != 1]

# ------------------------------------------------------------------------------
# 4. Stochastic Transcript Simulation Logic
# ------------------------------------------------------------------------------
simulate_for_gene <- function(gene_id, first_exons, other_exons, named_transcriptID) {
  
  # Sample one initiating exon
  gene_first_exons <- first_exons[first_exons$gene_id == gene_id, ]
  initiated_exon   <- gene_first_exons[sample(seq_along(gene_first_exons), 1)]
  gene_strand      <- as.character(strand(initiated_exon))
  
  # Filter valid downstream exons strictly based on strand progression
  gene_other_exons <- other_exons[other_exons$gene_id == gene_id, ]
  if (gene_strand == "+") {
    valid_other_exons <- gene_other_exons[gene_other_exons$exon_number > initiated_exon$exon_number, ]
  } else {
    valid_other_exons <- gene_other_exons[gene_other_exons$exon_number < initiated_exon$exon_number, ]
  }
  
  # Determine number of downstream exons to append (stochastic length)
  num_other_exons <- ifelse(length(valid_other_exons) < 5, sample(2:3, 1), sample(3:5, 1))
  num_other_exons <- min(num_other_exons, length(valid_other_exons))
  
  # Sample and reduce downstream exons
  sampled_idx <- sample(seq_len(length(valid_other_exons)), num_other_exons)
  sampled_other_exons <- reduce(valid_other_exons[sampled_idx])
  
  # Combine and format metadata
  object <- reduce(sort(c(initiated_exon, sampled_other_exons))[, c("type", "gene_id", "gene_name", "transcript_id")])
  object$transcript_id <- named_transcriptID
  object$type          <- "exon"
  object$gene_id       <- gene_id
  object$gene_name     <- initiated_exon$gene_name
  
  # Create parent transcript feature
  transcript <- GRanges(
    seqnames      = seqnames(range(object)),
    ranges        = ranges(range(object)),
    strand        = strand(range(object)),
    type          = "transcript",
    gene_id       = gene_id,
    transcript_id = named_transcriptID,
    gene_name     = initiated_exon$gene_name
  )
  
  return(c(transcript, object))
}

# ------------------------------------------------------------------------------
# 5. Execute Simulation & Export
# ------------------------------------------------------------------------------
message("Simulating 1,000 unique transcripts...")

# Sample exactly 1,000 genes to act as the ground truth framework
available_gene_ids <- unique(first_exons$gene_id)
sampled_gene_ids   <- sample(available_gene_ids, 1000)

simulate_data <- lapply(seq_along(sampled_gene_ids), function(index) {
  named_transcriptID <- paste0("MSTRG.", index)
  simulate_for_gene(
    gene_id = sampled_gene_ids[index],
    first_exons = first_exons,
    other_exons = other_exons,
    named_transcriptID = named_transcriptID
  )
})

combined_data <- sort(do.call(c, simulate_data))

export(combined_data, output_gtf)

