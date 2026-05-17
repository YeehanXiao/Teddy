#' @title Convert Probability Matrix to Position Weight Matrix
#' @param pcm A numeric matrix with 4 rows representing nucleotide probabilities for 'A', 'C', 'G', and 'T'.
#' Each column represents a position in the motif. The values should be proportions and each column should sum up to 1.
#' @return A PWM object that can be used in further analysis with motif matching functions. 
#' @export
pcmFunction <- function(pcm){
  pcm <- t(pcm)
  rownames(pcm) <- c("A","C","G","T")
  pcm <- sweep(pcm, 2, colSums(pcm), FUN="/")
  t <- matrix(as.integer(pcm * 1000), nrow = 4)
  rownames(t) <- c("A","C","G","T")
  adjustment <- function(col) {
    diff <- 1000 - sum(col)
    if (diff != 0) {
      col[which.max(col)] <- col[which.max(col)] + diff
    }
    return(col)
  }
  t <- apply(t, 2, adjustment)
  t <- matrix(as.integer(t), nrow = 4)
  rownames(t) <- c("A","C","G","T")
  pwm <- PWM(t)
  return(pwm)
}

#' @title Search the motif binding sites through TE-chimeric transcripts
#' @param object A Granges object of target TE-chimeric transcripts.
#' @param te A Granges object of TE location annotations, typically derived from annotated TE reference.
#' @param pwm The pwm matrix for the target motif.
#' @param minoverlap A non-negative integer specifying the minimum number for ranges to be considered overlapping. Default is 15.
#' @param genome The genome to use, from the BSgenome package (e.g., BSgenome.Mmusculus.UCSC.mm10).
#' @param min.score A numeric value specifying the minimum score for a PWM match to be considered significant. Default is 0.9.
#' @param filter A non-negative integer specifying the threshold of matches present for a TE-chimeric transcripts. This allows for filtering out transcripts with a low number of significant motif matches, focusing the analysis on more likely candidates. Default is 0.
#' @importFrom GenomicRanges findOverlaps
#' @importFrom Biostrings getSeq matchPWM
#' @import BSgenome.Mmusculus.UCSC.mm10
#' @export
MotifSearch <- function(object,
                        te,
                        pwm,
                        minoverlap = 15,
                        genome,
                        min.score = 0.9,
                        filter = 0) {
  if (missing(genome)) {
    stop("Please provide a genome object matching the annotation build.")
  }
  overlaps <- GenomicRanges::findOverlaps(
    object,
    te,
    minoverlap = minoverlap
  )
  target_te <- te[S4Vectors::subjectHits(overlaps)]
  match_results <- lapply(
    Biostrings::getSeq(x = genome, target_te),
    function(seq) {
      Biostrings::matchPWM(
        pwm,
        seq,
        min.score = min.score
      )
    }
  )
  n_matches <- vapply(
    match_results,
    function(x) length(x@ranges),
    numeric(1)
  )
  keep <- which(n_matches > filter)
  reattach_overlap <- GenomicRanges::findOverlaps(
    object,
    target_te[keep],
    minoverlap = minoverlap
  )
  reattach <- object[S4Vectors::queryHits(reattach_overlap)]
  return(list(
    match_results = match_results,
    overlap_results = reattach,
    target_te = target_te[keep]
  ))
}

#' Export transcript sequences
#'
#' Export spliced transcript sequences by retrieving exon sequences from an
#' exon-level GTF `GRanges` object. Exons can either be concatenated into full
#' transcript sequences or exported separately.
#'
#' @param txid Character vector. Transcript IDs to export.
#' @param GTF A `GRanges` object containing exon-level transcript annotation,
#' such as the output of `extractGTF(type = "exon")`.
#' @param genome Genome sequence object accepted by `Biostrings::getSeq()`,
#' such as a `BSgenome` object or `FaFile`.
#' @param outdir Character. Output directory for FASTA files.
#' Default is `"transcript_fasta"`.
#' @param write Logical. If `TRUE`, write sequences to FASTA files and invisibly
#' return output file paths. If `FALSE`, return a `DNAStringSet`.
#' @param collapse Logical. If `TRUE`, concatenate exon sequences into one
#' spliced transcript sequence. If `FALSE`, return or write exon sequences
#' separately.
#'
#' @return If `write = TRUE`, invisibly returns output FASTA file paths.
#' If `write = FALSE`, returns a `DNAStringSet`.
#'
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges strand start
#' @importFrom Biostrings getSeq xscat DNAStringSet writeXStringSet reverseComplement
#'
#' @export
#'
export_transcript <- function(txid,
                              GTF,
                              genome,
                              outdir = "transcript_fasta",
                              write = TRUE,
                              collapse = TRUE) {
  if (!inherits(GTF, "GRanges")) {
    stop("GTF must be a GRanges object.")
  }
  
  if (!"transcript_id" %in% colnames(S4Vectors::mcols(GTF))) {
    stop("GTF must contain a transcript_id column.")
  }
  
  export_one <- function(one_txid) {
    gr <- GTF[GTF$transcript_id == one_txid]
    
    if (length(gr) == 0) {
      message("Transcript ", one_txid, " not found in GTF; skipped.")
      return(NULL)
    }
    
    tx_strand <- unique(as.character(GenomicRanges::strand(gr)))
    tx_strand <- tx_strand[tx_strand != "*"]
    
    if (length(tx_strand) == 1 && tx_strand == "-") {
      gr <- gr[order(GenomicRanges::start(gr), decreasing = TRUE)]
    } else {
      gr <- gr[order(GenomicRanges::start(gr))]
    }
    
    seqs <- Biostrings::getSeq(genome, gr)
    
    if (length(tx_strand) == 1 && tx_strand == "-") {
      seqs <- Biostrings::reverseComplement(seqs)
    }
    
    if (isTRUE(collapse)) {
      seq_concat <- do.call(Biostrings::xscat, as.list(seqs))
      seq_out <- Biostrings::DNAStringSet(seq_concat)
      names(seq_out) <- one_txid
    } else {
      seq_out <- Biostrings::DNAStringSet(seqs)
      names(seq_out) <- paste0("exon", seq_along(seq_out))
    }
    
    if (isTRUE(write)) {
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
      outfile <- file.path(outdir, paste0(one_txid, ".fa"))
      Biostrings::writeXStringSet(seq_out, filepath = outfile)
      return(outfile)
    }
    
    seq_out
  }
  
  res <- lapply(txid, export_one)
  res <- Filter(Negate(is.null), res)
  
  if (isTRUE(write)) {
    return(invisible(unlist(res)))
  }
  
  do.call(c, res)
}


