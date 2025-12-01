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
                        minoverlap =15,
                        genome = BSgenome.Mmusculus.UCSC.mm10,
                        min.score = 0.9,
                        filter = 0) {
  overlaps <- findOverlaps(object, te, minoverlap = minoverlap)
  target_te <- te[subjectHits(overlaps)]
  results <- lapply(
    getSeq(x = BSgenome.Mmusculus.UCSC.mm10, target_te), function(seq) matchPWM(pwm, seq, min.score = min.score))
  lengths_of_results <- sapply(results, function(result) length(result@ranges))
  reattach_overlap <- findOverlaps(
    object, 
    target_te[which(lengths_of_results > filter)],
    minoverlap = minoverlap
  )
  reattach <- object[queryHits(reattach_overlap)]
  return(list(match_results = results, overlap_results = reattach,target_te[which(lengths_of_results > filter)]))
}





