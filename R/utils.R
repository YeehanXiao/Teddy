#' Standardize transposon annotations for TEDDY
#'
#' This helper prepares a transposon annotation object for TEDDY by converting
#' chromosome naming to NCBI style and ensuring that a TE-name metadata column
#' named `names` is available. If `names` is absent, common TE-name columns such
#' as `name`, `repName`, `repname`, or `TE_name` are used to create it.
#'
#' @param transposon A \link[GenomicRanges:GRanges-class]{GRanges} object
#'   containing transposon annotations.
#' @param replace_name Logical. If `TRUE` and the object contains a `name`
#'   column but not a `names` column, the `name` column is renamed to `names`.
#'   If `FALSE`, a new `names` column is created while preserving the original
#'   `name` column. Default is `FALSE`.
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object with NCBI-style
#'   chromosome names and a character metadata column named `names`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' transposon <- NCBI_check(transposon)
#' transposon <- NCBI_check(transposon, replace_name = TRUE)
#' }
NCBI_check <- function(transposon, replace_name = FALSE) {
  if (!inherits(transposon, "GRanges")) {
    stop("`transposon` must be a GRanges object.")
  }
  
  mcols_names <- colnames(S4Vectors::mcols(transposon))
  
  if (!"names" %in% mcols_names) {
    candidate_cols <- c("name", "repName", "repname", "TE_name")
    hit <- candidate_cols[candidate_cols %in% mcols_names][1]
    
    if (!is.na(hit)) {
      if (replace_name && hit == "name") {
        colnames(S4Vectors::mcols(transposon))[
          colnames(S4Vectors::mcols(transposon)) == "name"
        ] <- "names"
      } else {
        transposon$names <- as.character(S4Vectors::mcols(transposon)[[hit]])
      }
    } else {
      stop("Required TE name column missing. Please provide one of: `names`, `name`, `repName`, `repname`, or `TE_name`.")
    }
  } else {
    transposon$names <- as.character(transposon$names)
  }
  
  GenomeInfoDb::seqlevelsStyle(transposon) <- "NCBI"
  GenomeInfoDb::seqlevels(transposon) <- sub("^chr", "", GenomeInfoDb::seqlevels(transposon))
  GenomeInfoDb::seqnames(transposon) <- sub("^chr", "", GenomeInfoDb::seqnames(transposon))
  
  transposon
}


.Stringtiebin <- function(args = "") {
  if (is.null(args) || args == "") {
    stop("The StringTie executables require additional arguments.")
  }
  args <- gsub("^ *| *$", "", args)
  args <- unlist(strsplit(args, split = " "))
  bin <- file.path(system.file(package = "Teddy"), "stringtie")
  output <- system2(command = bin, args = args)
}


.Stringtiebin3 <- function(args = "") {
  if (is.null(args) || args == "") {
    stop("The StringTie executables require additional arguments.")
  }
  args <- gsub("^ *| *$", "", args)
  args <- unlist(strsplit(args, split = " "))
  bin <- file.path(system.file(package = "Teddy"), "stringtie3", "stringtie")
  output <- system2(command = bin, args = args)
  return(output)
}

.gffcompareBin <- function(args = "") {
  if (is.null(args) || args == "") {
    stop("Gffcompare requires additional arguments.")
  }
  args <- gsub("^ *| *$", "", args)
  args <- unlist(strsplit(args, split = " "))
  bin <- file.path(system.file(package = "Teddy"), "gffcompare")
  output <- system2(command = bin, args = args)
}



.parseDots <- function(...) {
  if (...length() == 0) return("")
  dots <- list(...)
  args <- lapply(names(dots), FUN = function(x) {
    paste(x, dots[x], sep = " ")
  })
  args <- paste0(unlist(args), collapse = " ")
  return(args)
}

.ExtractTranscript <- function(gtf, transcriptGR) {
  gtfs <- rtracklayer::import.gff(gtf)
  
  transcript <- gtfs[gtfs$type == "transcript"]
  exprts <- GenomicRanges::mcols(transcript)[, c("cov", "FPKM", "TPM")]
  exprts <- apply(X = exprts, MARGIN = 2, as.numeric)
  exprts <- as.data.frame(exprts, row.names = transcript$transcript_id)
  
  express <- matrix(0, nrow = length(transcriptGR), ncol = 3,
                    dimnames = list(transcriptGR$transcript_id,
                                    c("cov", "FPKM", "TPM")))
  express <- as.data.frame(express)
  
  express[rownames(exprts), ] <- exprts
  express <- express[transcriptGR$transcript_id, ]
  
  gtf <- basename(gtf)
  sample <- gsub(pattern = ".gtf$", replacement = "", x = gtf)
  coldata <- data.frame(row.names = sample, sample = sample)
  
  
  cov <- express[, 'cov', drop = FALSE]
  colnames(cov) <- sample
  FPKM = express[, 'FPKM', drop = FALSE]
  colnames(FPKM) <- sample
  TPM = express[, "TPM", drop = FALSE]
  colnames(TPM) <- sample
  
  SE <- SummarizedExperiment::SummarizedExperiment(
    S4Vectors::SimpleList(cov = cov, FPKM = FPKM,
                          TPM = TPM),
                          rowRanges = transcriptGR,
                          colData = coldata)
  return(SE)
}



UnfoldCoefs <- function(varlist, fit, mm, mf){
  coefLists <- lapply(seq(ncol(varlist)), function(index) {
    termName <- colnames(varlist)[index]
    varsInTerm <- stringr::str_split(termName, stringr::fixed(":"))[[1]]
    stopifnot(all(varlist[varsInTerm, index] == 1))
    stopifnot(sum(varlist[ , index]) == length(varsInTerm))
    coefNames <- colnames(mm)[index+1]
    varLevels <- lapply(varsInTerm, function(v) levels(factor(mf[[v]])))
    coefIndices <- array(0, dim = sapply(varLevels, length), dimnames = varLevels)
    lvlTbl <- stringr::str_match(coefNames, stringr::str_c( "^", stringr::str_c( varsInTerm, "([^:]*)", collapse=":" ), "$" ) )[ , -1, drop=FALSE]
    coefIndices[lvlTbl] <- coefficients(fit)[coefNames]
    coefIndices
  })
  names(coefLists) <- colnames(varlist)
  a <- array(c(`(Intercept)` = "(Intercept)"))
  coefLists <- c(list(a), coefLists)
  names(coefLists)[1] <- "intercept"
  coefLists[1] <- coefficients(fit)[1]
  return(coefLists)
}

rmDepCols <- function(m) {
  q <- qr(m)
  if (q$rank < ncol(m)) 
    m[, -q$pivot[(q$rank + 1):ncol(m)]]
  else m
}

getEffectsForGene <- function(gene,
                              groups,
                              findex,
                              frm,
                              mf,
                              disp,
                              otherCounts,
                              chimericCounts,
                              N_samples,
                              samples,
                              features,
                              object){
  idx <- groups %in% gene & findex
  sub_chimericCounts <- chimericCounts[idx, ,drop = FALSE]
  sub_otherCounts <- otherCounts[idx, ,drop = FALSE]
  sub_disp <- disp[idx]
  names(sub_disp) <- features[idx]
  N_exons <- sum(idx)
  newMf <- as.data.frame(rep(samples, each = N_exons))
  colnames(newMf) <- "sample"
  chi_sizeFactors <- mf$sizeFactor[colData(object)$chimeric == "chimeric"]
  chi_condition <- mf$condition[colData(object)$chimeric == "chimeric"]
  newMf$count <- as.vector(t(sub_chimericCounts[, samples]))
  newMf$dispersion <- rep(sub_disp, N_samples)
  newMf$sizeFactors <- rep(chi_sizeFactors, each = N_exons)
  newMf$Chimeric <- "chimeric"
  newMf$condition <- rep(chi_condition, each = N_exons)
  othersMf <- newMf
  othersMf$Chimeric <- "others"
  othersMf$count <- as.vector(t(sub_otherCounts[, samples]))
  Mf <- rbind(othersMf, newMf)[,c("sample", "condition", "Chimeric", "sizeFactors", "dispersion", "count")]
  Mf$Chimeric <- factor(Mf$Chimeric, levels = c("others", "chimeric"))
  mm <- model.matrix(frm, Mf)
  fit <- glmnb.fit(mm, Mf$count, dispersion = Mf$dispersion, offset = log(Mf$sizeFactors), tol=0.01)
  varlist <- attr(terms(frm), "factors")
  coef <- UnfoldCoefs(varlist = varlist, fit = fit, mm = mm, mf = Mf)
  return(coef)
}


