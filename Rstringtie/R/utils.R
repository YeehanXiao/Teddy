.Stringtiebin <- function(args = "") {
  if (is.null(args) || args == "") {
    stop("The stringtie binaries nedd to called",
         " with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  args <- unlist(strsplit(args, split = " "))
  bin <- file.path(system.file(package = "Rstringtie"), "stringtie")
  output <- system2(command = bin, args = args)
}

.gffcompareBin <- function(args = "") {
  if (is.null(args) || args == "") {
    stop("The stringtie binaries nedd to called",
         " with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  args <- unlist(strsplit(args, split = " "))
  bin <- file.path(system.file(package = "Rstringtie"), "gffcompare")
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

