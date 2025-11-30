#' @title R wrapper to Run stringtie
#' @description Run stringtie for transcriptome assembly 
#' @param bam Character, BAM file (sorted by genomic coordinates)
#' @param reference Character, reference annotation GTF (optional)
#' @param outfile Character, output GTF file name
#' @param params Character, extra parameters (default "")
#' @param longRead Logical, whether the data is long-read (PacBio/ONT). Default FALSE (short-read).

#'
#' @export
stringtieAssembly <- function(bam, reference = NULL, outfile, params = "", longRead = FALSE) {
  if (!all(file.exists(bam))) stop("BAM file not found: ", bam)
  
  args <- c(bam)
  
  if (longRead) {
    args <- c(args, "-L")
  }
  

  if (!is.null(reference)) {
    args <- c(args, "-G", reference)
  }
  
  args <- c(args, "-o", outfile)
  

  if (nzchar(params)) {
    args <- c(args, strsplit(params, "\\s+")[[1]])
  }

  if (longRead) {
    rc <- .Stringtiebin3(paste(args, collapse = " "))
  } else {
    rc <- .Stringtiebin(paste(args, collapse = " "))
  }
  
  invisible(rc)
}

#' @title R wrapper to Run stringtie tool
#' @description Merge transcripts
#'
#' @param reference Use a reference annotation file to guide assembly process
#' @param gtfFiles Character, GTF files assembled from transcripts by stringtie
#' @param outfile Character, the output of the merged GTF
#' @param params Other parameters
#' @export
stringtieMerge <- function(reference, gtfFiles, outfile, params = "") {
  programtags <- "--merge"
  reference <- paste("-G", reference, sep = " ")
  gtffile <- paste(gtfFiles, collapse = " ")
  outfile <- paste("-o", outfile, sep = " ")
  cmd <- sprintf("%s %s %s %s %s",
                 programtags,
                 reference,
                 gtffile,
                 outfile,
                 params)
  return(invisible(lapply(cmd, .Stringtiebin)))
}

#' @title R wrapper to Run gffcompare
#' @description The function to compare and merge accuracy of one or more GFF files (the “query” files),
#' when compared with a reference annotation (also provided as GFF).
#' @param reference Use a reference annotation file to guide compare assembly process.
#' @param gtffile GTF files with gffcompare annotation transcripts.
#' @param outfile The name of the output annotated merged GTF.
#' @param params Other parameters
#' @export

gffcompareAnno <- function(reference, gtffile, outfile, params = "") {
  reference <- paste("-r", reference, sep = " ")
  outfile <- gsub(pattern = "[.]gtf$", replacement = "", x = outfile)
  gtffile <- paste(gtffile, collapse = " ")
  outfile <- paste("-o", outfile, sep = " ")
  cmd <- sprintf("%s %s %s %s",
                 reference,
                 outfile,
                 gtffile,
                 params)
  return(invisible(lapply(cmd, .gffcompareBin)))
}

#' @title Preparing the genome annotation object
#' @description Flatten exon appearing multiple times among different transcripts in GTF file
#' @param gtffile GTF file.
#' @param singleGens Whether to allocate the exon overlapping with two genes to a single gene. Default is TRUE.
#' @param transposon A GRanges object with transposon data
#' @param minoverlap Minimum overlap for \code{\link[IRanges]{findOverlaps}}. Default is 10.
#' @param cores Number of cores used for parallel processing. Default is 1.
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom rtracklayer import.gff
#' @importFrom GenomicFeatures exonicParts makeTxDbFromGRanges
#' @importFrom GenomicRanges strand
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors subjectHits queryHits mcols
#' @export
prepareAnno <- function(gtffile, singleGens = TRUE, transposon = NULL, minoverlap = 10,  cores = 4) {
  gtfGr <- rtracklayer::import.gff(con = gtffile)
  message("Remove transcripts missing strand information.")
  gtfGr <- gtfGr[!GenomicRanges::strand(gtfGr) == "*"]
  txdb <- GenomicFeatures::makeTxDbFromGRanges(gr = gtfGr)
  exonicParts <- GenomicFeatures::exonicParts(txdb = txdb, 
                                              linked.to.single.gene.only = singleGens)
  exonrank <- split(x = exonicParts$exonic_part, 
                    f = exonicParts$gene_id, drop = TRUE)
  gestrand <- split(x = GenomicRanges::strand(exonicParts), 
                    f = exonicParts$gene_id, drop = TRUE)
  
  message("Using ", cores, " cores for exon rank ordering.")
  exonicpart <- BiocParallel::bplapply(names(exonrank), BPPARAM = BiocParallel::MulticoreParam(cores), FUN = function(x) {
    if (unique(as.character(gestrand[[x]])) == "-") {
      order(as.integer(exonrank[[x]]), decreasing = TRUE)
    } else {
      order(as.integer(exonrank[[x]]))
    }
  })

  names(exonicpart) <- names(exonrank)
  exonicpart <- exonicpart[unique(exonicParts$gene_id)]
  exonicParts$exonic_part <- unlist(exonicpart)

  if (!is.null(transposon)) {
    if (!c("names") %in% colnames(S4Vectors::mcols(transposon)) || !is(transposon, "GRanges")) {
      stop("Transposone must be a Granges object that includes a column named 'names'.")
    }
    overlaps <- IRanges::findOverlaps(query = exonicParts, subject = transposon, 
                                      minoverlap = minoverlap)
    repeats <- split(x = transposon$names[S4Vectors::subjectHits(overlaps)], 
                     f = S4Vectors::queryHits(overlaps))
    repeats <- lapply(X = repeats, FUN = function(x) paste(x, collapse = ","))
    exonicParts$transposon <- "none"
    exonicParts$transposon[as.numeric(names(repeats))] <- unlist(repeats)
  }

  return(exonicParts)
}

#' @title Consolidation of information on transcript assembly of multiple samples
#' @description Transcript-quantification is prerequisite for many downstream investigations.
#' Several metrics have been proposed for measuring abundance in transcript level based on RNA-seq data.
#' @param reference Compared with a reference annotation file
#' @param bamFiles Character vector. BAM files (sorted by coordinate).
#' @param gtfFiles Character vector. Output GTF files (will be overwritten after re-quantification).
#' @param params Other parameters for StringTie
#' @param longRead Logical scalar or vector: TRUE for ONT/PacBio
#' @param cores Number of cores to use in parallel (default: 1)
#' @export
stringtieCombine <- function(reference = NULL, bamFiles = NULL, gtfFiles = NULL,
                             params = "", longRead = FALSE, cores = 1) {
  if (is.null(reference)) stop("Please provide the reference annotation file.")
  if (is.null(bamFiles) || is.null(gtfFiles)) stop("Please provide both `bamFiles` and `gtfFiles`.")
  if (length(bamFiles) != length(gtfFiles)) stop("`bamFiles` and `gtfFiles` must have the same length.")

  ## step 1: Re-quantify with -e (overwrite requested) - parallel version
  rq_params <- paste("-e", params)
  message("Re-quantifying transcript abundance using ", cores, " core(s)...")
  parallel::mclapply(seq_along(bamFiles), function(i) {
    Teddy::stringtieAssembly(
      bam       = bamFiles[i],
      reference = reference,
      outfile   = gtfFiles[i],
      params    = rq_params,
      longRead  = longRead
    )
  }, mc.cores = cores)

  ## step 2: Preprocess gtf
  gtfGR <- rtracklayer::import.gff(reference)
  index <- which(gtfGR$type == "transcript")
  gtfGR$gene_name <- rep(gtfGR$gene_name[index], c(index[2:length(index)], length(gtfGR) + 1) - index)
  gtfGR$gene_name <- ifelse(is.na(gtfGR$gene_name), gtfGR$gene_id, gtfGR$gene_name)
  transcriptGR <- gtfGR[index]

  ## step 3: Extract quantification results
  SElist <- lapply(X = gtfFiles,
                   FUN = .ExtractTranscript,
                   transcriptGR = transcriptGR)

  ## step 4: Create SummarizedExperiment object
  SE <- do.call(IRanges::cbind, SElist)
  S4Vectors::metadata(SE) <- list(gtf = gtfGR)
  return(SE)
}

#' @title Counting reads on the exon
#' @description Counting the number of reads that fall into each exon
#' bin defined in the flattened GTF/GFF object.
#' @param annotation A GRanges object, typically from prepareAnno().
#' @param bamfile Character vector of BAM file(s).
#' @param isPairedEnd Logical scalar or vector, whether the library is paired-end.
#' @param strandSpecific Integer: 0 (unstranded), 1 (forward), or 2 (reverse).
#' @param allowMultiOverlap Logical, whether to allow reads to overlap multiple features.
#' @param isLongRead Logical, whether this is long-read data (e.g. ONT, PacBio).
#' @param nthreads Number of threads to use.
#' @param ... Additional arguments passed to featureCounts().
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Rsubread featureCounts
#' @export
countAnno <- function(annotation, bamfile,
                      isPairedEnd = TRUE,
                      strandSpecific = 0,
                      allowMultiOverlap = TRUE,
                      isLongRead = FALSE,
                      nthreads = 1,
                      ...) {
  
  annframe <- as.data.frame(annotation)
  names(annframe)[names(annframe) == "seqnames"] <- "Chr"
  annframe$GeneId <- paste(annotation$gene_id, annotation$exonic_part, sep = ":")
  
  fc_args <- list(
    files = bamfile,
    annot.ext = annframe,
    isGTFAnnotationFile = FALSE,
    useMetaFeatures = FALSE,
    isPairedEnd = isPairedEnd,
    strandSpecific = strandSpecific,
    allowMultiOverlap = allowMultiOverlap,
    nthreads = nthreads,
    ...
  )
  
  if (isLongRead) {
    fc_args$isLongRead <- TRUE
    if ("maxMOp" %in% names(fc_args)) {
      fc_args$maxMOp <- NULL
      warning("Long-read mode detected: Ignoring 'maxMOp'. Only 'isLongRead=TRUE' will be used.")
    }
  }
  
  fc_result <- do.call(Rsubread::featureCounts, fc_args)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    list(counts = as.matrix(fc_result$counts)),
    rowRanges = annotation
  )
  return(se)
}



#' Extract GTF information
#'
#' This function filters and extracts GTF information from a `SummarizedExperiment` object 
#' based on the specified element type (exon, transcript, or both) and an FPKM threshold. 
#'
#' @param combineSE A \code{SummarizedExperiment} object, typically the output from the \code{\link{combineSE}} function.
#' It should contain GTF metadata in its \code{@metadata$gtf} slot and FPKM values in one of its assays.
#' @param filter A numeric value specifying the minimum FPKM threshold for transcripts to be included in the output. 
#' Default is 1.
#' @param type A character string specifying the type of genetic elements to extract. 
#' Valid options are "exon", "transcript", or "both". Default is "exon".
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @return A filtered subset of the GTF metadata from the `combineSE` object, 
#' including only the specified types of genetic elements that meet the FPKM threshold.
#' @export
extractGTF <- function(combineSE, filter = 1, type = c("exon", "transcript", "both")){
  if (!inherits(combineSE, "SummarizedExperiment")) {
    stop("combineSE must be a SummarizedExperiment object")
  }
  type <- match.arg(type)
  Transcript_FPKM_index <- rowSums(assays(combineSE)[["FPKM"]] > filter) >= 1
  sub_combineSE <- subset(combineSE, subset = Transcript_FPKM_index)
  GTF <- sub_combineSE@metadata$gtf
  if (type == "exon") {
    type_index <- GTF$type == "exon"
    GTF <- GTF[type_index, ]
  } else if (type == "transcript") {
    type_index <- GTF$type == "transcript"
    GTF <- GTF[type_index, ]
  }
}



#' @title Testing for differential TE-chimeric exon usage
#' @description Detect to what extent TE-chimeric exon affect the expression of the corresponding transcript
#' @param SEobject An object of RangedSummarizedExperiment class, out from \bold{countAnno}.
#' @param condition Vector, indicating the experimental condition of the samples.
#' @param annotation The genome annotation object, out from \bold{prepareAnno}.
#' @param maxit control parameter: maximum number of iterations to allow for convergence when calculating dispersion.
#' @param niter whether to print messages at each step
#' @param quiet number of times to iterate between estimation of means and estimation of dispersion
#' @param warning whether to print warning at each step
#' @import SummarizedExperiment 
#' @import DESeq2
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom MatrixGenerics rowVars
#'
#' @export

ChimericDrivenTest <- function(SEobject,
                               condition = NULL,
                               annotation = NULL,
                               maxit=100,
                               niter=10,
                               quiet=FALSE,
                               warning=FALSE){
  if (!warning) {
    options(warn=-1)
  }
  if (is.null(condition)) {
    stop("Pleses set the condition for the SummarizedExperiment object.")
  }
  if (is.null(SEobject@assays@data[['counts']])) {
    stop("The count in the SummarizedExperiment object is needed.")
  }
  if (is.null(annotation)) {
    annotation <- rowRanges(SEobject)
  }
  design <- formula(x = "~ sample + chimeric + condition:chimeric")
  reducedModel <- formula(x = "~ sample + chimeric")
  
  colData(SEobject)$condition <- condition
  mcols(SEobject)$chimeric <- annotation$transposon
  meta <- as.data.frame(mcols(SEobject))
  sampleInfo <- colData(SEobject)
  featureID <- sprintf("E%3.3d", meta$exonic_part)
  groupID <- as.character(meta$gene_id)
  transcripts <- as.list(meta$tx_name)
  count <- as.matrix(SEobject@assays@data[['counts']])
  Nrow <- nrow(count)
  if(length(groupID) != Nrow)
    stop( "The length of 'groupID' must be the same as the number of rows in count matrix!", 
          call.=FALSE )
  if(length(featureID) != Nrow)
    stop( "The length of 'featureID' parameter must be the same as the number of rows of countData!", 
          call.=FALSE )
  modelInfo <- cbind.data.frame(sample = rownames(sampleInfo), sampleInfo)
  modelInfo <- rbind.data.frame(cbind(modelInfo, chimeric = "chimeric"), 
                                cbind(modelInfo, chimeric= "others"))
  rownames(modelInfo) <- NULL
  vars <- all.vars(design)
  if (any(!vars %in% colnames(modelInfo))) {
    stop("The variable 'sample' in the design formula is not in the 'modelInfo'!")
  }
  gene_exons <- split(seq(nrow(count)), as.character(groupID))
  identify_chimeric <- unlist(lapply(gene_exons, function(i) {
    subchimeric <- meta$chimeric[i]
    ifelse(any(subchimeric != "none"), TRUE, FALSE)
  }))
  chiobject <- SEobject[unlist(gene_exons[identify_chimeric]), ]
  chi_meta <- as.data.frame(mcols(chiobject))
  chi_featureID <- sprintf("E%3.3d", chi_meta$exonic_part)
  chi_groupID <- as.character(chi_meta$gene_id)
  chi_TEclass <- as.character(chi_meta$chimeric)
  chi_transcripts <- as.list(chi_meta$tx_name)
  chi_count <- as.matrix(chiobject@assays@data[['counts']])
  rownames(chi_count) <- paste(chi_groupID, chi_featureID, sep=":")
  chi_Nrow <- nrow(chi_count)
  chi_gene_exons <- split(seq(nrow(chi_count)), as.character(chi_groupID))
  
  #chi_exon with others
  others <- lapply(chi_gene_exons, function(i) {
    transposon <- chi_meta$chimeric[i]
    chi_idx <- which(transposon != "none")
    subSE <- chi_count[i, , drop = FALSE]
    other_idx <- which(transposon == "none")
    sumothers <- t(vapply(seq(length(chi_idx)), function(r) colSums(subSE[other_idx, , drop = FALSE]), numeric(dim(subSE)[2])))
    rownames(sumothers) <- rownames(subSE)[chi_idx]
    sumothers
    })
  others <- do.call(rbind, others)
  chi_exon_ids <- chi_meta$chimeric != "none"
  chimeric <- chi_count[chi_exon_ids, ]
  finalcount <- cbind(chimeric, others)
  chi_se <- SummarizedExperiment(finalcount, colData = modelInfo)
  mcols(chi_se)$featureID <- chi_featureID[chi_exon_ids]
  mcols(chi_se)$groupID <- chi_groupID[chi_exon_ids]
  mcols(chi_se)$TEclass <-chi_TEclass[chi_exon_ids]
  mcols(chi_se)$exonBaseMean <- rowMeans(finalcount)
  mcols(chi_se)$exonBaseVar <- rowVars(finalcount)
  normFactors <- estimateSizeFactorsForMatrix(count, median)
  chi_sizeFactors <- rep(normFactors, 2)
  modelInfo$sizefactor <- chi_sizeFactors
  normalizeSEcount <- function(object,sizefactors){
    t(t(assay(object)) / sizefactors)
  }
  allZero <- unname(rowSums(finalcount) == 0 |
                      rowSums(normalizeSEcount(chi_se,chi_sizeFactors)[, modelInfo$chimeric=="others"]) == 0)
  mcols(chi_se)$allZero <- allZero
  chi_dds <- DESeqDataSet(chi_se, design, ignoreRank = TRUE)
  colData(chi_dds)$sizeFactor <- chi_sizeFactors
  fullModelM <- rmDepCols(model.matrix(design, modelInfo))
  reducedModelM <- rmDepCols(model.matrix(reducedModel, modelInfo))
  chi_dds <- DESeq2::estimateDispersionsGeneEst(chi_dds, maxit = maxit, quiet = quiet,
                                        modelMatrix = fullModelM,
                                        niter = niter)
  chi_dds <- DESeq2::estimateDispersionsFit(chi_dds, fitType = "parametric", minDisp = 1e-08, quiet = quiet)
  dispersion_priorvar <- DESeq2::estimateDispersionsPriorVar(chi_dds, minDisp = 1e-08, modelMatrix = fullModelM)
  chi_dds <- DESeq2::estimateDispersionsMAP(chi_dds, dispPriorVar = dispersion_priorvar, modelMatrix = fullModelM)
  chi_test <- DESeq2::nbinomLRT(chi_dds, reduced = reducedModelM, full = fullModelM)
  return(chi_test)
}


#' @title Extract the result from the differentially expressed TE-chimeric exon test 
#' @param object An object of DE TE-chimeric exon test, out from \bold{ChimericDrivenTest}.
#' @export
extractTest <- function(object){
  #if(Filter){
  #  findex <- !is.na(results(object)$padj)
  #}else{
  #  findex <- rep(TRUE, nrow(object))
  #}
  #findex <- findex & !mcols(object)$allZero
  LRTout <- DESeq2::results(object)
  chimericCounts_N <- counts(object, normalized = TRUE)[, colData(object)$chimeric == "chimeric"]
  LRTout$exonExpr <- rowMeans(chimericCounts_N)
  LRTout$featureID <- mcols(object)$featureID
  LRTout$groupID <- mcols(object)$groupID
  LRTout$TEclass <- mcols(object)$TEclass
  LRTout$dispersion <- mcols(object)$dispersion
  LRTout$padj <- p.adjust(LRTout$pvalue, method = "BH")
  LRTout <- LRTout[, c("groupID", "featureID", "exonExpr", "TEclass", "dispersion", "stat", "pvalue", "padj")]
  # LRTout <- LRTout[findex,]
  mcols(LRTout)[1:5, 1] <- "input"
  mcols(LRTout)[1:5, 2] <- c("TranscriptID", "ExonID", "Expression of the exon", "TEclass", "Dispersion estimate among the transcript")
  return(LRTout)
}

#' @title Calculating the fold change
#' @description TE-chimeric exon usage fold changes are calculated based on the coefficients of the GLM fit
#' @param object An object of test
#' @param genes Specify the genes for the calculation of foldchanges
#' @param corssVar Default as "condition"
#' @importFrom BiocParallel bplapply
#' @import statmod methods 
#' @importFrom stats as.formula
#' @export
calculateFoldchange <- function(object,
                                genes,
                                crossVar = "condition"){
  frm <- stats::as.formula(paste("count ~", crossVar, "*Chimeric"))
  findex <- !is.na(DESeq2::results(object)$padj) & !mcols(object)$allZero
  groups <- mcols(object)$groupID
  features <- mcols(object)$featureID
  modelInfo <- colData(object)
  N_samples <- nrow(modelInfo) / 2
  disp <- DESeq2::dispersions(object)
  disp[is.na(object)] <- 1e-6
  chimericCounts <- counts(object)[,colData(object)$chimeric == "chimeric"]
  otherCounts <- counts(object)[,colData(object)$chimeric == "others"]
  samples <- unique(colData(object)$sample)
  mfSmall <- as.data.frame(colData(object))
  if (!is.null(genes)) {
    alleffects <- BiocParallel::bplapply( genes,
                            getEffectsForGene,
                            groups = groups, 
                            otherCounts = otherCounts,
                            chimericCounts = chimericCounts,
                            disp = disp,
                            mf = mfSmall, 
                            frm = frm,
                            findex = findex,
                            N_samples = N_samples,
                            samples = samples,
                            features = features,
                            object = object)
    names(alleffects) <- genes
  }
  return(alleffects)
}








#' Process GTF and SummarizedExperiment for TE overlap analysis
#'
#' This function processes a GTF file or a SummarizedExperiment object to identify overlaps with transposable elements (TE).
#' It can handle both external (local) GTF files and in-memory SummarizedExperiment objects. 
#' The function extracts exon and transcript information, finds overlaps with TE annotations, and filters out transcripts with zero expression.
#'
#' @param gtfFile Character. Path to the GTF file (optional, required if combineSE is NULL).
#' @param te GRanges object. Transposable element annotations, required for overlap analysis.
#' The `te` object must contain at least three metadata columns: `names`, `family`, and `class`. 
#' @param combineSE SummarizedExperiment. Pre-processed RNA-seq quantification object (optional, required if gtfFile is NULL).
#' @param minoverlap Integer. Minimum required overlap (in base pairs) between GTF exons and TE annotations.
#' @param threads Integer. Number of threads to use for parallel processing. Default is 4.
#' 
#' @return A GRanges object with transcripts overlapping TEs, annotated with TE information.
#' 
#' @import BiocParallel
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' 
#' @export
processGTF <- function(te, gtf_path = NULL, combineSE = NULL, minoverlap = 0,threads = 9) {
  
  # Load GTF data from either combineSE or GTF path
  if (!is.null(combineSE)) {
    GTF <- combineSE@metadata$gtf
  } else if (!is.null(gtf_path)) {
    #GTF <- readRDS(gtf_path)
    GTF <- rtracklayer::import(gtf_path)
  } else {
    stop("Please provide either a GTF path or a combineSE object.")
  }
  
  # Filter to include only exons
  exon_index <- GTF$type == "exon"
  GTF <- GTF[exon_index, ]
  
  # Filter out '*' strand
  GTF <- GTF[GenomicRanges::strand(GTF) != "*"]
  
  # Rank exons based on strand direction
  t_exonRank <- split(x = GTF$exon_number, f = GTF$transcript_id)
  t_exonRank <- t_exonRank[unique(GTF$transcript_id)]
  
  t_gestrand <- split(x = GenomicRanges::strand(GTF), f = GTF$transcript_id, drop = TRUE)
  t_gestrand <- t_gestrand[unique(GTF$transcript_id)]
  
  bp_param <- MulticoreParam(workers = threads)  
  
  exoninTx <- bplapply(names(t_exonRank), function(x) {
    if (unique(as.character(t_gestrand[[x]])) == "-") {
      order(as.integer(t_exonRank[[x]]), decreasing = TRUE)
    } else {
      order(as.integer(t_exonRank[[x]]))
    }
  }, BPPARAM = bp_param)
  
  
  
  names(exoninTx) <- unique(GTF$transcript_id)
  GTF$tx_exon_rank <- unlist(exoninTx, recursive = FALSE, use.names = TRUE)
  
  
  # Find overlaps between GTF and TE
  overlaps <- findOverlaps(GTF, te, minoverlap = minoverlap)
  
  chi_GTF <- GTF[queryHits(overlaps)]
  chi_te <- te[subjectHits(overlaps)]
  
  # Add TE details to GTF
  chi_GTF$name <- chi_te$name
  chi_GTF$family <- chi_te$family
  chi_GTF$class <- chi_te$class
  
  class <- split(x = te$class[S4Vectors::subjectHits(overlaps)], f = S4Vectors::queryHits(overlaps))
  class <- lapply(class, function(x) paste(x, collapse = ","))
  
  repeats <- split(x = te$names[S4Vectors::subjectHits(overlaps)], f = S4Vectors::queryHits(overlaps))
  repeats <- lapply(repeats, function(x) paste(x, collapse = ","))
  
  GTF$TE_name <- "none"
  GTF$TE_name[as.numeric(names(repeats))] <- unlist(repeats)
  
  GTF$TE_class <- "none"
  GTF$TE_class[as.numeric(names(class))] <- unlist(class)
  
  # Filter to keep only rows with TE information
  CHI_GTF <- GTF[GTF$TE_name != "none", ]
  
  return(CHI_GTF)
}
