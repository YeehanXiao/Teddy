#' @title R wrapper to Run stringtie
#' @description R wrapper to run stringtie, enabling the reconstruction of
#' transcriptome from RNA-seq reads
#' @param bam Character, a bam file from RNA-Seq reads sorted by genomic location
#' @param reference Character, a genome reference to guide transcript assembly
#' @param gtfFiles Output GTF files
#' @param params Other parameters
#'
#' @export
stringtieAssembly <- function(bam, reference, gtfFiles, params = "") {
  reference <- paste("-G", reference, sep = " ")
  outfile <- paste0("-o", gtfFiles, sep = " ")
  cmd <- sprintf("%s %s %s %s",
                 bam,
                 reference,
                 outfile,
                 params)
  return(invisible(lapply(cmd, .Stringtiebin)))
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
#' @importFrom rtracklayer import.gff
#' @importFrom GenomicFeatures exonicParts
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#' @importFrom GenomicRanges strand
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors mcols
#' @export
prepareAnno <- function(gtffile, singleGens = TRUE, transposon = NULL, minoverlap = 10) {
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
  exonicpart <- base::lapply(names(exonrank), FUN = function(x) {
    if(unique(as.character(gestrand[[x]])) == "-") {
      exonrank[[x]] <- order(as.integer(exonrank[[x]]),decreasing = TRUE)
    } else {
      exonrank[[x]] <- exonrank[[x]]
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

#' @title Consolidataion of information on transcript assembly of multiple samples
#' @description Transcript-quantification is prerequisite for many downstream investigations.
#'  Several metrics have been proposed for measuring abundance in transcript level based on
#'  RNA-seq data.
#' @param reference Compared with a reference annotation files
#' @param bamfile Bamfile, must be of SAM/BAM/CRAM format, sorted by their genomic location.
#' @param gtfFiles Output GTF files
#' @param params Other parameters
#' @export
#' 
stringtieCombine <- function(reference = NULL, bamfile = NULL, gtfFiles = NULL, params = "") {
  # step 1: Quantify
  params = paste("-e", params)
  stringtieAssembly(bam = bamfile, 
                    reference = reference, 
                    gtfFiles = gtfFiles, 
                    params = params)
  
  # step 2: Preprocess gtf
  gtfGR <- rtracklayer::import.gff(reference)
  index <- which(gtfGR$type == "transcript")
  gtfGR$gene_name <- rep(gtfGR$gene_name[index], c(index[2:length(index)], length(gtfGR) + 1) - index)
  gtfGR$gene_name <- ifelse(is.na(gtfGR$gene_name), gtfGR$gene_id, gtfGR$gene_name)
  transcriptGR <- gtfGR[index]
  
  # step 3: Extract quantification results
  SElist <- lapply(X = gtfFiles, 
                   FUN = .ExtractTranscript, 
                   transcriptGR = transcriptGR)
  
  
  # step 4: Create SummarizedExperiment object
  SE <- do.call(IRanges::cbind, SElist)
  S4Vectors::metadata(SE) <- list(gtf = gtfGR)
  return(SE)
}

#' @title Counting reads on the exon
#' @description Counting the number of reads that falls into each of the exon
#' counting bins defined in the flattened GFF object.
#' @param annotation Flattened GFF object.
#' @param bamfile Name of the BAM file(s).
#' @param strandSpecific An integer vector indicating if strand-specific read
#' counting should be performed. See \code{\link[Rsubread]{featureCounts}}
#' @param allowMultiOverlap Logical indicating if a read is allowed to be
#' assigned to more than one feature (or meta-feature) if it is found to
#' overlap with more than one feature (or meta-feature). TRUE by default.
#' @param isPairedEnd A logical scalar or a logical vector, indicating whether
#' libraries contain paired-end reads or not.
#' @param nthreads An integer specifying the number of threads to be used.
#' @param ... Additional arguments. See \code{\link[Rsubread]{featureCounts}}
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Rsubread featureCounts
#'
#' @export
countAnno <- function(annotation, bamfile,
                      isPairedEnd = TRUE,
                      strandSpecific = 0,
                      allowMultiOverlap = TRUE,
                      nthreads = 1,
                      ...) {
  #bam <- BamFileList(bamfile)
  #count <- GenomicAlignments::summarizeOverlaps(
  #          features = annotation,
  #          reads = bam,
  #          singleEnd = singleEnd,
  #          ...)
  #return(count)
  annframe <- as.data.frame(annotation)
  names(annframe)[names(annframe) == "seqnames"] <- "Chr"
  annframe$GeneId <- paste(annotation$gene_id, annotation$exonic_part, sep = ":")
  hints <- Rsubread::featureCounts(files = bamfile,
                                   annot.ext = annframe,
                                   strandSpecific = strandSpecific,
                                   allowMultiOverlap = allowMultiOverlap,
                                   isPairedEnd = isPairedEnd,
                                   useMetaFeatures = FALSE,
                                   isGTFAnnotationFile = FALSE,
                                   nthreads = nthreads,
                                   ...)
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = as.matrix(hints$counts)),
                                                   rowRanges = annotation)
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
#' @import statmod methods stats
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
