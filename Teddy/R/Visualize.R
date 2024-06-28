#' @title Plot the formation of a specific transcript
#' @param GTF A Granges object, extracted from the output of stringtieCombine\code{\link{stringtieCombine}},
#' filtered to rows with the type "exon". Use the \code{loadData} function to view the specific format.
#' @param geneName The gene name to which the transcript belongs.
#' @param TEname The TE element involved in the fusion
#' @param txid The transcript ID of the TE-chimeric transcript for which the formation is to be plotted.
#' @param rank Int. The index of the exon in which the TE element is located in the chimeric event. 
#' This information can be obtained by referring to \item{annotation} or the output of \code{\link{stringtieCombine}}.
#' @importFrom ggplot2 ggplot geom_rect geom_segment geom_text
#' @export
#'   
formPlot <- function(GTF,txid = NULL,rank = 1,geneName = "Gene", TEname = "MT2B2") {
  pic_table <- GTF[GTF$transcript_id == txid]
  if(length(pic_table) == 0) {
    print("The provided Transcript ID was not found in the GTF.")
    return()
  }
  N_exon <- length(pic_table)
  fill <- rep("#ECD1D6",N_exon)
  fill[rank] <- "#3A3B4F"
  picdata <- data.frame(
    x = seq(1,1*N_exon,1),
    y = rep(0.5,N_exon),
    width = rep(0.5,N_exon),
    height = rep(0.1, N_exon),
    fill =fill
  )
  median <- (picdata$x[1:N_exon-1]+picdata$x[2:N_exon])/2
  p <- ggplot(picdata) +
    geom_rect(aes(xmin = x-(width/2), xmax = x+(width / 2),
                  ymin = y-(height/2), ymax = y+(height / 2)),
              fill = fill, color = "black")+
    xlim(0, N_exon+1) + ylim(0, 1) +   
    theme_void()
  pic_segment_1 <- data.frame(
    x = picdata$x[1:N_exon-1]+0.25,
    xend = median,
    y = rep(0.55,N_exon-1),
    yend = rep(0.6,N_exon-1)
  )
  pic_segment_2 <- data.frame(
    x = median,
    xend = picdata$x[2:N_exon]-0.25,
    y = rep(0.6,N_exon-1),
    yend = rep(0.55,N_exon-1)
  )
  p <- p + geom_segment(data = pic_segment_1, 
                        aes(x = x, 
                            y = y,
                            xend = xend, 
                            yend = yend),
                        color = "black") +
    geom_segment(data = pic_segment_2,
                 aes(x = x, 
                     y = y,
                     xend = xend, 
                     yend = yend),
                 color = "black") 
  p <- p + geom_text(x = picdata$x[rank], y = 0.4, label = TEname,
                     color = "black")
  p <- p + geom_text(x = 0.7, y = 0.7, label = paste(geneName,txid,sep = ":"),
                     color = "black", hjust = 0)
  return(p)
}


#' @title Plot the gene model in the control and reference groups
#' @param count A data matrix of bin-based counts，extracted from the output of countAnno\code{\link{countAnno}}
#' @param conditions A factor vector indicating the conditions (e.g., control and reference) of each sample
#' @param annotation A Granges object for annotated transcripts
#' @param idx The index of the gene.
#' @param labels A character vector, containing labels for the control and reference groups.
#' @param gene_name The gene name.
#' @param Title The title of plot, e.g., gene name.
#' @param gtf The dir for merged GTF file as the output of \code{\link{gffcompareAnno}}
#' @param txid The transcript id of a specific transcript for the gene
#' @param chi_test DESeq2 object. The output from the ChimericDrivenTest, extracted by \code{\link{extractTest}}
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @export
#' 
diffBinPlot <- function(count,
                        conditions,
                        annotation,
                        idx,
                        labels,
                        Title,
                        gtf,
                        txid,
                        chi_test,
                        gene_name,
                        ...) {
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = data.frame(condition = conditions),
                                design = ~ condition)
  dds_vst <- vst(dds)
  count_vst <- assay(dds_vst)
  subcounts <- count_vst[idx,]
  #
  texty <- "Expression"
  linecolor <- "#ECD1D6"
  label <- paste0("E0",annotation[idx]$exonic_part)
  colors <- c("#9BB0CB", "#8c3837")
  #
  layout(matrix(c(1,1,2), nrow = 3))
  p <-plot.new()
  par(mar = c(0,5,3,2))
  p <-plot.window(xlim=c(0, 1), ylim=c(0, max(subcounts)+1))
  intervals<-(0:nrow(subcounts))/nrow(subcounts)
  middle <- (intervals[1:length(intervals)-1]+intervals[2:length(intervals)])/2
  label <- paste0("E0",seq_len(nrow(subcounts)))
  axis(1, at=middle, labels=label)
  axis(2)
  mtext(texty,side=2,line=1.8)
  level=3
  result <- as.data.frame(subcounts) %>%
    dplyr::reframe(untreat = rowMeans(across(1:level)),
                   treat = rowMeans(across((level+1):ncol(count))))
  abline(v=middle, lty="dotted", col="#ECD1D6")
  xflatsegEnd <- intervals[1:length(intervals)-1]+0.8*(1/nrow(subcounts))
  segments(x0 = intervals[1:length(intervals)-1],
           x1 = xflatsegEnd,
           y0 = result[,1],
           y1 = result[,1],
           col  = "#9BB0CB",
           lwd = 1.8)
  segments(x0 = intervals[1:length(intervals)-1],
           x1 = xflatsegEnd,
           y0 = result[,2],
           y1 = result[,2],
           col  = "#8c3837",
           lwd = 1.8)
  segments(x0 = xflatsegEnd[1:length(xflatsegEnd)-1],
           x1 = intervals[2:(length(intervals)-1)],
           y0 = result[1:(dim(result)[1]-1),1],
           y1 = result[2:dim(result)[1],1],
           col="#9BB0CB", 
           lty="dotted",
           lwd = 1.2) 
  segments(x0 = xflatsegEnd[1:length(xflatsegEnd)-1],
           x1 = intervals[2:(length(intervals)-1)],
           y0 = result[1:(dim(result)[1]-1),2],
           y1 = result[2:dim(result)[1],2],
           col="#8c3837", 
           lty="dotted",
           lwd = 1.2) 
  #
  gtf <- rtracklayer::import.gff(gtf)
  form_GR <- gtf[gtf$transcript_id==txid,]
  form_GR <- form_GR[form_GR$type=="exon"]
  overlaps <- findOverlaps(form_GR,annotation[idx,])
  map_GR <- annotation[idx,][subjectHits(overlaps)]
  exonicparts <- as.data.frame(annotation)[idx,"exonic_part"]
  exonicCHI <- as.data.frame(annotation)[idx,"transposon"]
  tx_loc <- exonicparts[subjectHits(overlaps)]
  tx_CHI <- exonicCHI[subjectHits(overlaps)]
  names(tx_loc)<- tx_CHI
  subset <- chi_test[complete.cases(chi_test$geneName), ]
  sig_loc <- as.numeric(gsub(".*E0*(\\d+)$", "\\1", subset[subset$geneName == gene_name, "featureID"]))
  sig_x <- middle[sig_loc]
  sig_y <- result[sig_loc,2]+0.3
  TE_y <- result[sig_loc,2]-0.2
  points(sig_x, sig_y, pch = "*", cex = 1.5)
  text(sig_x, TE_y, sapply(strsplit(subset[subset$geneName == gene_name, "TEclass"], ","), function(x) (x[1])), cex = 0.5)
  #
  p <-plot.new()
  par(mar = par(mar = c(0,0,5,0)))
  starts <- start(map_GR)
  ends <- end(map_GR)
  widths <- ends - starts
  plot.window(xlim=c(0,2*sum(ends-starts)), ylim=c(0, 1))
  rango <- seq_len(length(starts))
  gaplength <- starts[2:length(starts)]-ends[1:(length(ends)-1)]
  scalegap <- gaplength/sum(gaplength)*sum(ends-starts)
  cum_width <- cumsum(widths)
  cum_gap <- cumsum(scalegap)
  scalestarts <- c(0,cum_width[1:length(cum_width)-1]+cum_gap)
  scaleends <- cum_width+c(0,cum_gap)
  scalemid <- (scaleends[1:length(scaleends)-1]+scalestarts[2:length(scalestarts)])/2
  col <- ifelse(tx_CHI!="none", "#ECD1D6", "#3A3B4F")
  rect(scalestarts[rango], 0.2, scaleends[rango], 0.8,col = col)
  segments(scaleends[1:length(scaleends)-1], 0.4, scalemid, 0.5)
  segments(scalemid, 0.5, scalestarts[2:length(scalestarts)], 0.4)
  text(x=sum(ends-starts),y=0.1,label = txid)
  scaleintervals <- intervals*2*sum(ends-starts)
  segments(x0 = 0.5*(scalestarts+scaleends),
           x1 = middle[tx_loc]*(2*sum(ends-starts)),
           y0 = 0.8,
           y1 = 1,
           lty="dotted")
  
}