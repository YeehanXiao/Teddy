% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{stringtieCombine}
\alias{stringtieCombine}
\title{Consolidataion of information on transcript assembly of multiple samples}
\usage{
stringtieCombine(
  reference = NULL,
  bamfile = NULL,
  gtfFiles = NULL,
  params = ""
)
}
\arguments{
\item{reference}{Compared with a reference annotation files}

\item{bamfile}{Bamfile, must be of SAM/BAM/CRAM format, sorted by their genomic location.}

\item{gtfFiles}{Output GTF files}

\item{params}{Other parameters}
}
\description{
Transcript-quantification is prerequisite for many downstream investigations.
 Several metrics have been proposed for measuring abundance in transcript level based on
 RNA-seq data.
}
