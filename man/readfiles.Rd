% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loaddata.R
\name{readfiles}
\alias{readfiles}
\title{Read data}
\usage{
readfiles(paths = NULL, format = "10x", ...)
}
\arguments{
\item{paths}{a path or a list of paths to the csv file/s, or the directory location/s containing the 10X files: \cr
\enumerate{
\item matrix.mtx, sparse matrix format file
\item genes.tsv, gene names tab-separated file
\item barcodes.tsv, barcodes tab-separated file
}; with format="txt", specify the location of the csv matrix file.}

\item{format}{One of c("10X", "txt"). Defaults to "10X".}

\item{...}{pass specific arguments to the \code{read.csv()} function when \code{format = "txt"}.}
}
\value{
SingleCellExperiment object.
}
\description{
Read 10X or csv data,
}
\details{
Loads expression data from 10X format directory or csv file as specified. The rows in the matrix are considered genes; the columns represent samples
}
