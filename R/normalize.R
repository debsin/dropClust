#' Normalize a SingleCellExperiment object using total counts
#'
#' Compute normalized expression values from count data in a SingleCellExperiment object, using the median normalized total count stored in the object.
#'
#' @param object A SingleCellExperiment object.
#' @param return_log Logical scalar, should normalized values be returned on the log2 scale?
#' If \code{TRUE}, output is stored as \code{"logcounts"} in the returned object; if \code{FALSE} output is stored as \code{"normcounts"}.#'
#' @details
#' Normalized expression values are computed by dividing the counts for each cell by median normalized total count of that cell.
#' If \code{log=TRUE}, log-normalized values are calculated by adding a pseudocount of \code{1} to the normalized count and performing a log2 transformation.
#'
#' @return A SingleCellExperiment object containing normalized expression values in \code{"normcounts"} if \code{log=FALSE},
#' and log-normalized expression values in \code{"logcounts"} if \code{log=TRUE}.
#' @importFrom SingleCellExperiment counts normcounts<- normcounts logcounts<-
#' @export
#' @examples
#' library(SingleCellExperiment)
#' ncells <- 100
#' ngenes <- 2000
#' x <- matrix(rpois(ncells*ngenes, lambda = 10), ncol=ncells, nrow=ngenes, byrow=TRUE)
#' rownames(x) <- paste0("Gene", seq_len(ngenes))
#' colnames(x) <- paste0("Cell", seq_len(ncells))
#' sce <- SingleCellExperiment(list(counts=x))
#' sce <- CountNormalize(sce)
#'
CountNormalize<- function(object, return_log = FALSE) {

  ## Compute normalized expression values.
  # x_filt = assay(object, i = exprs_values, withDimnames=FALSE)
  x_filt = Matrix::t(counts(object))
  # cs <- Matrix::colSums(x_filt)
  # cs_med <- cs/stats::median(cs)
  # assay(object, "normcounts") <- x_filt/cs_med
  rs<-Matrix::rowSums(x_filt)
  rs_med<-stats::median(rs)
  normcounts(object)<-Matrix::t(x_filt/(rs/rs_med))
  if(return_log)
    logcounts(object)<-Log2Normalize(normcounts(object))

  object@metadata[["dropClust"]] = c(unlist(object@metadata[["dropClust"]]), "CountNormalize")

  ## return object
  return(object)
}

#' Compute log normalized expression values from count data in a SingleCellExperiment object.
#'
#' @param x either a sparse or a dense matrix.
#' @param return.sparse Logical when \code{TRUE} return type is same as input, when
#' \code{FLASE} a dense matrix is returned.
#' @details
#' The input matrix is transformed into the log base 2 scale after addition of a pseudo count \code{1}.
#' @return log normalized matrix same as the input dimensions.
#'
#'
#' @export
#' @examples
#' a = matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#'
#' lnorm <- Log2Normalize(a)
Log2Normalize<- function(x, return.sparse=FALSE) {

  ## Compute normalized expression values.
  if(return.sparse)
    return(base::log2(x+1))
  else
    return(as.matrix(base::log2(x+1)))

}

