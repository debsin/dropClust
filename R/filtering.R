#' Filter out poor quality cells
#'
#' Keep only those cells (columns) expressing atleast count = \code{min_count} in the number of genes specified within the quantile range between \code{ql_th} and \code{qh_th}.
#'
#' @param object the SingleCellExperiment object to filter
#' @param min_count integer threshold for expression count
#' @param ql_th quantile at probability with values in [0,1] for lower limit
#' @param qh_th quantile at probability with values in [0,1] for upper limit
#' @return SingleCellExperiment object with the bad cells removed
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#' object <- SingleCellExperiment(assays = list(counts = counts))
#' FilterCells(object)
FilterCells <- function(object, min_count=3, ql_th = 0.001, qh_th = 1) {

  try(if(qh_th < ql_th) stop("Lower limit must be smaller than upper limit."))


  if (is(object, "SingleCellExperiment")) {
    cS <- Matrix::colSums(SingleCellExperiment::counts(object)>=min_count)

  } else {
    cS <- Matrix::colSums(object>=min_count)
  }
  l1 = stats::quantile(cS,probs = ql_th)
  l2 = stats::quantile(cS,probs = qh_th)
  keep_cells = intersect(which(cS>=l1), which(cS<=l2))

  cat(length(cS)-length(keep_cells), "bad cells removed.\n")

  return(object[, keep_cells])
}




#' Filter out poor quality genes
#'
#' Genes (rows) with count greater than \code{min_count} and not expressed in atleast \code{min_cell} cells are removed.
#'
#' @param object the SingleCellExperiment  to filter genes
#' @param min_count integer threshold for expression count
#' @param min_cell integer threshold for number of cells expressing a particular gene
#' @return SingleCellExperiment object with the bad genes removed
#' @export
#' @examples
#' library(SingleCellExperiment)
#' counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#' object <- SingleCellExperiment(assays = list(counts = counts))
#' FilterGenes(object)
FilterGenes <- function(object, min_count=2, min_cell=3) {

  mat = SingleCellExperiment::counts(object)

  rS <- Matrix::rowSums(mat > min_count)

  keep_genes = which(rS > min_cell)

  cat(length(rS)-length(keep_genes), "genes filtered out,", length(keep_genes),"genes remaining.\n")

  # object@metadata[["dropClust"]] = c(unlist(object@metadata[["dropClust"]]), "FilterGenes")

  return(object[keep_genes, ])
}
