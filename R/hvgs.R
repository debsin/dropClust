#' Find highly variable genes
#' @description Get variable genes from normalized UMI counts using Fano Factor metric.
#' @details Compute Fano Factor metric for each gene. The metric computes the median absolute deviation of dispersion across multiple bins for each gene.
#' @param object A SingleCellExperiment object containing normalized expression values in \code{"normcounts"}.
#' @param ngenes_keep integer to return top ranking \code{ngenes_keep} number of genes.
#' @return  A SingleCellExperiment object with an additional column named \code{HVG} in \code{rowData} column.
#' The column stores a a logical value against each gene to indicate if it has been ranked within the top
#' \code{ngenes_keep}. It also generates an additional column \code{dispersion_norm} in \code{rowData} to
#' store the dispersion metric against each gene.
#' @importFrom SingleCellExperiment normcounts rowData
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
#' sce <- RankGenes(sce)
RankGenes<-function(object, ngenes_keep=1000) {

  avg<-Matrix::rowMeans(normcounts(object))
  dispersion <- ColDispersion(Matrix::t(normcounts(object)))

  avg_bin<-cut(avg,breaks=c(-Inf,stats::quantile(avg,seq(0.1,1,0.05)),Inf))
  df=as.data.frame(cbind(avg,dispersion,avg_bin))

  variance_by_bin<-plyr::ddply(df,"avg_bin",function(x) {
    data.frame(bin_median=stats::median(x$dispersion),
               bin_mad=stats::mad(x$dispersion))
  })
  df$bin_disp_median<-
    variance_by_bin$bin_median[match(df$avg_bin, variance_by_bin$avg_bin)]
  df$bin_disp_mad<-
    variance_by_bin$bin_mad[match(df$avg_bin, variance_by_bin$avg_bin)]
  df$dispersion_norm<-
    with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)

  SummarizedExperiment::rowData(object)$dispersion_norm <- df$dispersion_norm

  cat("Sort Top Genes...\n")
  disp_cut_off<-sort(df$dispersion_norm,decreasing=TRUE)[ngenes_keep]
  cat("Cutoff Genes...\n")
  SummarizedExperiment::rowData(object)$HVG<-df$dispersion_norm >= disp_cut_off

  object@metadata[["dropClust"]] = c(unlist(object@metadata[["dropClust"]]), "RankGenes")


  object
}

# ---------------------------------------------
# Select PCA based Genes
# ---------------------------------------------
#' Select PCA based Genes
#' @description Performs gene selection on sampled cells based on PCA loadings
#' @details Genes are ranked for selection in 3 steps:\cr
#' \enumerate{
#' \item First 50 principal components are obtained using Singular value Decomposition is used as implemented in the \code{irlba} R package.
#' \item Among the first 50 components, top 10 components are selected in the order of their modality.
#' \item Genes are ordered based on their average loadings in the rotation matrix containing the top 10 components.
#' }
#' @param object A SingleCellExperiment object containing normalized expression values in \code{"normcounts"}.
#' @param top integer specifying to number of genes to return in their order of ranking
#' @return A SingleCellExperiment object with an additional column named \code{PCAGenes} in \code{rowData} column.
#' The column stores a a logical value against each gene to indicate if it has been ranked within the \code{top}.
#' @export
#' @examples
#' library(SingleCellExperiment)
#' ncells <- 1000
#' ngenes <- 2000
#' x <- matrix(rpois(ncells*ngenes, lambda = 10), ncol=ncells, nrow=ngenes, byrow=TRUE)
#' rownames(x) <- paste0("Gene", seq_len(ngenes))
#' colnames(x) <- paste0("Cell", seq_len(ncells))
#' sce <- SingleCellExperiment(list(counts=x))
#' sce <- CountNormalize(sce)
#' sce <- RankGenes(sce)
#' sce <- Sampling(sce)
#' sce <- RankPCAGenes(sce)
RankPCAGenes<-function(object,top=200){
  if(is.null(rowData(object)$HVG) || is.null(colData(object)$Sampling)){
    stop("SingleCellExperiment object does not contain high variable genes or Sampling() is not called prior.")
  }

  subsamples = SingleCellExperiment::colData(object)$Sampling
  hvg_genes = SingleCellExperiment::rowData(object)$HVG

  mat = Log2Normalize(normcounts(object)[hvg_genes,subsamples], return.sparse = F)

  mat = t(mat)
  cat("Find best PCA components...")

  # mat = assay(object,i = "logcounts")[rowData(object)$HVG, colData(object)$Sampling]
  nc = 50
  invisible(print(dim(mat)))
  PRR <- irlba::prcomp_irlba(mat, nc)
  S = c()
  for(i in 1:nc)
  {
    K = PRR$x[,i]

    # if(prev){
    #   # L = mclust::Mclust(K)
    #   # S = c(S,L$G)
    #   peaks= NbClust::NbClust(K, distance = "euclidean",
    #                  min.nc = 2, max.nc = 10,
    #                  method = "kmeans", index ="ball")$Best.nc[1]
    # }
    # else{
    #   peaks = all_peaks(K)
    # }

    peaks = all_peaks(K)
    S = c(S,peaks)
    # cat(i,"\t")
  }
  graphics::barplot(S)

  df<-as.data.frame(cbind("MOD" = S,
                          "VAR" = apply(PRR$x[,1:nc],2, stats::var)))
  rownames(df)<-paste0("PC",seq(nc))
  #gene_pca <-abs(PRR$v[,which(S>=3)])
  best_components = utils::head(with(df, order(-MOD, -VAR)), 10)
  gene_pca <-abs(PRR$rotation[,best_components])
  gene_max<-apply(gene_pca, 1, max)
  rank_gene <- order(gene_max, decreasing=TRUE)


  top_pc_genes = utils::head(rank_gene,top)

  # gene.ids = which(rowData(object)$HVG==TRUE)
  SummarizedExperiment::rowData(object)$PCAGenes <- rep(FALSE, nrow(object))
  SummarizedExperiment::rowData(object)$PCAGenes[which(hvg_genes==T)[top_pc_genes]] <- TRUE
  cat(top,"genes selected.\n")

  object@metadata[["dropClust"]] = c(unlist(object@metadata[["dropClust"]]), "RankPCAGenes")

  return(object)
}

all_peaks <- function(x) {
  den.fit <- stats::density(x, kernel=c("gaussian"))
  den.spline <- stats::smooth.spline(den.fit$x, den.fit$y,
                                     all.knots=TRUE, spar = 0.5)
  d1 <- stats::predict(den.spline, den.spline$x, deriv=1)
  peak.count <- length(rle(den.sign <- sign(d1$y))$values)/2
  if (is.na(peak.count) == TRUE) { peak.count <- 0 }
  return(peak.count)
}
