# ---------------------------------------------
# Clustering of samples using
# Hierarchical Clustering
# ---------------------------------------------
#' Clustering of samples
#' @description Performs clustering on sampled cells and Post-hoc Cluster Assignment.
#' @details Clustering is carried out in two alternate approaches on the sampled cells.
#' For the default setting or quick identification of the existing broad clusters,
#' a louvain based partition is employed. Otherwise for fine-tuned clustering with outliers,
#' hierarchical clustering is used with \code{cutreeDynamic} for dendrogram cut. Also, Assigns cluster membership to unsampled cells by using cluster membership information of the nearest neighbours.
#' An approximate nearest neighbour graph is constructed out of the samples population using the \code{find_ann()} module.
#' Some cells are left un assigned when its neighbour's cluster membership doesn't form a majority as specified by the \code{conf} parameter.
#' Unassigned cells (\code{NA}) are excluded in the plot or further downstream analysis.
#' @param data normalized expression matrix as returned by \code{matrix.subset()}.
#' The matrix columns represent the subset of genes obtained from \code{pc_genes} module.
#' Subsetting the return matrix by \code{matrix.subset()} module with samples from primary clusters and PCA selected genes.
#' @param sp.samples transcriptome indices obtained from the Structure Preserving Sampling step.
#' @param default logical, when \code{TRUE}, louvain based partition is used and other clustering parameters are ignored.
#' When \code{FALSE}, hierarchical clustering is used.
#' @param minClusterSize integer, specifies the size of the smallest cluster; ignored when \code{default = TRUE}.
#' @param deepSplit integer, level of dendrogram split [0-4], higher value produces finer clusters; ignored when \code{default = TRUE}.
#' @param conf numeric [0-1], defines the expected confidence of majority for a consensus. Cells remain unassigned when majority is below \code{conf}.
#' @param k_nn integers, specifies number of nearest neighbours, defaults to 10.
#' @return List of:\cr
#' \enumerate{
#' \item \code{cluster.ident} vector cluster identifiers ranging from 1 to the number of clusters for respective data points.\cr
#' \item \code{nn.ids}  matrix, each row corresponds to a cell, whose columns depict cluster membership of its neighbours; as returned by the \code{find_ann()} module. \cr
#' }
#' Unassigned samples are represented by\code{NA} values.
#' @export
cluster.cells<-function(data, sp.samples, default = TRUE, minClusterSize = 20, deepSplit = 3, conf=0.75, k_nn=10){

  subsamples_louvain = sp.samples

  ## ------------------------------------------------------------------------
  #### Hierarchical Clustering on subsamples
  # Adjust Minimum cluster size with argument minClusterSize (default = 20)
  # Adjust tree cut with argument level deepSplit (default = 3), higher value produces more clusters.
  # ss_sel_genes_mat<-as.matrix(data[subsamples_louvain,top_pc_genes])
  ss_clusters<-ss_clustering(as.matrix(data[subsamples_louvain,]),
                             default = default,
                             minClusterSize = minClusterSize,
                             deepSplit = deepSplit)

#
#   # Collect non-outliers
#   if(sum(is.na(ss_clusters$outliers))>0)
#     clustered_samples = subsamples_louvain else {
#       clustered_samples = subsamples_louvain[-ss_clusters$outliers]}


  ## ------------------------------------------------------------------------
  # ----------------------------
  # Find K Nearest Neighbours among sub-samples
  # ----------------------------
  gc<-gc(verbose=FALSE)
  INDEX <- find_ann(as.matrix(data),sp.samples,n_nn = k_nn)


  # ----------------------------
  # UNANNOTATED Class Assignment
  # ----------------------------

  cat("Assign Cluster Ids...\n")

  clust_col <- cluster_assign(INDEX,ss_clusters,conf)
  cat(paste("Unassigned Cells",length(clust_col[is.na(clust_col)])),"\n")

  residue = factor(clust_col[sp.samples],exclude = NA)

  keep_labels = which(clust_col %in% levels(residue))

  all.labels = clust_col
  all.labels[-keep_labels]<-NA


  lev = 1:length(unique(all.labels[stats::complete.cases(all.labels)]))
  cc.f<-factor(all.labels,exclude = NA)
  levels(cc.f)<-lev

  cat(paste("Number of Predicted Clusters:",length(levels(cc.f))),"\n")

  return(list("cluster.ident" = cc.f, "nn.ids" = INDEX))

}


ss_clustering<-function(ss_sel_genes, default = TRUE , minClusterSize = 20,
                        deepSplit = 3){
  # ss_sel_genes = ss_sel_genes_mat

  if (default==TRUE){
    hc_labs_clean = ann_partition(ss_sel_genes)[,2]

    return(list("labels"=hc_labs_clean))
  }
  cat("Perfom Hierarchical Clustering...\n")
  dim(ss_sel_genes)
  d = rdist::rdist(ss_sel_genes, metric = "angular")
  hc<-stats::hclust(d,method = "average")
  # par()
  # plot(hc)

  hc_labs<-dynamicTreeCut::cutreeDynamic(dendro = hc, cutHeight = NULL,
                                         minClusterSize = minClusterSize,
                                         method = "hybrid",
                                         deepSplit = deepSplit,
                                         pamStage = TRUE,
                                         distM = as.matrix(d),
                                         verbose = 0)

  names(hc_labs) = NULL


  return(list("labels"=hc_labs))
}



cluster_assign<-function(INDEX, ss_clusters,conf){

  nn = dim(INDEX)[2]
  clust_col = list()


  for( i in seq(nrow(INDEX)))
  {
    clust_col[[i]] = NA
    tab = table(ss_clusters$labels[INDEX[i,]])
    max_id = which.max(tab)

    if(tab[max_id]/nn > conf){
      clust_col[[i]]=names(max_id)
    }

  }


  clust_col = as.numeric(unlist(clust_col, use.names = FALSE))


  return (clust_col)
}
