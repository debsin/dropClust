#' Clustering of samples
#' @description Performs clustering on sampled cells and Post-hoc Cluster Assignment.
#' @details Clustering is carried out in two alternate approaches on the sampled cells.
#' For the default setting or quick identification of the existing broad clusters,
#' a Louvain based partition is employed. Otherwise for fine-tuned clustering with outliers,
#' hierarchical clustering is used with \code{cutreeDynamic} for dendrogram cut. Also, Assigns cluster membership to unsampled cells by using cluster membership information of the nearest neighbours.
#' An approximate nearest neighbour graph is constructed out of the samples population using the \code{find_ann()} module.
#' Some cells are left un-assigned when its neighbour's cluster membership doesn't form a majority as specified by the \code{conf} parameter.
#' Unassigned cells (\code{NA}) are excluded in the plot or further downstream analysis.
#' @param object A SingleCellExperiment object containing normalized expression values in \code{"normcounts"}.
#' @param use.subsamples uses the indices obtained from the sampling method.
#' @param method  character, one of c("default", "hclust", "kmeans").,
#' In the default mode, louvain based partition is used.
#' When \code{hclust}, hierarchical clustering is used.
#' @param use.reduced.dims optional, when used, the name of the \code{reducedDim()} table to be used for clustering.
#' @param minClusterSize integer, specifies the size of the smallest cluster; works when \code{method = "hclust"}.
#' @param deepSplit integer, level of dendrogram split [0-4], higher value produces finer clusters; ignored when
#' \code{method = "hclust"}.
#' @param k_nn integers, specifies number of nearest neighbours, defaults to 10.
#' @param conf numeric [0-1], defines the expected confidence of majority for a consensus. Cells remain unassigned when majority is below \code{conf}.
#' @param use.previous optional, when \code{TRUE}, the clustering step is skipped and only
#' the post-hoc clusering is repeated with the new \code{conf}.
#' @param ... For the \code{kmeans} method option, the argument \code{center} must be passed which specifies the number of clusters.
#' @return List of:\cr
#' \enumerate{
#' \item \code{cluster.ident} vector cluster identifiers ranging from 1 to the number of clusters for respective data points.\cr
#' \item \code{nn.ids}  matrix, each row corresponds to a cell, whose columns depict cluster membership of its neighbours; as returned by the \code{find_ann()} module. \cr
#' }
#'
#' Unassigned samples are represented by\code{NA} values.
#' @importFrom methods is
#' @importFrom SingleCellExperiment reducedDim  normcounts colData
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
#' sce <- Cluster(sce,  use.subsamples=FALSE, conf=0.1)
Cluster<-function(object,
                  use.subsamples= TRUE,
                  method = "default",
                  use.reduced.dims = NULL,
                  minClusterSize = 20,
                  deepSplit = 3,
                  k_nn=10,
                  conf=0.75,
                  use.previous = FALSE, ...){

  invisible(gc())
  if(any(names(reducedDims(object))=="CComponents"))
    use.reduced.dims = "CComponents"

  # If using integrated embeddings for clustering
  if(!is.null(use.reduced.dims)){
    mat = reducedDim(object, use.reduced.dims)

    cat("Clustering on embedded dimensions...")
    ss_clusters<- .ssClustering(ss_sel_genes = mat,
                                method = method,
                                minClusterSize = minClusterSize,
                                deepSplit = deepSplit,
                                trees = NULL, ...)
    SummarizedExperiment::colData(object)$ClusterIDs<-as.factor(ss_clusters$labels)
    cat("Done.\n")

    prev.method = object@metadata$dropClust["Clustering"]
    object@metadata$dropClust["Clustering"] = method


    return(object)
  }

  prev.method = ifelse(is.na(object@metadata$dropClust["Clustering"]), "", object@metadata$dropClust["Clustering"])

  if(use.previous == FALSE || is.null(colData(object)$ANN) || !prev.method==method){

    if(use.subsamples){
      if(is.null(colData(object)$Sampling))
        stop("Subsamples not found.")
      subsamples_louvain = SingleCellExperiment::colData(object)$Sampling
    } else{
      subsamples_louvain = seq(ncol(object))
    }


      if(!is.null(SingleCellExperiment::rowData(object)$PCAGenes))
        select_genes = which(SingleCellExperiment::rowData(object)$PCAGenes==TRUE)
      else
        select_genes = SingleCellExperiment::rowData(object)$HVG

      mat = Log2Normalize(normcounts(object)[select_genes,],return.sparse = FALSE)


    s_mat = t(mat[,subsamples_louvain])
    cat(nrow(s_mat), "samples and", ncol(s_mat), "genes used for clustering.\n")

    ann_trees = .buildAnnTrees(s_mat)

    # Clustering on subsamples-------------------------------------------------
    ss_clusters<-.ssClustering(ss_sel_genes = s_mat,
                               method = method,
                               minClusterSize = minClusterSize,
                               deepSplit = deepSplit,
                               trees = ann_trees, ...)
    SummarizedExperiment::colData(object)$Sample_ClusterIDs<-NA
    SummarizedExperiment::colData(object)$Sample_ClusterIDs[subsamples_louvain] <- ss_clusters$labels


    # Find K Nearest Neighbours among sub-samples ----
    cat("Find nearest neighbours among sub-samples...")
    SummarizedExperiment::colData(object)$ANN <- NA
    n_nn = k_nn + 1
    INDEX = mat.or.vec(ncol(object),n_nn)

    for(i in seq(ncol(object))){
      v = mat[,i]
      INDEX[i,] = ann_trees$getNNsByVector(v, n_nn)
    }
    cat("Done.\n")
    SummarizedExperiment::colData(object)$ANN <- apply(INDEX,1, list)
  }
  else{
    nn_length = length(unlist(colData(object)$ANN[[1]]))
    INDEX  = matrix(unlist(colData(object)$ANN),ncol=nn_length,byrow = T)
  }
  ss_labels = colData(object)$Sample_ClusterIDs[subsamples_louvain]

  # INDEX <- find_ann(object,graph = graph_samples, n_nn = k_nn)

  # UNANNOTATED Class Assignment----

  cat("Post-hoc Cluster Assignment...")
  clust_col <- .cluster_assign(INDEX,ss_labels,conf)
  cat("Done.\n")

  cat(paste("Unassigned Cells",length(clust_col[is.na(clust_col)])),"\n")

  residue = factor(clust_col[subsamples_louvain],exclude = NA)
  keep_labels = which(clust_col %in% levels(residue))
  all.labels = clust_col
  all.labels[-keep_labels]<-NA


  lev = 1:length(unique(all.labels[stats::complete.cases(all.labels)]))
  cc.f<-factor(all.labels,exclude = NA)
  levels(cc.f)<-lev

  cat(paste("Number of Predicted Clusters:",length(levels(cc.f))),"\n")

  SummarizedExperiment::colData(object)$ClusterIDs <-  cc.f

  object@metadata$dropClust["Clustering"] = method

  # "cluster.ident" = cc.f,
  return(object)

}


.ssClustering<-function(ss_sel_genes, method,
                        minClusterSize,
                        deepSplit,
                        trees, ...){
  switch(method,
         hclust={
           cat("Perfom Hierarchical Clustering...", minClusterSize, deepSplit, "\n")
           d = rdist::rdist(ss_sel_genes, metric = "angular")
           hc<-stats::hclust(d,method = "average")
           hc_labs<-dynamicTreeCut::cutreeDynamic(dendro = hc, cutHeight = NULL,
                                                  minClusterSize = minClusterSize,
                                                  method = "hybrid",
                                                  deepSplit = deepSplit,
                                                  pamStage = TRUE,
                                                  distM = as.matrix(d),
                                                  verbose = 0)

           names(hc_labs) = NULL
           cluster_label = list("labels"= hc_labs)
           cat(length(unique(hc_labs)), "clusters...")
           cat("Done.\n")
         },

         kmeans={
           cat("Perfom k-means Clustering...")
           km = stats::kmeans(ss_sel_genes, ...)
           cluster_label = list("labels"= km$cluster)
           cat("Done.\n")
         },

         {
           # default clutering
           if(is.null(trees))
             trees = .buildAnnTrees(ss_sel_genes)
           cat("Louvain Partitioning...")
           indices = lapply(seq(nrow(ss_sel_genes)),function(x) trees$getNNsByItem(x, 10))
           G=igraph::simplify(igraph::graph_from_adj_list(indices,
                                                          mode="all",
                                                          duplicate = FALSE))

           partition = igraph::cluster_louvain(G)
           cluster_label = list("labels"=partition$membership)
           cat("Done.\n")
         }
  )
  return(cluster_label)
}



.cluster_assign<-function(INDEX, sample.labels, conf){

  nn = dim(INDEX)[2]
  clust_col = list()


  for( i in seq(nrow(INDEX)))
  {
    clust_col[[i]] = NA
    tab = table(sample.labels[INDEX[i,]])
    max_id = which.max(tab)

    if(tab[max_id]/nn > conf){
      finally = clust_col[[i]]=names(max_id)
    }

  }


  clust_col = as.numeric(unlist(clust_col, use.names = FALSE))



  return (clust_col)
}




.buildAnnTrees<-function(data){

  dim(data)
  # set.seed(0)
  f = ncol(data)
  t = methods::new(RcppAnnoy::AnnoyAngular,f)  # Length of item vector that will be indexed

  cat("Build Graph with",nrow(data),"samples...")

  for(i in seq(nrow(data))){
    v = as.vector(data[i,])
    t$addItem(i, v)
  }
  t$build(30)# 30 trees
  cat("Done.\n")


  return(t)
}
