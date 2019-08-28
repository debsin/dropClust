# --------------------------------------------------
# 2D Projection
# --------------------------------------------------
#' Two Dimensional Embedding for Visualization
#' @description 2-D embedding of transcriptomes for visualization
#' @details Performs 2D embedding on sampled cells followed by the post-hoc projection of 2D co-ordinated of unsampled cells.
#' An ad-hoc clean up is also carried out to keep away stray points away from the plot - NA value is assigned to such points.
#' @param object A SingleCellExperiment object containing normalized expression values in \code{"normcounts"}.
#' @param embedding Type of embedding to use, one of \code{c("tsne", "umap")}. defaults to "tsne".
#' The \code{uwot} implementation of \code{umap} is used.
#' Specific arguments may be passed with relevant umap argument names.
#' @param use.subsamples logical indicating if sub-samples should be used for the post-hoc projection of
#' the coordinates of remaining samples.
#' @param ... \code{umap} arguments may be passed.
#' @return A SingleCellExperiment object with a new entry under the
#' \code{reducedDim()} containter to store the reduced dimension components
#' with the name from the argument \code{embedding}.
#' @importFrom methods is
#' @importFrom SingleCellExperiment colData reducedDim<- normcounts
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
#' sce <- PlotEmbedding(sce)
#' plot(reducedDim(sce, "tsne"))
PlotEmbedding<-function(object,embedding="tsne", use.subsamples=TRUE, ...){

  method.l = tolower(embedding)

  if(use.subsamples==TRUE && !is.null(colData(object)$Sampling)){
    # if(is.null(colData(object)$Sampling))
    #   stop("Subsamples not found.")
    subsamples_louvain = SingleCellExperiment::colData(object)$Sampling
  } else{
    subsamples_louvain = seq(ncol(object))
    use.subsamples = FALSE
  }

  if(!is.null(SingleCellExperiment::rowData(object)$PCAGenes))
    select_genes = SingleCellExperiment::rowData(object)$PCAGenes
  else
    select_genes = SingleCellExperiment::rowData(object)$HVG

  s_mat = dropClust::Log2Normalize(normcounts(object)[select_genes,subsamples_louvain], return.sparse = FALSE)


  # s_mat = LogNormalize(assay(object, "normcounts")[select_genes,subsamples_louvain],return.sparse = FALSE)
  # s_mat = mat[,subsamples_louvain]
  # s_mat = as.matrix(assay(object, "logcounts")[select_genes,subsamples_louvain])

  switch(method.l,
         tsne={
           cat("tSNE Embedding...")

           embeddings = Rtsne::Rtsne(t(s_mat),
                                     perplexity = 15,
                                     dims = 2,
                                     check_duplicates = FALSE)$Y
           cat("Done.\n")
         },
         umap={
           cat("UMAP Embedding...")
           embeddings = uwot::umap(t(s_mat), ...)
           cat("Done.\n")
         },

         {
           # default clutering
           stop(paste("Method",embedding,"not available."))
         }
  )

  if(use.subsamples){
    cat("Post-hoc co-ordinate assignment...")
    nn_length = length(unlist(colData(object)$ANN[[1]]))

    INDEX  = matrix(unlist(colData(object)$ANN),ncol=nn_length,byrow = T)
    clust_col = colData(object)$ClusterIDs
    clust_col[clust_col==0]<-NA

    sample.labels = clust_col[subsamples_louvain]
    PROJ = matrix(NA, nrow=nrow(INDEX), ncol=2)

    # nnn = ncol(INDEX)

    proj.idx = 1:nrow(INDEX)

    for( i in proj.idx)
    {
      if(is.na(clust_col[i])){next}
      maj_id = clust_col[i]
      v = INDEX[i,]
      maj_col = base::which(sample.labels[v] == maj_id)
      # if(length(maj_col)/nnn > 0.5){
      # PROJ[i,] = apply(tt[INDEX[i,maj_col],],2,function(x) mean(x, na.rm = T))
      PROJ[i,1] = mean(embeddings[v[maj_col],1], na.rm = TRUE)
      PROJ[i,2] = mean(embeddings[v[maj_col],2], na.rm = TRUE)
      # }
    }

    reducedDim(object, embedding)<-PROJ
    cat("Done.\n")
  }
  else
    reducedDim(object, embedding)<-embeddings[,1:2]
  return(object)

}
