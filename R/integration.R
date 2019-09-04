#' Correct for batch effects
#' @description Correct the merged count data based on rank values to obtain a set of reduced and corrected dimensions.
#' @details Concatenate the expression counts of all cells from different batches into one expression count object.
#' The merging is done on the set of union of DE genes obtained from the clustering of each batch.
#' @param object A list of SingleCellExperiment objects, each representing a  SingleCellExperiment object from a single batch.
#' @param method character, one of c("default","fastmnn").
#' \code{default} mode performs the dropClust based correction followed by UMAP embedding.
#' The \code{fastmnn} option performs the mutual neighbourhood based correction which is implemented in the \code{batchelor} package.
#' when \code{FALSE} the batches are merged on the set of common genes across batches.
#' @return  A SingleCellExperiment object with two new entry under the
#' \code{reducedDim()} container to store the reduced dimension components
#' with the name \code{"CComponents"} and the rank expression matrix named \code{"RankMat"}.
#' @param close_th for the method = default, specifies the value at which the
#' expression values of two genes will be considered as close pairs.
#' @param cells_th for the method default, specifies the value to
#' determine what proportion of total number of cells have close pairs.
#' @param components number of reduced dimensions to return.
#' @param ... \code{umap} arguments may be passed.
#' @importFrom SingleCellExperiment reducedDim normcounts
#' @export
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' ncells <- 100
#' ngenes <- 1200
#' lambda <-abs(rnorm(ngenes))
#' counts.1 <- matrix(rpois(ncells*ngenes, lambda =  lambda), ncol=ncells, nrow=ngenes, byrow=TRUE)
#' rownames(counts.1) <- paste0("Gene", seq_len(ngenes))
#' colnames(counts.1) <- paste0("Cell", seq_len(ncells))
#' sce.1 <- SingleCellExperiment(assays = list(counts = counts.1))
#' rowData(sce.1)$Symbol <- paste0("Gene", seq_len(ngenes))
#'
#' lambda <-abs(rnorm(ngenes))
#' counts.2 <- matrix(rpois(ncells*ngenes, lambda =  lambda), ncol=ncells, nrow=ngenes, byrow=TRUE)
#' rownames(counts.2) <- paste0("Gene", seq_len(ngenes))
#' colnames(counts.2) <- paste0("Cell", seq_len(ncells))
#' sce.2 <- SingleCellExperiment(assays = list(counts = counts.2))
#' rowData(sce.2)$Symbol <- paste0("Gene", seq_len(ngenes))
#'
#' mixed_sce <- Merge(list(sce.1, sce.2), use.de.genes =TRUE)
#' mixed_sce <- Correction(mixed_sce, close_th=0.1, cells_th=0.2)
#' }
Correction<-function(object, method="default", close_th=0.1, cells_th=0.1, components=2, ...){

  batchids = object$Batch
  temp<- CountNormalize(object)

  if(method=="fastmnn"){
    mat =  Log2Normalize(normcounts(temp))
    mnn.obj  = batchelor::fastMNN(mat,batch = batchids,  k=20, pc.input = FALSE,...)
    reducedDim(temp,"CComponents") = reducedDim(mnn.obj,"corrected")
    SummarizedExperiment::assay(temp, "reconstructed")<-SummarizedExperiment::assay(mnn.obj, "reconstructed")
  }  else{
    if(!is.null(rowData(temp)$CommonDEGenes))
      mat = Log2Normalize(normcounts(temp)[rowData(temp)$CommonDEGenes,])
    else
      mat = Log2Normalize(normcounts(temp))
    corrected = as.matrix(.batcheffect.dc(t(mat), close_th = close_th, cells_th = cells_th))

    cat("Embedding with UMAP...")
    if(!is.null(temp$Sampling)){
      umap_model = uwot::umap(corrected[temp$Sampling,], metric = 'cosine', n_components = components, ret_model=TRUE, ...)
      umap_proj_dc = uwot::umap_transform(corrected, umap_model)
    } else
      umap_proj_dc = uwot::umap(corrected, metric = 'cosine', n_components = components, ...)
    cat("Done\n")

    reducedDim(temp, "RankMat")<-corrected
    reducedDim(temp, "CComponents")<-umap_proj_dc
  }


  object@metadata[["dropClust"]] = c(object@metadata[["dropClust"]], "Correction")

  return(temp)
}


.batcheffect.dc<-function(amatrix, close_th = 0.01, cells_th = 0.05){
  ### Gene Removal for Ranking
  cat("Batch correcting...\n")
  #removes gene-columns and returns a reduced expression matrix
  #alpha marks the percentage of cells which need to be close for one of the genes to be removed
  input <- list(matrix  = as.matrix(amatrix),
                close_threshold = close_th,
                cells_threshold = cells_th,
                nreduce = 1000)
  remove_these <- reduce_genes_cpp(input)

  reduced_matrix <- amatrix[, which(remove_these==0)]
  cat("from",ncol(amatrix), "to",ncol(reduced_matrix),"genes.\n" )

  # mute genes
  # muted_matrix <- .mutegenes(reduced_matrix)

  ### Converting the given matrix to rank vectors
  rank_mat <- .get_rank_mat_r(as.matrix(reduced_matrix))
  dim(rank_mat)

  whole = Matrix::Matrix(rank_mat, sparse = TRUE)
  return(whole)
}

#
# .mutegenes<-function(matrix){
#
# }

#

### Converting the given matrix to rank vectors
.get_rank_mat_r <- function(matrix) {

  # cat("Constructing Rank matrix...\n")
  # Gets the gene names from the data frame
  columns <- names(matrix)
  # Creates a new data frame the same no of rows and columns
  # This data frame will be returned
  new_matrix <- matrix(nrow = dim(matrix)[1], ncol = dim(matrix)[2])
  # Setting the columns of the new dataframe as the genes from the original matrix
  # names(new_matrix) <- columns

  # Iterating over the matrix and getting the rank vector for each cell
  for (iter in 1:dim(matrix)[1]) {
    # Gives the names of the genes in sorted order based on the gene-expression value
    um_vec = matrix[iter,]
    # nth_max_gene = order(um_vec,decreasing = T)[100]
    # um_vec[which(um_vec<nth_max_gene)] = 0

    new_matrix[iter,]<-rank(replace(um_vec, um_vec==0, NA), na.last = 'keep' )

  }
  new_matrix[is.na(new_matrix)] <- 0
  # Return the rank-vectors dataframe
  return(new_matrix)
}
#
# reduce_genes_r<-function(input){
#   mat = input$matrix
#   close_threshold = input$close_threshold
#   cells_threshold = input$cells_threshold
#
#   print(dim(mat))
#   print(close_threshold)
#   print(cells_threshold)
#   ncells = nrow(mat)
#   ngenes = ncol(mat)
#
#   vote_matrix = mat.or.vec(nr = ngenes,nc = ngenes)
#
#   for(i in 1:ngenes){
#     for(j in 1:ngenes){
#       delta_ratios = abs(mat[,i]-mat[,j])/pmax(mat[,i],mat[j])
#       print(delta_ratios)
#       if(sum(delta_ratios<close_threshold) > (cells_threshold * ncells)) {
#         vote_matrix[i,j] = vote_matrix[i,j] - 1
#       }
#       else
#         vote_matrix[i,j] = vote_matrix[i,j] + 1
#     }
#   }
#
#   rs = rowSums(vote_matrix)
#   ordered = order(rs,decreasing = T)
#   qualified = setdiff(ordered, which(seq(ngenes)>0))
#
#   return(qualified)
# }

