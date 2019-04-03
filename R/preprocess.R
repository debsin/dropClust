# --------------------------------------------------
# Batch-Effect Removal
# --------------------------------------------------
#' Batch-Effect Removal
#' @description Minimise batch effect by representing normalized expression values with suitable ranks
#' @details Minimising batch effect relies on an assumption that relative ordering of genes with
#' regards to their expressions are preserved across batches.
#' @param amatrix matrix, each row corresponds to a cell, whose columns depict cluster membership of its neighbours; as returned by \code{find_ann()} module.
#' @param close_th numeric, [0,1] percentage of proximity between two genes.
#' Two genes are considered close when \code{(v2 - v1)/v1 < close_thresh}.
#' @param cells_th numeric, [0,1] percentage of cells which need to be close for one of the genes to be removed.
#' @return numeric matrix where expression values are replaced with ranks.
#' @export
batcheffect.removal<-function(amatrix, close_th = 0.01, cells_th = 0.05){
  ### Gene Removal for Ranking
  cat("Batch correcting...\n")
  #removes gene-columns and returns a reduced expression matrix
  #alpha marks the percentage of cells which need to be close for one of the genes to be removed
  input <- list(matrix  = as.matrix(amatrix),
                close_threshold = close_th,
                cells_threshold = cells_th,
                nreduce = 1000)
  remove_these <- reduce_genes_cpp(input)

  # print(remove_these)
  reduced_matrix <- amatrix[, which(remove_these == 0)]

  ### Converting the given matrix to rank vectors
  rank_mat <- get_rank_mat(as.matrix(reduced_matrix))
  dim(rank_mat)

  whole = Matrix::Matrix(rank_mat, sparse = TRUE)
  return(whole)
}



### Converting the given matrix to rank vectors
get_rank_mat <- function(matrix) {

  cat("Constructing Rank matrix...\n")
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
    sorted_genes <- order(unlist(matrix[iter, ]))
    # sorted_genes <- grr::order2(unlist(matrix[iter, ]))
    # Assigns the rank to each gene and modifies the corresponding vector in the new dataframe
    new_matrix[iter,]<-match(colnames(matrix),colnames(matrix)[sorted_genes])
    # for (ind in 1:length(sorted_genes)) {
    #   new_matrix[iter, sorted_genes[ind]] <- ind
    # }
  }
  # Return the rank-vectors dataframe
  return(new_matrix)
}




#' Rare Genes and Cells
#' @description Identifying potential genes co-expressed in only a minor set of cell populations. This helps to boost rare cell discovery.
#' @details Highly co-expressed genes (between 5-250 genes) are searched for within only a small percentage of cell populations (between 10-200 cells).
#' The searched is performed recursively until converged.
#' @param x named list containing 3 objects: \cr
#' \enumerate{
#' \item \code{mat} sparse numeric matrix with count expression data, with genes in columns. \cr
#' \item \code{gene_symbols} string vector with unique gene names for each column in \code{mat}. \cr
#' \item \code{barcodes} string vector of length equal to number of rows in \code{mat} identifying each transcriptome per row . \cr
#' }
#' @return named list containing two character vectors: gene names in \code{rare_genes} and barcodes in \code{rare_cells}.
#' @import RcppEigen
#' @export
get_rare_genes <-function(x) {

  cat("Finding co-expressed rare genes...\n")
  th=0.0001
  tol=99
  mat = x$mat
  rare_genes = x$gene_symbols
  rare_cells = x$barcodes
  iter = 0
  while(tol>th && stats::complete.cases(tol))
  {
    iter = iter+1

    cs <- Matrix::colSums(mat>=5)
    genes_filtered = intersect(which(cs>=10), which(cs<=200))

    rs <- Matrix::rowSums(mat[,genes_filtered]>=5)
    cells_with_coex_gene = intersect(which(rs>=5), which(rs<=250))

    gene_mat = mat[cells_with_coex_gene,genes_filtered]

    tol = (nrow(mat)-nrow(gene_mat))/nrow(mat)

    rare_cells = rare_cells[cells_with_coex_gene]
    rare_genes = rare_genes[genes_filtered]
    mat = gene_mat
  }
  # dim(gene_mat)

  return(list("rare_genes"=rare_genes,"rare_cells"=rare_cells))
}



#' Cell Filter
#' @description Filter poor quality cells.
#' @details Keep only those cells expressing atleast count = \code{min_count} in the number of genes specified within the quantile range between \code{ql_th} and \code{qh_th}.
#' @param data named list containing 3 objects: \cr
#' \enumerate{
#' \item \code{mat} sparse numeric matrix with count expression data, with genes in columns. \cr
#' \item \code{gene_symbols} string vector with unique gene names for each column in \code{mat}. \cr
#' \item \code{barcodes} string vector of length equal to number of rows in \code{mat} identifying each transcriptome per row. \cr
#' }
#' @param min_count integer threshold for expression count
#' @param ql_th quantile at probability with values in [0,1] for lower limit
#' @param qh_th quantile at probability with values in [0,1] for upper limit
#' @return updated input \code{data} list appended with character vector named \code{keep_cells} containing indices of good cells.
#' @export
filter_cells<-function(data, min_count=3, ql_th = 0.001, qh_th = 1){

  try(if(qh_th < ql_th) stop("Lower limit greater than upper limit!!"))

  rs<-Matrix::rowSums(data$mat>=min_count)

  l1 = stats::quantile(rs,probs = ql_th)
  l2 = stats::quantile(rs,probs = qh_th)


  keep_cells = intersect(which(rs>=l1), which(rs<=l2))
  #Return good cells based on IQR rule
  #q<-quantile(rs)
  #iqr<-IQR(rs)

  #l1 = q[2] - 1.5*(iqr)
  #l2 = q[4] + 1.5*(iqr)

  #keep_cells = intersect(which(rs > l1), which(rs <l2))
  cat(paste(length(rs)-length(keep_cells), "bad cells present.\n"))
  data$mat = data$mat[keep_cells,]
  data$barcodes = data$barcodes[keep_cells]
  keep_cells = list(keep_cells)
  names(keep_cells)="keep_cells"
  data = append(data, keep_cells)
  return(data)
}


#' UMI Count Normalization with Gene Filter
#' @description Matrix normalization by UMI count, also remove poorly expressed genes.
#' @details Normalizes count value using: \code{total count/median} of respective cells.
#' Genes with count greater than \code{min.count} and not expressed in atleast \code{min.cell} cells are removed.
#' @param x named list containing 3 objects: \cr
#' \enumerate{
#' \item \code{mat} sparse numeric matrix with count expression data, with genes in columns. \cr
#' \item \code{gene_symbols} string vector with unique gene names for each column in \code{mat}. \cr
#' \item \code{barcodes} string vector of length equal to number of rows in \code{mat} idendifying each transcriptome per row. \cr
#' }
#' @param min.count integer threshold for expression count
#' @param min.cell integer threshold for number of cells expressing a particular gene
#' @return list containing normalized matrix \code{m} and a character vector named \code{use_genes} of length = dim(m)[2] containing the gene names .
#' @export
normalize<-function(x, min.count=2, min.cell=3) {

  cat("Discard poor genes...\n")
  mat  = x$mat
  cs <- Matrix::colSums(mat>min.count)
  x_use_genes <- which(cs > min.cell)

  x_filt<-mat[,x_use_genes]
  gene_symbols = x$gene_symbols[x_use_genes]
  cat("Dimensions of filtered Matrix: ")
  cat(dim(x_filt),"\n")

  cat("Median normalize matrix...\n")

  rs<-Matrix::rowSums(x_filt)
  rs_med<-stats::median(rs)
  x_norm<-x_filt/(rs/rs_med)

  list(m=x_norm,use_genes=gene_symbols)
}



#' Top Dispersed Genes
#' @description Get variable genes from normalized UMI counts using Fano Factor metric.
#' @details Compute Fano Factor metric for each gene. The metric computes the median absolute deviation of dispersion across multiple bins for each gene.
#' @param normalized_data object returned by \code{normalize()} module.
#' @param ngenes_keep integer to return top ranking \code{ngenes_keep} number of genes.
#' @return vector of gene names of length \code{ngenes_keep} .
#' @export
dispersion_genes<-function(normalized_data, ngenes_keep=1000) {
  avg<-Matrix::colMeans(normalized_data$m)
  dispersion = ColDispersion(normalized_data$m)

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


  cat("Sort Top Genes...\n")
  disp_cut_off<-sort(df$dispersion_norm,decreasing=TRUE)[ngenes_keep]
  cat("Cutoff Genes...\n")
  df$used<-df$dispersion_norm >= disp_cut_off
  features = utils::head(order(-df$dispersion_norm),ngenes_keep)
  dispersed_genes = normalized_data$use_genes[features]

  return(dispersed_genes)
}




#' Subset Normalized Matrix
#' @description Subset Normalized matrix to include highly dispersed genes and suspected rare cell genes.
#' @details Subset input normalized matrix to include only \code{dp_genes} and \code{rare_genes}.
#' @param normalized_data object returned by \code{normalize()} module.
#' @param dp_genes gene names returned by \code{dispersion_genes()} module.
#' @param rare_data optional, return object obtained from the \code{get_rare_genes()} module.
#' @param batch_correction logical, when \code{TRUE} performs batch correction with default parameters.
#' @return numeric matrix.
#' @export
matrix.transform<-function(normalized_data,
                           dp_genes,
                           rare_data=NULL,
                           batch_correction=FALSE){

  if(is.null(rare_data)){
    rare_genes=c()
  }else{
    rare_genes= rare_data$rare_genes
  }

  dp_ids = list()
  for(i in 1:length(dp_genes)){
    dp_ids[[i]] = which(normalized_data$use_genes==dp_genes[i])
  }
  dp_ids=unlist(dp_ids, use.names = FALSE)

  #
  rare_feature_ids = which(normalized_data$use_genes %in% rare_genes)


  #Final Data
  m_n_whole<-normalized_data$m[,union(dp_ids,rare_feature_ids)]
  m_filt<-Matrix::Matrix(log2(m_n_whole+1),sparse = TRUE)

  colnames(m_filt) = normalized_data$use_genes[union(dp_ids,rare_feature_ids)]

  if(batch_correction==TRUE){
    m_filt<-batcheffect.removal(m_filt)
  }

  return(m_filt)

}

