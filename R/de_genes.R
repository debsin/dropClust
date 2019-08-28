# --------------------------------------------------
# Sample from each cluster
# --------------------------------------------------

find.cluster.samples<-function(labels,sample_length=100){
  # set.seed(0)
  fixed_samples = c()
  for(clust_id in levels(as.factor(labels))){
    if(clust_id==0) next;
    ids = which(labels==clust_id)
    fixed_samples = c(fixed_samples,
                      sample(ids,min(sample_length,length(ids))))
  }

  return(fixed_samples)

}
# --------------------------------------------------
# Reduce large matrix for DE analysis
# --------------------------------------------------
#' Reduce large matrix for DE analysis
#' @description Construct submatrix by sampling from each cluster
#' @param norm.mat numeric (or character) vector having length same as the number of rows in the normalized matrix.
#' @param clust.list list object as returned by \code{cluster.cells} module.
#' @param sample_length integer to specify a maximum length of samples from each cluster.
#' @return integer vector of row  identifiers
#' @export
.reduce_mat_de<-function(norm.mat, clust.list, sample_length=100){
  # Filter unassigned transcriptomes.
  select_sample_ids = which(is.na(clust.list)==FALSE)

  # Normalize by umi counts (median)
  norm_mat<-as.matrix(norm.mat[select_sample_ids,])
  colnames(norm_mat)<-colnames(norm.mat)
  rownames(norm_mat)<-rownames(norm.mat)[select_sample_ids]
  dim(norm_mat)

  # Read predicted cluster IDs
  pred_labels = clust.list[select_sample_ids]

  # build sub-matrix for DE gene computation and heatmap
  fixed_samples <-find.cluster.samples(pred_labels, sample_length)

  label = pred_labels[fixed_samples]
  mat_samples = norm_mat[fixed_samples,]
  return(list(mat_samples = mat_samples, labels = label))
}


# --------------------------------------------------
# Find DE Genes
# --------------------------------------------------
#' Find DE genes
#' @description Find cluster specific differentially expressed genes
#' @details Performs significance tests in one-vs-all manner to determine cluster specific genes: two-sample Wilcoxon ('Mann-Whitney') test followed by 'fdr' adjustment and log-fold-change.
#' @param object A SingleCellExperiment object containing normalized expression values in \code{"normcounts"}.
#' @param selected_clusters vector of selected cluster identifier to be considered. When unspecified (=NA), defaults to all predicted clusters.
#' @param lfc_th numeric, [0,1] log2 fold change threshold value to define significance.
#' @param q_th numeric, [0,1] fdr adjusted p-value to define significance.
#' @param nDE integer, specifies the number of DE genes to return per cluster default = 30.
#' @importFrom SingleCellExperiment normcounts rowData
#' @return list containing:\cr
#' \itemize{
#' \item \code{DE_up} vector of unique DE gene names\cr
#' \item \code{DE_res} list of data frames for respective cluster specific DE genes, each row of a data frame mentions the gene significant pvalues, qvalues and LFC values}
#' @export
FindMarkers <- function(object, selected_clusters=NA, lfc_th, q_th, nDE=30)
{
  if(is.null(object$ClusterIDs)) stop("ClusterIDs not found.")

  if(is.null(rowData(object)$Symbol))
    SummarizedExperiment::rowData(object)$Symbol = rownames(object)

  de_data<- .reduce_mat_de(Matrix::t(normcounts(object)),object$ClusterIDs)

  raw_data = t(de_data$mat_samples)

  rownames(raw_data) = rowData(object)$Symbol
  labels = de_data$labels

  if(any(is.na(selected_clusters)))
    selected_clusters = unique(labels)


  if(!all(selected_clusters %in% labels))
    stop(paste("Cluster Ids must be among:", paste(unique(labels),collapse = " ")))

  cell.ids = which(labels %in% selected_clusters)


  DE_list = list()

  cat("Computing for DE genes:\n")

  for(i in unique(selected_clusters)){

    IND_a = which(labels == i)
    IND_b = which(labels != i)
    cat(paste("Cluster",i,":\n"))

    # getting wilcoxon p values
    cat("\tComputing Wilcoxon p values...")
    PVAL<-apply(raw_data,1,function(x) stats::wilcox.test(x[IND_a], x[IND_b])$p.value)
    # fdr
    QVAL <- stats::p.adjust(PVAL,method="fdr")
    cat("Done.\n")
    # prepare final result
    DE_res <- data.frame(cbind("pvalues"=PVAL,"qvalues"=QVAL))

    rownames(DE_res) <- rownames(raw_data)


    # log fold change
    cat("\tComputing Log fold change values...")
    LFC = log2(Matrix::rowMeans(raw_data[,IND_a])/
                 Matrix::rowMeans(raw_data[,IND_b]))


    cat("Done.\n")


    sig = DE_res[intersect(which(LFC > 0),
                           intersect(which(abs(LFC) >= lfc_th), which(DE_res$qvalues <= q_th))), ]

    # sig = DE_res[intersect(which(abs(LFC)>=lfc_th),
    #                        which(DE_res$qvalues<=q_th)),]


    rN = rownames(sig)[order(sig$qvalues)]


    DE_list[[paste0(i)]] =
      data.frame(gene = rN,
                 qvalues = DE_res$qvalues[match(rN,rownames(DE_res))],
                 fc = LFC[match(rN,rownames(LFC))])
  }

  # names(DE_list) = unique(labels)


  DE_up<-.top.de.genes(DE_list, nDE)

  genes_df <-.listToDF(DE_up)
  colnames(genes_df)<-paste("cluster",names(DE_up),sep='_')

  RES = list(genes.df = genes_df, DE_res = DE_list)
  return(RES)
}


# --------------------------------------------------
# Convert List to Data Columns
# --------------------------------------------------
.listToDF<-function(l){
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  mat <- as.data.frame(sapply(l, "[", i = seq.max))
  gene.df <- sapply(mat, as.character)
  gene.df[is.na(gene.df)] <- ""
  return(gene.df)
}



# --------------------------------------------------
# Fetch top DE genes
# --------------------------------------------------
.top.de.genes<-function(l,nDE=30){
  DE_up=list()
  for(type in names(l)){
    ordered = order(l[[type]]$qvalues,abs(l[[type]]$fc))
    row = as.character(utils::head(l[[type]]$gene[ordered],nDE))
    DE_up[[type]]= row
  }

  # all_ct_genes = unique(unlist(DE_up))
  names(DE_up)<-names(l)
  return(DE_up)
}
