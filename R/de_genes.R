
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
#' @param log.norm.mat numeric (or character) vector having length same as the number of rows in the normalized matrix.
#' @param clust.list list object as returned by \code{cluster.cells} module.
#' @param sample_length integer to specify a maximum length of samples from each cluster.
#' @return integer vector of row  identifiers
#' @export
reduce_mat_de<-function(log.norm.mat, clust.list, sample_length){
  # Filter unassigned transciptomes.
  select_sample_ids = which(is.na(clust.list$cluster.ident)==FALSE)

  # Normalize by umi counts (median)
  norm_mat<-as.matrix(log.norm.mat$m[select_sample_ids,])
  colnames(norm_mat)<-log.norm.mat$use_genes
  rownames(norm_mat)<-rownames(log.norm.mat$m)[select_sample_ids]
  dim(norm_mat)

  # Read predicted cluster IDs
  pred_labels = clust.list$cluster.ident[select_sample_ids]

  # build sub-matrix for DE gene computation and heatmap
  fixed_samples <-find.cluster.samples(pred_labels)

  label = pred_labels[fixed_samples]
  mat_samples = norm_mat[fixed_samples,]
  return(list(mat_samples = mat_samples, labels = label))
}


# --------------------------------------------------
# Find DE Genes
# --------------------------------------------------
#' Find DE genes
#' @description Find cluster specific differentially expressed genes
#' @details Performs significance tests in one-vs-all mnner to determine cluster specific genes: two-sample Wilcoxon ('Mann-Whitney') test followed by 'fdr' adjustment and log-fold-change.
#' @param de_data list containing (1) matrix subset, each row corresponds to a transcriptome sample, columns represent genes; (2) predicted labels of the samples in matrix
#' @param selected_clusters vector of selected cluster identifer to be considered. When unspecified (=NA), defaults to all predicted clusters.
#' @param lfc_th numeric, [0,1] log2 fold change threshold value to define significance.
#' @param q_th numeric, [0,1] fdr adjusted p-value to define significance.
#' @param nDE integer, specifies the number of DE genes to return per cluster.
#' @return list containing:\cr
#' \itemize{
#' \item \code{DE_up} vector of unique DE gene names\cr
#' \item \code{DE_res} list of data frames for respective cluster specific DE genes, each row of a data frame mentions the gene significant pvalues, qvalues and LFC values}
#' @export
DE_genes <- function(de_data,selected_clusters=NA,lfc_th,q_th,nDE=30 )
{



  raw_data = de_data$mat_samples
  labels = de_data$labels




  if(any(is.na(selected_clusters))==TRUE)
    selected_clusters = unique(labels)


  if(!all(selected_clusters %in% labels))
    stop(paste("Cluster Ids must be among:", paste(unique(labels),collapse = " ")))

  cell.ids = which(labels %in% selected_clusters)


  DE_list = list()

  cat(paste("\nComputing for", length(unique(selected_clusters)),"clusters:\n"))

  for(i in unique(selected_clusters)){

    IND_a = which(labels == i)
    #IND_b = setdiff(1:ncol(data),IND_a)
    IND_b = which(labels != i)
    cat(paste("\nCluster",i,":\n"))

    # getting wilcoxon p values
    cat("\nComputing Wilcoxon p values...\n")
    PVAL<-apply(raw_data,2,function(x) stats::wilcox.test(x[IND_a],
                                                          x[IND_b])$p.value)

    # fdr
    fdr <- stats::p.adjust(PVAL,method="fdr")

    # prepare final result
    DE_res <- data.frame(cbind(pvalues=PVAL,qvalues=fdr))
    rownames(DE_res) <- colnames(raw_data)


    # log fold change
    cat("Computing Log fold change values...\n")
    LFC = log2(Matrix::colMeans(raw_data[IND_a,])/
                 Matrix::colMeans(raw_data[IND_b,]))


    cat("\nCompleted successfully.\n")


    sig = DE_res[intersect(which(LFC>0),
                           intersect(which(abs(LFC)>=lfc_th),
                                     which(DE_res$qvalues<=q_th))),]
    # sig = DE_res[intersect(which(abs(LFC)>=lfc_th),
    #                        which(DE_res$qvalues<=q_th)),]


    rN = rownames(sig)[order(sig$qvalues)]


    DE_list[[paste0(i)]] =
      data.frame(gene = rN,
                 q_val = DE_res$qvalues[match(rN,rownames(DE_res))],
                 fc = LFC[match(rN,names(LFC))])
  }

  # names(DE_list) = unique(labels)


  DE_up<-top.de.genes(DE_list, nDE)

  genes_df <-list_to_df(DE_up)
  colnames(genes_df)<-paste("cluster",names(DE_up),sep='_')

  RES = list(genes.df = genes_df, DE_res = DE_list)
  return(RES)
}


# --------------------------------------------------
# Convert List to Data Columns
# --------------------------------------------------
list_to_df<-function(l){
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
top.de.genes<-function(l,nDE=30){
  DE_up=list()
  for(type in names(l)){
    ordered = order(l[[type]]$q_val,abs(l[[type]]$fc))
    row = as.character(utils::head(l[[type]]$gene[ordered],nDE))
    DE_up[[type]]= row
  }

  # all_ct_genes = unique(unlist(DE_up))
  names(DE_up)<-names(l)
  return(DE_up)
}
