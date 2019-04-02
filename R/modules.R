#' Intitial 1/3rd sampling if very big dataset
#' @param data_list list containing \code{read10x()} like structure
#' @param rare_data optional, return object obtained from the \code{get_rare_genes()} module.
#' @return two column integer matrix, first column represents sample id, second column contain cluster membership id
#' @export
initial.samples<-function(data_list,rare_data=NULL){

  if(is.null(rare_data)){
    rare_cells=c()
  }else{
    rare_cells= rare_data$rare_cells
  }
  # set.seed(0)

  no_samples = dim(data_list$mat)[1]

  i = ifelse(no_samples>20000, min(20000,round(no_samples/3)), no_samples)
  # i=min(20000, round(no_samples/3))
  sample_ids = sample(1:no_samples, i)

  sample_ids = union(sample_ids,which(data_list$barcodes %in% rare_cells))

  return(sample_ids)
}



#' ANN Graph Partition
#' @param data numeric matrix with \code{n} samples in rows
#' @return two column integer matrix, first column represents sample id, second column contain cluster membership id
#' @export
ann_partition<-function(data){
  dim(data)
  # set.seed(0)
  f = dim(data)[2]
  cat("Build graph...\n")
  t = methods::new(RcppAnnoy::AnnoyAngular,f)  # Length of item vector that will be indexed
  for(i in seq(nrow(data))){
    v = data[i,]
    t$addItem(i, v)
  }
  t$build(30)# 30 trees

  indices = lapply(seq(nrow(data)),function(x) t$getNNsByItem(x, 6))

  G=igraph::simplify(igraph::graph_from_adj_list(indices,
                                                 mode="all",
                                                 duplicate = FALSE))
  cat("Louvain Partition...\n")

  partition = igraph::cluster_louvain(G)

  dataMatrix =as.data.frame(cbind(seq(nrow(data)),partition$membership))
  return(dataMatrix)
}


#' ANN Search
#' @param data numeric matrix with \code{n} samples in rows
#' @param idx integers to specify specific sample indices. The ANN graph is fit on these samples only.
#' @param n_nn integers, specifies number of nearest neighbours.
#' @return integer matrix of dimension \code{n x n_nn}
#' @export
find_ann<-function(data,idx, n_nn = 10){
  n_nn = n_nn+1
  sub_mat = data[idx,]
  dim(sub_mat)
  cat("Build Graph...\n")
  f = dim(sub_mat)[2]
  t = methods::new(RcppAnnoy::AnnoyAngular,f)  # Length of item vector that will be indexed
  # set.seed(100)

  for(i in seq(nrow(sub_mat))){
    v = sub_mat[i,]
    t$addItem(i, v)
  }
  t$build(30) # 30 trees
  cat("Find ANN...\n")
  indices = mat.or.vec(nrow(data),n_nn)

  for(i in seq(nrow(data))){
    indices[i,] = t$getNNsByVector(data[i,], n_nn)
  }
  return(indices)
}






#' Cluster Validity Index
#' @description Validate clusters with known annotations when available.
#' @details Computes cluster external validity metrics in terms of Rand Index, Adjusted Rand Index and Purity. \cr
#' Purity of a cluster is measured as a recall of the largest known cell type within a predicted cluster. Total purity across all clusters is returned.
#' @param predicted_clusters vector of predicted cluster IDs.
#' @param known_types vector of annotated IDs corresponding to \code{predicted_clusters} in the same order.
#' @param show_table logical, When TRUE, shows the contingency table of concordance between the IDs.
#' @return numeric a vector of three metrics.
#' @export
sc_metric<-function(predicted_clusters ,known_types ,show_table = TRUE){
  dm<-table(known_types,predicted_clusters)
  if(show_table==TRUE){print(dm)}
  r.i. = flexclust::randIndex(dm, correct=TRUE, original=TRUE)
  return(round(c(r.i.,
                 "Purity" = ClusterPurity(predicted_clusters,
                                          known_types)), 3))
}

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(clusters,classes), 2, max)) / length(clusters)
}



# ---------------------------------------------
# Optimize parameter 'pinit' for
# exponential decay sampling
# ---------------------------------------------
optimized_param<-function(partition, nsamples = 500){
  if(nsamples > nrow(partition)) stop("Expected number of samples exceeded total population.")
  oldseed = .Random.seed
  # set.seed(0)
  global.min <- 0
  tol <- 1e-3

  max.time= 20
  lower <- c(0.05, 0.9,500)
  upper <- c(0.1,0.95,4000)

  params<-NULL

  pin_find<-function(params, MAX_c = nsamples){
    output <- partition
    # oldseed = .Random.seed
    # set.seed(0)
    cluster_freq = table(output[,2])
    pinit = params[1]
    pfin = params[2]
    K=params[3]
    prop = round((pinit -
                    exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
    t(rbind(cluster_freq,prop))
    prop = reshape2::melt(prop)$value
    subsamples_louvain<-c()
    for(k in 1:length(prop)){
    subsamples_louvain = c(subsamples_louvain,
                             sample(output[which(output$V2==(k)),1],
                                    size = prop[k],replace = FALSE))
    }
    # .Random.seed = oldseed

    return(abs(MAX_c-length(subsamples_louvain)))
  }


  out <- GenSA::GenSA(par=params,fn = pin_find, lower = lower, upper = upper,
               control=list(max.time = max.time,
                            threshold.stop=global.min+tol,
                            verbose=TRUE))

  .Random.seed = oldseed

  #out[c("value","par","counts")]

  names(out$par) <-c("pinit", "pfin", "K")
  return(out$par)
}


# ---------------------------------------------
# Sampling Primary Clusters
# ---------------------------------------------
#' Sampling Primary Clusters
#' @description Performs sampling from the primary clusters in an inverse exponential order of cluster size.
#' @details Sampling in inverse proportion of cluster size following a exponential decay equation. To ensure selection of sufficient representative transcriptomes from small clusters, an exponential decay function  is used to determine the proportion of transciptomes to be sampled from each cluster. For $i^{th}$ cluster, the proportion of expression profiles $p_i$ was obtained as follows.\cr
#' \ifelse{html}{\out{p<sub>i</sub> = p<sub>l</sub> <font face="symbol">-</font> e<sup><font face="symbol">-</font>(S<sub>i</sub>)/(K)</sup>}}{\deqn{Lp_{i} = p_{l} - e^{-\frac{S_i}{K}} (p_{l} - p_{u})}{ASCII}}
#' where \eqn{S_i} is the size of cluster \eqn{i}, \eqn{K} is a scaling factor, \eqn{p_i} is the proportion of cells to be sampled from the $i^{th}$ Louvain cluster. $p_l$ and $p_u$ are lower and upper bounds of the proportion value respectively.
#' @references {
#' \insertRef{sengupta2013reformulated}{dropClust}
#' }
#' @param data numeric matrix with \code{n} samples in rows
#' @param optm_parameters logical, when TRUE the parameters (\code{pinit, pfin, K}) are optimized such that exactly \code{nsamples} are returned. Optimization is performed using simulated anneling
#' @param nsamples integer, total number of samples to return post sampling; ignored when \code{optm_parameters = FALSE}.
#' @param pinit numeric [0,0.5], minimum probability of that sampling occurs from a cluster, ignored when \code{optm_parameters = TRUE}.
#' @param pfin numeric [0.5,1], maximum probability of that sampling occurs from a cluster, ignored when \code{optm_parameters = TRUE}.
#' @param K numeric, scaling factor analogous to Boltzman constant, ignored when \code{optm_parameters = TRUE}.
#' @return numeric vector containing the indices of sampled cells.
#' @export
sampling<-function(data, optm_parameters=FALSE,
                   nsamples=500,
                   pinit=0.195,
                   pfin = 0.9,
                   K=500){

  partition<-ann_partition(as.matrix(data))

  output<-partition
  if(optm_parameters==TRUE){
    param = optimized_param(partition, nsamples)
    pinit = param[1]
    pfin = param[2]
    K = param[3]
    cat("Optimized parameters:\n", param,"\n")
  }

  oldseed = .Random.seed
  # set.seed(0)
  cluster_freq = table(partition[,2])
  prop = round((pinit -
                  exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
  t(rbind(cluster_freq,prop))
  prop = reshape2::melt(prop)$value
  subsamples_louvain<-c()
  for(k in 1:length(prop)){
    subsamples_louvain = c(subsamples_louvain,
                           sample(partition[which(partition[,2]==(k)),1],
                                  size = prop[k],replace = FALSE))
  }
  .Random.seed = oldseed
  return(subsamples_louvain)
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
#' @param mat normalized expression matrix containing the subset of samples. Subsetting the return matrix by \code{matrix.subset()} module with samples from primary clusters.
#' @param top integer specifying to number of genes to return in their order of ranking
#' @return integer vector containing the indices of \code{top} number of genes in the input matrix
#' @export
pc_genes<-function(mat,top=200){
  nc = 50
  cat("Find best PCA components...\n")
  PRR <- irlba::prcomp_irlba(as.matrix(mat),nc)
  S = c()
  for(i in 1:nc)
  {
    K = PRR$x[,i]
    # L = mclust::Mclust(K)
    # S = c(S,L$G)
    peaks = all_peaks(K)
    # modes= NbClust::NbClust(K, distance = "euclidean",
    #                min.nc = 2, max.nc = 10,
    #                method = "kmeans", index ="ball")$Best.nc[1]
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



  return(top_pc_genes)
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


# --------------------------------------------------
# 2D Projection
# --------------------------------------------------
#' Two Dimensional Embedding for Visualization
#' @description 2-D embedding of transcriptomes for visualization
#' @details Perms t-SNE for 2D embedding on sampled cells followed by the post-hoc projection of 2D co-ordinated of unsampled cells.
#' An ad-hoc clean up is also carried out to keep away stray points away from the plot - NA value is assigned to such points.
#' @param data normalized expression matrix containing the subset of samples and PCA based selected genes.
#' @param sp.samples transcriptome indices obtained from the Structure Preserving Sampling step.
#' @param clust.list list object as returned by \code{cluster.cells} module.
#' @return numeric two-column matrix containing the 2D embedding of all data points.
#' @export
compute_2d_embedding<-function(data, sp.samples, clust.list){
  # set.seed(0)
  # data = as.matrix(whole[sample_ids[subsamples_louvain[-ss_clusters$outliers]],top_pc_genes])
  # print("Compute Tsne Projection using PCA top genes...")
  ts<-Rtsne::Rtsne(as.matrix(data[sp.samples,]),
                   perplexity = 15,
                   dims = 2,
                   check_duplicates = FALSE)

  # print("Projecting TSNE co-ordinates...")
  tt<-ts$Y

  INDEX = clust.list$nn.ids
  clust_col = clust.list$cluster.ident
  clust_col[clust_col==0]<-NA

  sample.labels = clust_col[sp.samples]
  PROJ = matrix(NA, nrow=nrow(INDEX), ncol=2)

  nnn = ncol(INDEX)

  proj.idx = 1:nrow(INDEX)

  for( i in proj.idx)
  {
    if(is.na(clust_col[i])){next}
    maj_id = clust_col[i]
    maj_col = which(sample.labels[INDEX[i,]] == maj_id)
    # if(length(maj_col)/nnn > 0.5){
      # PROJ[i,] = apply(tt[INDEX[i,maj_col],],2,function(x) mean(x, na.rm = T))
      PROJ[i,1] = mean(tt[INDEX[i,maj_col],1], na.rm = TRUE)
      PROJ[i,2] = mean(tt[INDEX[i,maj_col],2], na.rm = TRUE)
    # }
  }


  # PROJ[sp.samples,] = tt

  return(PROJ)

}



# ---------------------------------------------
# Read  data
# ---------------------------------------------
#' Read data
#' @description Read 10X or csv data,
#' @details Loads expression data from 10X format directory or csv file as specified. The rows in the matrix are considerd genes; the columns represent samples
#' @param data.path When format="10X", specify direcory location contining 10X files: \cr
#' \enumerate{
#' \item matrix.mtx, sparse matrix format file
#' \item genes.tsv, gene names tab-separated file
#' \item barcodes.tsv, barcodes tab-separated file
#' }; with format="txt", specify the location of the csv matrix file.
#' @param format One of c("10X", "txt"). Defaults to "10X".
#' @param header TRUE/FALSE. Accept the header row to contain gene names within the csv matrix, defaults to TRUE. When FALSE, columns are labelled automatically.
#' @param sep mention the character separator, defaults to ','.
#' @return list containing the 3 input files.
#' @export
read_data<-function(data.path, format="10X", header=TRUE, sep=',', quote="\""){




  data = switch(format,
                "10X" = read10X(data.path),
                "txt" = input_csv(data.path, header, sep, quote),
                stop('Invalid "format" specified, should be one of c("10X","txt")'))



  return(data)
}



input_csv<-function(x,  header, sep, quote){


  data = read.csv(x,
                  header = header,
                  sep = sep,
                  quote = quote)


  if(!(has_rownames(data))) {
    gene_symbols=paste("gene",formatC(1:nrow(data), width=floor(log10(nrow(data)))+1, flag="0"), sep="_")
  } else gene_symbols = rownames(data)


  if(header==T && class(data[,1])=="factor"){
    mat = data[,-1]
    gene_symbols = data[,1]
  }

  barcodes = colnames(mat)
  mat = as.matrix(mat)

  if(header==F) barcodes=paste("cell",formatC(1:ncol(mat), width=floor(log10(ncol(mat)))+1, flag="0"), sep="_")

  mat = t(mat)

  rownames(mat)<- barcodes
  colnames(mat)<-gene_symbols


  return(list("mat" = mat, "barcodes" = barcodes , "gene_symbols" = gene_symbols))
}


## reused Seurat code
read10X<-function(data.dir=NULL){
  full_data <- list()
  for(i in seq_along(data.dir)){
    run <- data.dir[i]
    if (!dir.exists(run)){
      stop("Directory provided does not exist")
    }

    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }

    barcode.loc <- paste(run, "barcodes.tsv", sep ="")
    gene.loc <- paste(run, "genes.tsv", sep ="")
    matrix.loc <- paste(run, "matrix.mtx", sep ="")

    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (!file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }

    data <- Matrix::t(Matrix::readMM(matrix.loc))
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if(all(grepl("\\-1$", cell.names)) == TRUE) {
      cell.names <- as.vector(as.character(sapply(cell.names,
                                                  extract_field,
                                                  1,
                                                  delim = "-")))
    }
    colnames(data) <- make.unique(as.character(sapply(gene.names,
                                                      extract_field,
                                                      2, delim = "\\t")))

    if(is.null(names(data.dir))){
      if(i < 2){
        rownames(data) <- cell.names
      }
      else {
        rownames(data) <- paste0(i, "_", cell.names, sep = "")
      }
    } else {
      rownames(data) <- paste0(names(data.dir)[i],"_",cell.names)
    }

  }

  full_data <- list("mat" = data,
                    "barcodes" = rownames(data) ,
                    "gene_symbols" = colnames(data))
  return(full_data)

}

extract_field=function(string,field=1,delim="_") {
  fields=as.numeric(unlist(strsplit(as.character(field),",")))
  if (length(fields)==1)  return(strsplit(string,delim)[[1]][field])
  return(paste(strsplit(string,delim)[[1]][fields],collapse = delim))
}
