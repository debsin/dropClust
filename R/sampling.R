#' Sampling Primary Clusters
#' @description Performs sampling from the primary clusters in an inverse exponential order of cluster size.
#' @details Sampling in inverse proportion of cluster size following a exponential decay equation.
#' To ensure selection of sufficient representative transcriptomes from small clusters,
#' an exponential decay function  is used to determine the proportion of transcriptomes to be sampled from each
#' cluster. For $i^{th}$ cluster, the proportion of expression profiles $p_i$ was obtained as follows.\cr
#' \ifelse{html}{\out{p<sub>i</sub> = p<sub>l</sub> <font face="symbol">-</font> e<sup><font face="symbol">-</font>(S<sub>i</sub>)/(K)</sup>}}{\deqn{Lp_{i} = p_{l} - e^{-\frac{S_i}{K}} (p_{l} - p_{u})}{ASCII}}
#' where \eqn{S_i} is the size of cluster \eqn{i}, \eqn{K} is a scaling factor, \eqn{p_i} is the proportion of cells to be sampled from the $i^{th}$ Louvain cluster. $p_l$ and $p_u$ are lower and upper bounds of the proportion value respectively.
#' @references {
#' \insertRef{sengupta2013reformulated}{dropClust}
#' }
#' @param object A SingleCellExperiment object containing normalized expression values in \code{"normcounts"}.
#' @param nsamples integer, total number of samples to return post sampling; ignored when \code{optm_parameters = FALSE}.
#' @param method character, one of c("sps","random"). Structure Preserving Sampling (sps) selects proportional number of members from each cluster obtained from partitioning an approximate nearest neighbour graph.
#' @param optm_parameters logical, when TRUE the parameters (\code{pinit, pfin, K}) are optimized such that exactly \code{nsamples} are returned. Optimization is performed using simulated annealing
#' @param pinit numeric [0,0.5], minimum probability of that sampling occurs from a cluster, ignored when \code{optm_parameters = TRUE}.
#' @param pfin numeric [0.5,1], maximum probability of that sampling occurs from a cluster, ignored when \code{optm_parameters = TRUE}.
#' @param K numeric, scaling factor analogous to Boltzmann constant, ignored when \code{optm_parameters = TRUE}.
#' @return A SingleCellExperiment object with an additional column named \code{Sampling} in \code{colData} column.
#' The column stores a a logical value against each cell  to indicate if it has been sampled.
#' @importFrom SingleCellExperiment reducedDimNames normcounts colData
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
#' sce <- Sampling(sce)
Sampling<-function(object,
                   nsamples=500,
                   method = "sps",
                   optm_parameters=FALSE,
                   pinit=0.195,
                   pfin = 0.9,
                   K=500){

  if(nsamples>=ncol(object)) {
    SummarizedExperiment::colData(object)$Sampling=  rep(TRUE, ncol(object))
    return(object)
  }
  if(!any(method %in% c("random","sps")))
    stop("Method not found.")

  no_samples = ncol(object)
  init = ifelse(no_samples<20000, no_samples,  min(20000,round(no_samples/3)))
  sample_ids = sample(1:no_samples, init)

  if(method=="sps"){

    if(!any(reducedDimNames(object)=="CComponents"))
      data = Log2Normalize(normcounts(object)[SingleCellExperiment::rowData(object)$HVG, sample_ids],return.sparse = FALSE)
    else
      data = as.matrix(normcounts(object)[, sample_ids])
    # data = as.matrix(assay(object,i = "logcounts")[SingleCellExperiment::rowData(object)$HVG,idx])

    partition<-.annPartition(data)

    if(optm_parameters==TRUE){
      param = .optimized_param(partition, nsamples)
      pinit = param[1]
      pfin = param[2]
      K = param[3]
      cat("Optimized parameters:\n", param,"\n")
    }

    oldseed = .Random.seed
    cluster_freq = table(partition[,2])
    prop = round((pinit -
                    exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
    t(rbind(cluster_freq,prop))
    prop = reshape2::melt(prop)$value
    subsamples<-c()
    for(k in 1:length(prop)){
      subsamples = c(subsamples,
                     sample(partition[which(partition[,2]==(k)),1],
                            size = prop[k],replace = FALSE))
    }
    .Random.seed = oldseed
    SummarizedExperiment::colData(object)$Sampling = rep(FALSE, ncol(object))
    SummarizedExperiment::colData(object)$Sampling[sample_ids[subsamples]] =  TRUE
    cat(length(subsamples), "samples extracted.\n")
  }

  else if(method=="random"){
    oldseed = .Random.seed
    subsamples = sample(sample_ids, nsamples)
    .Random.seed = oldseed
    SummarizedExperiment::colData(object)$Sampling = rep(FALSE, ncol(object))
    SummarizedExperiment::colData(object)$Sampling[subsamples] =  TRUE
  }
  else{
    cat("Invalid sampling. Fallback to all samples.")
  }

  object@metadata[["dropClust"]] = c(unlist(object@metadata[["dropClust"]]),"Sampling")


  object
}



# ---------------------------------------------
# Optimize parameter 'pinit' for
# exponential decay sampling
# ---------------------------------------------
.optimized_param<-function(partition, nsamples = 500){
  if(nsamples > nrow(partition)) nsamples = nrow(partition)
  oldseed = .Random.seed
  # set.seed(0)
  global.min <- 0
  tol <- 1e-3

  max.time= 20
  lower <- c(0.05, 0.9,500)
  upper <- c(0.1,0.95,4000)

  params<-NULL

  pin_find<-function(params, MAX_c = nsamples){
    output <- as.numeric(partition)
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


#' ANN Graph Partition
#' @param data numeric matrix with \code{n} samples in rows
#' @return two column integer matrix, first column represents sample id, second column contain cluster membership id
.annPartition<-function(data){
  # set.seed(0)
  f = dim(data)[1]
  cat("Building graph with", ncol(data), "nodes...")
  t = methods::new(RcppAnnoy::AnnoyAngular,f)  # Length of item vector that will be indexed
  for(i in seq(ncol(data))){
    v = data[,i]
    t$addItem(i, v)
  }
  t$build(30)# 30 trees

  indices = lapply(seq(ncol(data)),function(x) t$getNNsByItem(x, 6))

  G=igraph::simplify(igraph::graph_from_adj_list(indices,
                                                 mode="all",
                                                 duplicate = FALSE))
  cat("Louvain Partition...")

  partition = igraph::cluster_louvain(G)
  cat("Done.\n")

  dataMatrix =as.data.frame(cbind(seq(ncol(data)),partition$membership))
  return(dataMatrix)
}
