#' Merge datasets into one by determining a common genes
#' @description Merge the count data from different batches into one object
#' @details Concatenate the expression counts of all cells from different batches into one expression count object.
#' The merging is done on the set of union of DE genes obtained from the clustering of each batch.
#' @param objects A list of SingleCellExperiment objects, each representing a  SingleCellExperiment object from a single batch.
#' @param use.de.genes logical to specify if DE genes shoud be computed,
#' when \code{FALSE} the batches are merged on the set of common genes across batches.
#' @return  A SingleCellExperiment object with an additional column named \code{Batch} in \code{colData} column.
#' The column stores the origin of the batch for each cell. The different batches are labeled in the order
#' as appearing in the input list.
#' @importFrom SingleCellExperiment rowData colData
#' @export
#' @examples
#' library(SingleCellExperiment)
#' ncells <- 500
#' ngenes <- 2000
#' counts.1 <- matrix(rpois(ncells*ngenes, lambda = 10), ncol=ncells, nrow=ngenes, byrow=TRUE)
#' rownames(counts.1) <- paste0("Gene", seq_len(ngenes))
#' colnames(counts.1) <- paste0("Cell", seq_len(ncells))
#' sce.1 <- SingleCellExperiment(assays = list(counts = counts.1))
#'
#' counts.2 <- matrix(rpois(ncells*ngenes, lambda = 15), ncol=ncells, nrow=ngenes, byrow=TRUE)
#' rownames(counts.2) <- paste0("Gene", seq_len(ngenes))
#' colnames(counts.2) <- paste0("Cell", seq_len(ncells))
#' sce.2 <- SingleCellExperiment(assays = list(counts = counts.2))
#'
#' mix_sce <- Merge(list(sce.1, sce.2), use.de.genes =FALSE)
Merge<-function(objects, use.de.genes = TRUE){

  if(length(objects)==1) stop("Multiple objects must be present to integrate.")


  if(is.null(rownames(objects[[1]])))
    stop("Objects must contain rownames.")

  if(is.null(rowData(objects[[1]])$Symbol))
    SummarizedExperiment::rowData(objects[[1]])$Symbol <- rownames(objects[[1]])

  common_genes_raw = rownames(objects[[1]])
  common_colData = names(colData(objects[[1]]))

  all_de_genes = list()
  merge_genes = list()

  for (i in 1:length(objects)){
    cat("\nProcesing Batch",i,"...\n")

    if(is.null(rownames(objects[[i]])))
       stop("Objects must contain rownames.")

    if(is.null(rowData(objects[[i]])$Symbol))
      SummarizedExperiment::rowData(objects[[i]])$Symbol <- rownames(objects[[i]])

    # rownames(objects[[i]]) = rowData(objects[[i]])$Symbol

    common_genes_raw = intersect(common_genes_raw, rownames(objects[[i]]))
    temp = objects[[i]]

    temp<- FilterCells(temp)
    temp<- FilterGenes(temp, min_count = 1, min_cell = 1)
    temp<- CountNormalize(temp)
    merge_genes[[i]] = rownames(temp)


    common_colData = intersect(common_colData, names(colData(temp)))

    temp<- RankGenes(temp,ngenes_keep = 1000)
    merge_genes[[i]] = rownames(temp)[rowData(temp)$HVG]


    if(use.de.genes){

      temp<- Sampling(temp)
      temp<- RankPCAGenes(temp)
      temp<- Cluster(temp, method = "default", conf = 0.1, use.previous = F)



      GRP = sort(unique(temp$ClusterIDs))
      de_res = FindMarkers(temp,selected_clusters = GRP, lfc_th = 1, q_th = 0.001, nDE = 20)
      batch.gene.symbol = as.vector(de_res$genes.df) #returns gene symbols
      all_de_genes[[i]] = rownames(temp)[match(batch.gene.symbol, rowData(temp)$Symbol,nomatch=NA)]
      all_de_genes[[i]] = all_de_genes[[i]][!is.na(all_de_genes[[i]])]
    }
  }

  hvg_genes  = base::Reduce(union, merge_genes)

  common_genes = intersect(common_genes_raw, hvg_genes)

  mix.data.counts  =  counts(objects[[1]])[common_genes, ]
  # mix.data.ncounts = normcounts(objects[[1]])[common_genes, ]
  col.data = colData(objects[[1]])[,common_colData, drop=F]
  batch = rep(1, ncol(mix.data.counts))

  row.data = SummarizedExperiment::rowData(objects[[1]])[common_genes, , drop=F]

  for(i in 2:length(objects)){
    temp  <-      counts(objects[[i]])[common_genes, ]
    mix.data.counts  <- SingleCellExperiment::cbind(mix.data.counts, temp)

    # temp2 <-  normcounts(objects[[i]])[common_genes, ]
    # mix.data.ncounts <- SingleCellExperiment::cbind(mix.data.ncounts, temp2)

    col.data = SingleCellExperiment::rbind(col.data, colData(objects[[i]])[,common_colData, drop=F])
    batch  = c(batch, rep(i, ncol(temp)))
  }


  mix.data =  SingleCellExperiment(assays = list(counts = mix.data.counts),
                                   rowData = row.data, colData = col.data)

  SummarizedExperiment::colData(mix.data)$Batch = batch

  if(use.de.genes){
    common_genes_de = intersect(unique(as.vector(unlist(all_de_genes))), unique(unlist(common_genes)))
    SummarizedExperiment::rowData(mix.data)$CommonDEGenes = rep(FALSE, nrow(mix.data))
    SummarizedExperiment::rowData(mix.data)$CommonDEGenes[which(rownames(mix.data) %in% common_genes_de)]<- TRUE
  }

  mix.data@metadata[["dropClust"]] = c(mix.data@metadata[["dropClust"]], "Merge")


  return(mix.data)
}
