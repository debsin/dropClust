# ---------------------------------------------
#' Read data
#' @description Read 10X or csv data,
#' @details Loads expression data from 10X format directory or csv file as specified. The rows in the matrix are considered genes; the columns represent samples
#' @param paths a path or a list of paths to the csv file/s, or the directory location/s containing the 10X files: \cr
#' \enumerate{
#' \item matrix.mtx, sparse matrix format file
#' \item genes.tsv, gene names tab-separated file
#' \item barcodes.tsv, barcodes tab-separated file
#' }; with format="txt", specify the location of the csv matrix file.
#' @param format One of c("10X", "txt"). Defaults to "10X".
#' @param ...  pass specific arguments to the \code{read.csv()} function when \code{format = "txt"}.
#' @return SingleCellExperiment object.
#' @importFrom utils read.csv
#' @importFrom DropletUtils read10xCounts
#' @importFrom SingleCellExperiment SingleCellExperiment rowData
#' @export
readfiles<-function(paths=NULL, format="10x",...){
  switch(format,
         "10x"={
           object <- DropletUtils::read10xCounts(samples = paths, col.names = T,type = 'auto')
           if(is.null(rowData(object)$Symbol))
             SummarizedExperiment::rowData(object)$Symbol<-rownames(object)
           SummarizedExperiment:: rowData(object)$Symbol<-make.unique(rowData(object)$Symbol)
         },
         "txt"={
           data<-.inputcsv(paths, header=TRUE, sep=',', quote="\"", ...)
           object = SingleCellExperiment::SingleCellExperiment(assays = list(counts = data))
         },
         {
           stop("Unknown format specified.")
         }
  )

  cat(dim(object)[2],"Cells and",dim(object)[1], "genes present.\n")
  object@metadata[["dropClust"]] = ""

  return(object)
}

#' @importFrom tibble has_rownames
.inputcsv<-function(matrices,  header, sep, quote){

  nsets <- length(matrices)
  full_data <- vector("list", nsets)
  gene_info_list <- vector("list", nsets)
  cell_info_list <- vector("list", nsets)

  for (i in seq_len(nsets)) {
    x <- matrices[i]

    # cur.type <- .type_chooser(run, format)
    # if (cur.type=="sparse") {
    #   info <- .read_from_sparse(run, version=version)
    # } else {
    #   info <- .read_from_hdf5(run, genome=genome, version=version)
    # }

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

    colnames(mat) <- barcodes
    rownames(mat) <- gene_symbols

    full_data[[i]] <- mat
    gene_info_list[[i]] <- gene_symbols
    cell.names <- barcodes
    cell_info_list[[i]] <- S4Vectors::DataFrame(Sample = rep(x, length(cell.names)), Barcode = cell.names, row.names=NULL)
  }

  # Checking gene uniqueness.
  if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
    stop("gene information differs between sets")
  }
  gene_info <- gene_info_list[[1]]
  rownames(gene_info) <- gene_info

  # Forming the full data matrix.
  full_data <- do.call(cbind, full_data)

  # Adding the cell data.
  cell_info <- do.call(rbind, cell_info_list)

    if (nsets == 1L) {
      cnames <- cell_info$Barcode
    } else {
      sid <- rep(seq_along(cell_info_list), vapply(cell_info_list, nrow, 1L))
      cnames <- paste0(sid, "_", cell_info$Barcode)
    }
    colnames(full_data) <- cnames


  SingleCellExperiment(list(counts = full_data), rowData = gene_info, colData = cell_info)
}



extract_field=function(string,field=1,delim="_") {
  fields=as.numeric(unlist(strsplit(as.character(field),",")))
  if (length(fields)==1)  return(strsplit(string,delim)[[1]][field])
  return(paste(strsplit(string,delim)[[1]][fields],collapse = delim))
}
