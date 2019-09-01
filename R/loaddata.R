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
#' @importFrom SingleCellExperiment SingleCellExperiment rowData colData
#' @export
readfiles<-function(paths=NULL, format="10x",...){
  switch(format,
         "10x"={
           data.10x <- read10X(data.dir = paths)
           object <- SingleCellExperiment::SingleCellExperiment(list(counts = data.10x))
           SummarizedExperiment::rowData(object)$Symbol = rownames(data.10x)
           SummarizedExperiment::colData(object)$Barcode = colnames(data.10x)
         },
         "txt"={
           data<-.inputcsv(paths, header=TRUE, sep=',', quote="\"", ...)
           object = SingleCellExperiment::SingleCellExperiment(assays = list(counts = data))
           if(is.null(rowData(object)$Symbol))
             SummarizedExperiment::rowData(object)$Symbol<-rownames(object)
           SummarizedExperiment:: rowData(object)$Symbol<-make.unique(rowData(object)$Symbol)
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

## reused Seurat code
read10X<-function(data.dir=NULL,  gene.column = 2, unique.features = TRUE){
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing")
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- Matrix::readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      sep="",
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, ])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }

}



extract_field=function(string,field=1,delim="_") {
  fields=as.numeric(unlist(strsplit(as.character(field),",")))
  if (length(fields)==1)  return(strsplit(string,delim)[[1]][field])
  return(paste(strsplit(string,delim)[[1]][fields],collapse = delim))
}
