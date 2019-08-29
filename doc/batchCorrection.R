## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval=FALSE,
  comment = "#>"
)

## ------------------------------------------------------------------------
#  
#  rca.data <- read.csv('C:/Projects/data/GSE81861_Cell_Line_COUNT.csv',header = T,row.names = 1,stringsAsFactors = F)
#  gene.symbols = make.unique(unlist(lapply(rownames(rca.data), function(x) toupper(unlist(strsplit(x, split = "_"))[2]))))
#  
#  labels <- colnames(rca.data)
#  ## Stripping cell id and batch no from each label
#  annotations <- unlist(lapply(labels, function(x) toupper(unlist(strsplit(x, split = "[_]"))[3])))
#  cellnames = unlist(lapply(labels, function(x) toupper(unlist(strsplit(x, split = "[_]"))[1])))
#  
#  counts = Matrix:: Matrix(as.matrix(rca.data), sparse = T)
#  
#  batch = rep(1, length(labels))
#  batch[grep('B2',labels)]<-2
#  
#  
#  sce1 <- SingleCellExperiment(assays = list(counts = counts[,which(batch ==1)]))
#  rowData(sce1)$Symbol = gene.symbols
#  colData(sce1)$cell_line = annotations[which(batch ==1)]
#  colData(sce1)$Batch = batch[which(batch ==1)]
#  
#  sce2 <- SingleCellExperiment(assays = list(counts = counts[,which(batch ==2)]))
#  rowData(sce2)$Symbol = gene.symbols
#  colData(sce2)$cell_line = annotations[which(batch ==2)]
#  colData(sce2)$Batch = batch[which(batch ==2)]
#  
#  
#  set.seed(0)
#  objects  = list(sce1,sce2)
#  
#  all.objects = objects
#  merged_data<-Merge(all.objects)
#  
#  
#  
#  annotations = mix.data$cell_line
#  batch.id = mix.data$Batch
#  
#  
#  dc.corr <- Correction(merged_data, close_th = 0.1, cells_th = 0.1,
#                         components = 10, n_neighbors = 30, init = "spca")
#  
#  
#  
#  cc = Cluster(dc.corr,method = "kmeans",centers = 3)
#  
#  
#  de<-FindMarkers(cc,q_th = 0.001, lfc_th = 1.2)

