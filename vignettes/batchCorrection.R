## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval=FALSE,
  comment = "#>"
)

## ------------------------------------------------------------------------
#  
#  library(dropClust)
#  load(url("https://raw.githubusercontent.com/LuyiTian/CellBench_data/master/data/sincell_with_class.RData"))
#  
#  set.seed(0)
#  objects = list()
#  
#  objects[[1]] = sce_sc_10x_qc
#  
#  objects[[2]] = sce_sc_CELseq2_qc
#  
#  objects[[3]] = sce_sc_Dropseq_qc
#  

## ------------------------------------------------------------------------
#  
#  all.objects = objects
#  merged_data<-Merge(all.objects)
#  
#  annotations = merged_data$cell_line
#  batch.id = merged_data$Batch

## ------------------------------------------------------------------------
#  
#  dc.corr <- Correction(merged_data, close_th = 0.1, cells_th = 0.1,
#                         components = 10, n_neighbors = 30, init = "spca")

## ------------------------------------------------------------------------
#  cc = Cluster(dc.corr,method = "kmeans",centers = 3)

## ------------------------------------------------------------------------
#  de<-FindMarkers(cc,q_th = 0.001, lfc_th = 1.2)

