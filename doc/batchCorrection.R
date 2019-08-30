## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval=TRUE,
  comment = "#>"
)

## ----warning=FALSE,message=FALSE-----------------------------------------

library(dropClust)
load(url("https://raw.githubusercontent.com/LuyiTian/CellBench_data/master/data/sincell_with_class.RData"))

objects = list()

objects[[1]] = sce_sc_10x_qc

objects[[2]] = sce_sc_CELseq2_qc

objects[[3]] = sce_sc_Dropseq_qc


## ---- echo=FALSE---------------------------------------------------------

all.objects = objects
merged_data<-Merge(all.objects)


## ------------------------------------------------------------------------
set.seed(1)
dc.corr <-  Correction(merged_data,  method="default", close_th = 0.1, cells_th = 0.1,
                       components = 10, n_neighbors = 20,  min_dist = 0.5)

## ------------------------------------------------------------------------
dc.corr = Cluster(dc.corr,method = "kmeans",centers = 3)

## ---- fig.asp=1, fig.cap="Batch corrected dropClust based Clustering."----
ScatterPlot(dc.corr, title = "Clusters")


## ---- message=FALSE, warning=FALSE---------------------------------------
merged_data.fastmnn<-Merge(all.objects,use.de.genes = FALSE)

set.seed(1)
mnn.corr <-  Correction(merged_data.fastmnn,  method="fastmnn", d = 10)

mnn.corr = Cluster(mnn.corr,method = "kmeans",centers = 3)

ScatterPlot(mnn.corr, title = "Clusters")


## ------------------------------------------------------------------------
de<-FindMarkers(dc.corr,q_th = 0.001, lfc_th = 1.2,nDE = 10)
de$genes.df

