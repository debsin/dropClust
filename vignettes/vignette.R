## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval=TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(dropClust)
set.seed(0)

## ------------------------------------------------------------------------
# Load Data
# path contains decompressed files 
sce <-readfiles(path = "C:/Projects/dropClust/data/pbmc3k/hg19/")


## ------------------------------------------------------------------------
# Filter poor quality cells.  A threshold th corresponds to the total count of a cell.
sce<-FilterCells(sce)
sce<-FilterGenes(sce)

## ------------------------------------------------------------------------
sce<-CountNormalize(sce)

## ------------------------------------------------------------------------
sce<-RankGenes(sce, ngenes_keep = 1000)

## ------------------------------------------------------------------------
sce<-Sampling(sce)

## ------------------------------------------------------------------------
# Find PCA top 200 genes. This may take some time.
sce<-RankPCAGenes(sce)

## ------------------------------------------------------------------------
sce<-Cluster(sce, method = "default", conf = 0.5)

## ------------------------------------------------------------------------
sce<-PlotEmbedding(sce, embedding = "umap", spread = 10, min_dist = 0.1)
ScatterPlot(sce,title = "Clusters")

## ------------------------------------------------------------------------

DE_genes_all = FindMarkers(sce, selected_clusters=NA, lfc_th = 1, q_th =0.001, nDE=30)

write.csv(DE_genes_all$genes, 
          file = file.path(tempdir(),"ct_genes.csv"),
          quote = FALSE)

## ------------------------------------------------------------------------

marker_genes = c("S100A8", "GNLY", "PF4")

p<-PlotMarkers(sce, marker_genes)



## ------------------------------------------------------------------------
# Draw heatmap
#############################
p<-PlotHeatmap(sce, DE_res = DE_genes_all$DE_res,nDE = 10)

print(p)



