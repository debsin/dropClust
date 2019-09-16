library(DropletUtils)
library(SingleCellExperiment)
library("magrittr")


# set.seed(0)
sce <-readfiles(path = "C:/Projects/data/pbmc68K/hg19/")
annotations = read.csv("C:/Users/debajyoti/Desktop/68k_pbmc_barcodes_annotation.tsv",sep="\t",header = T)
# sce.copy<-sce

s  = Sys.time()

sce.1<-FilterCells(sce)
sce.1<-FilterGenes(sce.1)
sce.2<-CountNormalize(sce.1)

sce.2<-RankGenes(sce.2,ngenes_keep = 1000)

sce.3<-Sampling(sce.2)

sce.3<-RankPCAGenes(sce.3)

sce.3<-Cluster(sce.3, method = "default", conf = 0.8, use.previous = F, minClusterSize = 50, deepSplit = 3, center = 12)
sce.3@metadata$dropClust

sc_metric(sce.3$ClusterIDs, annotations$celltype[match(sce.3$Barcode, annotations$barcodes)])

sce.3<-PlotEmbedding(sce.3, embedding = "umap", spread = 10, min_dist = 0.1)
sce.3@metadata$dropClust

Sys.time() - s


plot_data = data.frame("Y1" = reducedDim(sce.3,"umap")[,1], Y2 = reducedDim(sce.3, "umap")[,2], color = annotations$celltype[match(sce.3$Barcode, annotations$barcodes)])
dropClust::all_plot(plot_data,title = "Annotations")

# de_res = FindMarkers(sce.3,lfc_th = 1,q_th = 0.001)




