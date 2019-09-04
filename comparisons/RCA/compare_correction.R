library(DropletUtils)
library(SingleCellExperiment)
library("magrittr")
## RCA dataset

rca.data <- read.csv('C:/Projects/data/GSE81861_Cell_Line_COUNT.csv',header = T,row.names = 1,stringsAsFactors = F)
gene.symbols = make.unique(unlist(lapply(rownames(rca.data), function(x) toupper(unlist(strsplit(x, split = "_"))[2]))))

labels <- colnames(rca.data)
## Stripping cell id and batch no from each label
annotations <- unlist(lapply(labels, function(x) toupper(unlist(strsplit(x, split = "[_]"))[3])))
cellnames = unlist(lapply(labels, function(x) toupper(unlist(strsplit(x, split = "[_]"))[1])))

counts = Matrix:: Matrix(as.matrix(rca.data), sparse = T)

batch = rep(1, length(labels))
batch[grep('B2',labels)]<-2

source("~/GitHub/dropClust/comparisons/RCA/plots.R")

sce1 <- SingleCellExperiment(assays = list(counts = counts[,which(batch ==1)]))
rowData(sce1)$Symbol = gene.symbols
colData(sce1)$cell_line = annotations[which(batch ==1)]
colData(sce1)$Batch = batch[which(batch ==1)]

sce2 <- SingleCellExperiment(assays = list(counts = counts[,which(batch ==2)]))
rowData(sce2)$Symbol = gene.symbols
colData(sce2)$cell_line = annotations[which(batch ==2)]
colData(sce2)$Batch = batch[which(batch ==2)]

library(dropClust)
set.seed(0)
objects  = list(sce1,sce2)

all.objects = objects
merged_data<-Merge(all.objects)


# dropClust Corrected
set.seed(1)
corrected_data <- Correction(merged_data, close_th = 0.07, cells_th = 0.1,
                       components = 10, n_neighbors = 10, min_dist = 0.5)

PROJ_c_dc = reducedDim(corrected_data, "CComponents")
plot_proj_df_true_dc<-data.frame(Y1 = PROJ_c_dc[,1],Y2 = PROJ_c_dc[,2],
                                 color = as.factor(corrected_data$cell_line),
                                 batch = as.integer(corrected_data$Batch))
batch_plot(plot_proj_df_true_dc,filename = NA, title = "Corrected DC ",  type=NULL)


# Uncorrected
mixed.UC <- FilterGenes(merged_data)
mixed.UC <- CountNormalize(mixed.UC)
mixed.UC <- RankGenes(mixed.UC, ngenes_keep = 1000)
uncorrected_mat <- t(normcounts(mixed.UC)[rowData(mixed.UC)$HVG,])
set.seed(1)
PROJ_uc = uwot::umap(as.matrix(uncorrected_mat),metric = 'cosine',n_neighbors =20)
plot_proj_df_true_uc<-data.frame(Y1 = PROJ_uc[,1],Y2 = PROJ_uc[,2],color = as.factor(mixed.UC$cell_line), batch = as.factor(mixed.UC$Batch))
batch_plot(plot_proj_df_true_uc,filename = NA, title = "Uncorrected",  type=NULL)



# MnnCorrect
source("~/GitHub/dropClust/comparisons/RCA/mnnCorrect.R")
dim(mnn_corr)
PROJ_mnn = uwot::umap(mnn_corr, n_neighbors = 20 )
plot_proj_df_true_mnn<-data.frame(Y1 = PROJ_mnn[,1],Y2 = PROJ_mnn[,2],color = as.factor(annotations), batch = as.factor(batch))
batch_plot(plot_proj_df_true_mnn,filename = NA, title = "Corrected Mnn ",  type=NULL)

# Seurat
source("~/GitHub/dropClust/comparisons/RCA/Seurat.R")

dim(seurat_corr)

# umap_proj_su = uwot::umap(seurat_corr, n_neighbors =20 )
PJ_su<-seurat_corr
plot_proj_df_true_surat<-data.frame(Y1 = PJ_su[,1],Y2 = PJ_su[,2],color = as.factor(annotations), batch = as.factor(batch))
batch_plot(plot_proj_df_true_surat,filename = NA, title = "Corrected Seurat",  type=NULL)


# Scanorama
source("~/GitHub/dropClust/comparisons/RCA/scanorama.R")
dim(scrama_corr)
umap_proj_scanorama = uwot::umap(scrama_corr,n_neighbors = 20 , metric = 'cosine')
PJ_scanorama<-umap_proj_scanorama
plot_proj_df_true_scanorama<-data.frame(Y1 = PJ_scanorama[,1],Y2 = PJ_scanorama[,2],
                                        color = as.factor(annotations), batch = as.factor(batch))
batch_plot(plot_proj_df_true_scanorama,filename = NA, title = "Corrected Scanorama",  type=NULL)





interested_types = c("All")

p<-list()
for(type in interested_types){
  if(type=="All") type = NULL
  p[[paste(type,"UC",collapse=" ")]]<- batch_plot(plot_proj_df_true_uc,filename = NA, title = "Uncorrected",  type=type)
  p[[paste(type,"DC")]]<- batch_plot(plot_proj_df_true_dc,filename = NA, title = "dropClust",  type=type)
  p[[paste(type,"Seurat")]]<- batch_plot(plot_proj_df_true_surat,filename = NA, title = "Seurat",  type=type)
  p[[paste(type,"Scanorama")]]<- batch_plot(plot_proj_df_true_scanorama,filename = NA, title = "Sconorama",  type=type)
  p[[paste(type,"Mnn")]]<- batch_plot(plot_proj_df_true_mnn,filename = NA, title = "MNNCorrect",  type=type)

}

p[[paste(type,"legend")]] <- cowplot::get_legend(legend_plot(plot_proj_df_true_uc,filename = NA, title = "Uncorrected",  type=NULL))

hlay <- rbind(c(1,1,2,2,3,3),
              c(1,1,2,2,3,3),
              c(1,1,2,2,3,3),
              c(NA,4,4,5,5,NA),
              c(NA,4,4,5,5,NA),
              c(NA,4,4,5,5,NA),
              c(6,6,6,6,6,6))

pdf("~/GitHub/dropClust/comparisons/RCA/all_types_complete_rca_mixed.pdf",width = 12, height = 8)
gridExtra::grid.arrange(grobs = p, layout_matrix = hlay)
dev.off()



unique(annotations)
interested_types = c("GM12878", "H1", "All")

p<-list()
for(type in interested_types){
  if(type=="All") type = NULL
  p[[paste(type,"UC",collapse=" ")]]<- batch_plot(plot_proj_df_true_uc,filename = NA, title = "Uncorrected",  type=type)
  p[[paste(type,"DC")]]<- batch_plot(plot_proj_df_true_dc,filename = NA, title = "dropClust",  type=type)
  p[[paste(type,"Seurat")]]<- batch_plot(plot_proj_df_true_surat,filename = NA, title = "Seurat",  type=type)
  p[[paste(type,"Scanorama")]]<- batch_plot(plot_proj_df_true_scanorama,filename = NA, title = "Sconorama",  type=type)
  p[[paste(type,"Mnn")]]<- batch_plot(plot_proj_df_true_mnn,filename = NA, title = "MNNCorrect",  type=type)

}

library(gridExtra)
# pdf("compararisons/RCA/all_types_complete.pdf",width = 14, height = 10)
# do.call("grid.arrange", c(p, ncol = 5))
# dev.off()

