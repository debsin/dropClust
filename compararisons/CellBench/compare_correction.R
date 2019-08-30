library(DropletUtils)
library(SingleCellExperiment)
library("magrittr")
## CellBench dataset

load("C:/Projects/data/Cellbench/sincell_with_class.RData")
load("~/GitHub/dropclust_benchmarks/analysis/bc/Cellbench/out_gene_ids.Rda")

source("compararisons/Cellbench/plots.R")


set.seed(0)
objects = list()
symbols = list("CelSeq" = out, "10x" = out2, "DropSeq" = out3)

na = match(symbols$`10x`$ensembl_gene_id, rownames(sce_sc_10x_qc))
objects[[1]] = sce_sc_10x_qc[na,]
rowData(objects[[1]])$Symbol = make.unique(symbols$`10x`$external_gene_name[match(symbols$`10x`$ensembl_gene_id, rownames(objects[[1]]))])

na = match(symbols$CelSeq$ensembl_gene_id, rownames(sce_sc_CELseq2_qc))
objects[[2]] = sce_sc_CELseq2_qc[na,]
rowData(objects[[2]])$Symbol = make.unique(symbols$CelSeq$external_gene_name[match(symbols$CelSeq$ensembl_gene_id, rownames(objects[[2]]))])

na = match(symbols$DropSeq$ensembl_gene_id, rownames(sce_sc_Dropseq_qc))
objects[[3]] = sce_sc_Dropseq_qc[na,]
rowData(objects[[3]])$Symbol = make.unique(symbols$DropSeq$external_gene_name[match(symbols$DropSeq$ensembl_gene_id, rownames(objects[[3]]))])


annotations = c(objects[[1]]$cell_line, objects[[2]]$cell_line, objects[[3]]$cell_line)
batch.id = c(rep(1, ncol(objects[[1]])), rep(2, ncol(objects[[2]])), rep(3, ncol(objects[[3]])))

all.objects = objects
merged_data<-Merge(all.objects)


# dropClust Corrected
set.seed(1)
corrected_data <- Correction(merged_data, close_th = 0.1, cells_th = 0.1,
                       components = 10, n_neighbors = 20,  min_dist = 0.5)

PROJ_c_dc = reducedDim(corrected_data, "CComponents")
plot_proj_df_true_dc<-data.frame(Y1 = PROJ_c_dc[,1],Y2 = PROJ_c_dc[,2],
                                 color = as.factor(corrected_data$cell_line),
                                 batch = as.integer(corrected_data$Batch))
batch_plot(plot_proj_df_true_dc,filename = NA, title = "Corrected DC ",  type=NULL)


# Uncorrected
merged_data<-Merge(all.objects, use.de.genes = FALSE)
mixed.UC <- CountNormalize(merged_data)
mixed.UC <- RankGenes(mixed.UC, ngenes_keep = 1000)
uncorrected_mat <- t(normcounts(mixed.UC)[rowData(mixed.UC)$HVG,])
PROJ_uc = uwot::umap(as.matrix(uncorrected_mat),n_neighbors =20 )
plot_proj_df_true_uc<-data.frame(Y1 = PROJ_uc[,1],Y2 = PROJ_uc[,2],color = as.factor(mixed.UC$cell_line), batch = as.factor(mixed.UC$Batch))
batch_plot(plot_proj_df_true_uc,filename = NA, title = "Uncorrected",  type=NULL)



# MnnCorrect
source("compararisons/CellBench/mnnCorrect.R")
dim(mnn_corr)
PROJ_mnn = uwot::umap(mnn_corr, n_neighbors = 20 )
plot_proj_df_true_mnn<-data.frame(Y1 = PROJ_mnn[,1],Y2 = PROJ_mnn[,2],color = as.factor(annotations), batch = as.factor(batch.id))
batch_plot(plot_proj_df_true_mnn,filename = NA, title = "Corrected Mnn ",  type=NULL)

# Seurat
source("compararisons/CellBench/Seurat.R")

dim(seurat_corr)

# umap_proj_su = uwot::umap(seurat_corr, n_neighbors =20 )
PJ_su<-seurat_corr
plot_proj_df_true_surat<-data.frame(Y1 = PJ_su[,1],Y2 = PJ_su[,2],color = as.factor(annotations), batch = as.factor(batch.id))
batch_plot(plot_proj_df_true_surat,filename = NA, title = "Corrected Seurat",  type=NULL)


# Scanorama
source("compararisons/CellBench/scanorama.R")
dim(scrama_corr)
umap_proj_scanorama = uwot::umap(scrama_corr,n_neighbors = 20 )
PJ_scanorama<-umap_proj_scanorama
plot_proj_df_true_scanorama<-data.frame(Y1 = PJ_scanorama[,1],Y2 = PJ_scanorama[,2],color = as.factor(annotations), batch = as.factor(batch.id))
batch_plot(plot_proj_df_true_scanorama,filename = NA, title = "Corrected Scanorama",  type=NULL)




unique(annotations)
interested_types = c(unique(annotations),"All")

p<-list()
for(type in interested_types){
  if(type=="All") type = NULL
  p[[paste(type,"UC",collapse=" ")]]<- batch_plot(plot_proj_df_true_uc,filename = NA, title = "Uncorrected",  type=type)
  p[[paste(type,"DC")]]<- batch_plot(plot_proj_df_true_dc,filename = NA, title = "dropClust",  type=type)
  p[[paste(type,"Seurat")]]<- batch_plot(plot_proj_df_true_surat,filename = NA, title = "Seurat",  type=type)
  p[[paste(type,"Scanorama")]]<- batch_plot(plot_proj_df_true_scanorama,filename = NA, title = "Sconorama",  type=type)
  p[[paste(type,"Mnn")]]<- batch_plot(plot_proj_df_true_mnn,filename = NA, title = "mnncorrect",  type=type)

}

library(gridExtra)
pdf("compararisons/CellBench/all_types_complete.pdf",width = 12, height = 11)
do.call("grid.arrange", c(p, ncol = 5))
dev.off()


###
# Alternate Corrections

# batchc = sva::ComBat(as.matrix(t(dataAll)), batch.id, mod = NULL, par.prior = F,  mean.only	=T)
# umap_proj_combat = umap::umap(t(batchc), method = "umap-learn", metric = 'cosine',n_neighbors = 10 )
# PROJ_combat = umap_proj_combat$layout
# plot_proj_df_true_combat<-data.frame(Y1 = PROJ_combat[,1],Y2 = PROJ_combat[,2],color = as.factor(annotations), batch = as.factor(batch.id))
# batch_plot(plot_proj_df_true_combat,filename = NA, title = "Corrected UMAP",  type="H1")

# library("scran")
# mnn_mat  =  scran::fastMNN(log(as.matrix(t(wholeA))+1), log(as.matrix(t(wholeB))+1),k=20)
# umap_proj_c = umap::umap(mnn_mat$corrected, method = "umap-learn", metric = 'cosine',n_neighbors = 10 )
# PROJ_c = umap_proj_c$layout

# # Scatter Plot
# plot_proj_df_true_c<-data.frame(Y1 = PROJ_c[,1],Y2 = PROJ_c[,2],color = as.factor(annotations), batch = as.factor(batch.id))
# plot_proj_df_true_dc<-data.frame(Y1 = PROJ_c_dc[,1],Y2 = PROJ_c_dc[,2],color = as.factor(annotations), batch = as.factor(batch.id))

