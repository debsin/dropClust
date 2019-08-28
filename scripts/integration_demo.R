
load("C:/Projects/CellBench/CellBench_data-master/data/sincell_with_class.RData")

load("~/GitHub/dropclust_benchmarks/analysis/bc/Cellbench/out_gene_ids.Rda")

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

mix.data<-Merge(objects)


annotations = mix.data$cell_line#c(colData(objects[[1]])$cell_line_demuxlet, colData(objects[[2]])$cell_line_demuxlet, colData(objects[[3]])$cell_line_demuxlet)
batch.id = mix.data$Batch


dc.corr <- Integration(mix.data,close_th = 0.1, cells_th = 0.1,
                       components = 10, n_neighbors = 40, init = "spca", spread = 1)


cc = Cluster(dc.corr,method = "kmeans",centers = 3)


de<-FindMarkers(cc,q_th = 0.001, lfc_th = 1.2)

source("~/GitHub/dropclust_benchmarks/analysis/bc/Cellbench/plots.R")

PROJ_c_dc = reducedDim(cc, "iComponents")
plot_proj_df_true_dc<-data.frame(Y1 = PROJ_c_dc[,1],Y2 = PROJ_c_dc[,2], color = as.factor(cc$ClusterIDs), batch = as.integer(cc$Batch))
batch_plot(plot_proj_df_true_dc,filename = NA, title = "Corrected DC ",  type=NULL)

