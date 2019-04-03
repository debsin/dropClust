## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------

library(dropClust)
# -------------------------------------
# specify paths and load functions
# -------------------------------------
WORK_DIR = "E:/Projects/dropClust/"
DATA_DIR <- file.path(WORK_DIR,"data")       
FIG_DIR <-  file.path(WORK_DIR,"plots")         
REPORT_DIR  <- file.path(WORK_DIR,"report")    

dir.create(file.path(FIG_DIR),showWarnings = FALSE)
dir.create(file.path(REPORT_DIR),showWarnings = FALSE)

set.seed(0)

## ------------------------------------------------------------------------
# Load Data

pbmc.data<-read_data(file.path(DATA_DIR,'pbmc3k/hg19'), format = "10X")


## ------------------------------------------------------------------------
# Filter poor quality cells.  A threshold th corresponds to the total count of a cell.
filtered.data <- filter_cells(pbmc.data)

dim(filtered.data$mat)


## ------------------------------------------------------------------------
print("Fetch rare genes...")
rare_data <- get_rare_genes(filtered.data)




## ------------------------------------------------------------------------
# Filter poor genes
# Genes with UMI count greater than min.count = 2 in atleast min.cell = 3 cells is retained.
lnorm<-normalize(filtered.data, min.count=2, min.cell=3)


# Select Top Dispersed Genes by setting ngenes_keep.
dp_genes <- dispersion_genes(lnorm, ngenes_keep = 1000)


# Log Normalize Matrix with genes-subset, 
# perform batch effect removal operation when input contains batches
whole <- matrix.transform(lnorm,dp_genes ,rare_data, batch_correction = FALSE)



## ------------------------------------------------------------------------

sample_ids = initial.samples(filtered.data, rare_data)


# Structure preserving Sampling

samples_louvain<-sampling(whole[sample_ids,])
subsamples_louvain<-sample_ids[samples_louvain]



## ------------------------------------------------------------------------


# Find PCA top 200 genes. This may take some time.
top_pc_genes<-pc_genes(whole[subsamples_louvain,],top=200)


## ------------------------------------------------------------------------

# Adjust Minimum cluster size with argument minClusterSize (default = 20)
# Adjust tree cut with argument level deepSplit (default = 3), higher value produces more clusters.

clust.list<-cluster.cells(data = whole[,top_pc_genes], sp.samples = subsamples_louvain,
                          default = FALSE, minClusterSize = 30,deepSplit = 2, conf = 0.8)


## ------------------------------------------------------------------------
# ----------------------------
# Tsne & CO-ORDINATE Projection
# ----------------------------
PROJ <- compute_2d_embedding(data = whole[,top_pc_genes],
                             sp.samples = subsamples_louvain,
                             clust.list = clust.list)


# Scatter Plot

plot_proj_df_pred<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(clust.list$cluster.ident))

p<-all_plot(plot_proj_df_pred,filename = NA, title = "dropClust clusters")

print(p)



## ------------------------------------------------------------------------

# Construct sub matrix for DE analysis
de.mat<- reduce_mat_de(lnorm,clust.list)



# Pick Cell Type Specific Genes
#############################

# Cells of interest
GRP = levels(clust.list$cluster.ident)
# int_cells  = which(label %in% GRP)

DE_genes_nodes_all  <- DE_genes(de_data = de.mat, selected_clusters = GRP, lfc_th = 1,q_th = 0.001)

DE_genes_nodes_all$genes

write.csv(DE_genes_nodes_all$genes, 
          file = file.path(REPORT_DIR, "ct_genes.csv"),
          quote = FALSE)

## ------------------------------------------------------------------------

marker_genes = c("S100A8","GNLY","PF4" )

p<-plot_markers(de_data = de.mat, marker_genes)



## ------------------------------------------------------------------------
# Draw heatmap
#############################
p<-plot_heatmap(de_data = de.mat, DE_res = DE_genes_nodes_all$DE_res,nDE = 10)

print(p)



