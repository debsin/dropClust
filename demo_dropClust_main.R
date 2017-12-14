

#### Specify paths

sourceDir = getwd() # The directory hosting the current script 
setwd(sourceDir)
DATA_DIR <- file.path(sourceDir,"data")        
FIG_DIR <-  file.path(sourceDir,"plots")        
REPORT_DIR  <- file.path(sourceDir,"report")   
LOUVAIN_DIR <- file.path(sourceDir,"louvain")  

dir.create(file.path(FIG_DIR),showWarnings = F)
dir.create(file.path(REPORT_DIR),showWarnings = F)


# Load relevant libraries & Functions for clustering

suppressMessages(source("libraries.R") )
suppressMessages(source("all_functions.R"))

# Load Data and Annotations, Annotations may be omitted if unavailable.

mouse.data<-read10X(file.path(DATA_DIR,"es_mouse/"))
write.csv(x = mouse.data$gene_symbols, file = "gene_symbols.csv",quote = F,row.names =F)
annotation<- read.table(file.path(DATA_DIR,"annotations.csv"),sep = ',',header = T)
ref_id <- as.factor(annotation$x)
dim(mouse.data$mat)

# Filter poor quality cells.  A threshold th corresponds to the total count of a cell.
filtered.data = filter_cells(mouse.data,th = 5000)

dim(filtered.data$mat)

anno_labels= ref_id[filtered.data$keep_cells]
no_samples = dim(filtered.data$mat)[1]

# Filter poor genes
# Genes with UMI count greater than min.count = 2 in atleast min.cell = 3 cells is retained.
lnorm<-normalize_by_umi_2(filtered.data, min.count=2, min.cell=3)   

# The minimum betwen 20000 and 1/3 of the total number of samples is chosen
i=min(20000, round(no_samples/3))
sample_ids = sample(1:no_samples, i)


# Select Top Dispersed Genes by setting ngenes_keep. 
m_n_ngenes <- matrix.subset(lnorm, sample_ids, ngenes_keep = 1000)

# Structure preserving Sampling
call_lsh()

call_louvain()

#This step allows the user to fix a managbale size of samples for hierarchial clustering. Omit this step to use default values
opt_pinit = optimized_Pinit(nsamples = 500)

# Sub-sampling using obtained parameter
subsamples_louvain<-sampling(pinit = opt_pinit)
write.csv(x = subsamples_louvain, file = "subsamples_idx",quote = F,row.names =F)
write.csv(x = filtered.data$barcodes[sample_ids[subsamples_louvain]], file = "barcodes_subsamples.csv",quote = F,row.names =F)

# Check the distribution of Sampling across annotated cell types. Skip if annotation is unavailable.
table(anno_labels[sample_ids[subsamples_louvain]])

# Find PCA top 200 genes. This may take some time. 
top_pc_genes<-pc_genes(m_n_ngenes[subsamples_louvain,],top=200) 



write.csv(x = top_pc_genes, file = "pc_gene_ids_sub.csv",quote = F,row.names =F)

#### Hierarchical Clustering on subsamples
# Adjust Minimum cluster size with argument minClusterSize (default = 20)
# Adjust tree cut with argument level deepSplit (default = 3), higher value produces more clusters.

ss_sel_genes_mat<-as.matrix(m_n_ngenes[subsamples_louvain,top_pc_genes])
ss_clusters<-ss_clustering(ss_sel_genes_mat, minClusterSize = 20, deepSplit = 3) 

# Find KNN of remaining samples
system("python lsh/proj_neigh.py")


#### Un-annotated Cell Assignment
INDEX  = read.csv("neigh.txt",header=F,sep=" ")
clust_col<-cluster_assign(INDEX, ss_clusters)
write.csv(x = clust_col,file ="predicted.csv", quote = F,row.names = F)

# Skip if annotation is unavailable
sc_metric(clust_col, as.numeric(anno_labels),show_tab = F)


### 2D Vizualization:
#### T-sne & Projection of co-ordinates of all samples in 2D

PROJ = compute_2d_embedding(data = as.matrix(m_n_ngenes[subsamples_louvain,top_pc_genes]), ss_clusters, INDEX)

dropClust_df<-as.data.frame(cbind(PROJ, clust_col))
rownames(dropClust_df)<-filtered.data$barcodes
save(dropClust_df,file="demo_proj.Rda")



# 2D Vizualization: Known References to Days.
plot_proj_df<-data.frame("Y1" = PROJ[,1],"Y2" = PROJ[,2],color =as.factor(anno_labels))
all_plot(plot_proj_df,"known.jpg","Known Day")


# 2D Vizualization: Predicted Clusters
plot_proj_df_pred<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(clust_col))
all_plot(plot_proj_df_pred,"pred.jpg","Predicted clusters")