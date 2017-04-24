rm(list=ls()) # clear workspace

# -------------------------------------
# specify paths and load functions
# -------------------------------------
setwd("~/Projects/dropClust/")
DATA_DIR <- file.path(getwd(),"data/")         # SPECIFY HERE
FIG_DIR <-  paste0(getwd(),"/plots/")        # SPECIFY HERE
REPORT_DIR  <- paste0(getwd(),"/report/")       # SPECIFY HERE
LOUVAIN_DIR <- paste0(getwd(),"/louvain/")   # SPECIFY HERE

dir.create(file.path(FIG_DIR),showWarnings = F)
dir.create(file.path(REPORT_DIR),showWarnings = F)

# ----------------------------
# load relevant libraries
# ----------------------------
source("libraries.R") 
source("ext_functions.R") 

#############
# Load Data###
##############

#pbmc_68k<-read10X("~/Projects/filtered_matrices_mex/hg19/")
#m<-t(pbmc_68k)

pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
true_cls_id <- readRDS(file.path(DATA_DIR,'annotations_68k.rds'))

no_samples = length(pbmc_68k$all_data[[1]]$hg19$barcodes)

#Subsample for LSH training
set.seed(0)
i=min(20000, round(no_samples/3))
sample_ids = sample(1:no_samples, round(i))

# ----------------------------
# DropClust Sampling
# ----------------------------

subsamples_louvain<-dropClust_sampling(pbmc_68k$all_data[[1]]$hg19, sample_ids,true_cls_id)



# ----------------------------
# Find PC Genes
# ----------------------------

m_filt<-readMM("sub_matrix")

cat("FIND PCA top genes...\n")
system.time(top_pc_genes<-pc_genes(m_filt[subsamples_louvain,],200))

write.csv(x = top_pc_genes, file = "pc_gene_ids_sub.csv",quote = F,row.names =F)



# ----------------------------
# Hierarchical Clustering on subsamples
# ----------------------------

ss_sel_genes_mat<-as.matrix(m_filt[subsamples_louvain,top_pc_genes])
ss_clusters<-ss_clustering(ss_sel_genes_mat)


# ----------------------------
# Find K Nearest Neighbours among sub-samples
# ----------------------------
s = Sys.time()
system("python lsh/proj_neigh.py")
cat(paste("Projection LSH Time:",difftime(Sys.time(),s,units = "mins"),"...\n"))



# ----------------------------
# UNANNOTATED Class Assignment 
# ----------------------------

cat("Cast Cluster Ids...\n")
s=Sys.time()
INDEX  = read.csv("neigh.txt",header=F,sep=" ")
dim(INDEX)
INDEX = as.matrix(INDEX)+1

clust_col = list()

for( i in 1:nrow(INDEX))
{
  clust_col  = list(clust_col ,names(which.max(table(ss_clusters$labels[INDEX[i,]]))))
}
clust_col = as.numeric(unlist(clust_col))
cat(paste("Time:",difftime(Sys.time(),s,units = "mins"),"...\n"))
write.csv(x = clust_col,file ="predicted.csv", quote = F,row.names = F)
cat("Predicted Metric HCLUST:\n")
print(sc_metric(clust_col, as.numeric(true_cls_id),show_tab = F))


# ----------------------------
# Tsne & CO-ORDINATE Projection 
# ----------------------------

set.seed(0)
cat("Compute Tsne Projection using PCA top genes...\n")
system.time(ts<-Rtsne(as.matrix(m_filt[subsamples_louvain,top_pc_genes]),perplexity = 20,dims = 2))

cat("Projecting TSNE co-ordinates...\n")
s=Sys.time()
tt<-ts$Y[-ss_clusters$outliers,]
PROJ = matrix(NA, nrow=nrow(INDEX), ncol=2)

for( i in 1:nrow(INDEX))
{
  maj_id = clust_col[i]
  maj_col = which(ss_labels[INDEX[i,]]==maj_id)
  
  if(length(maj_col)==1){
    PROJ[i,] = Matrix::colMeans(tt[INDEX[i,],])
  }
  else{
    PROJ[i,] = Matrix::colMeans(tt[INDEX[i,maj_col],])
  }

}

dim(PROJ)
cat(paste("Projection Time:",difftime(Sys.time(),s,units = "mins"),"...\n"))

dropClust_df<-as.data.frame(cbind(PROJ, clust_col))
rownames(dropClust_df)<-pbmc_68k$all_data[[1]]$hg19$barcodes
save(dropClust_df,file="proj_68.Rda")


####################
# TRUE CLUSTERS
####################
plot_proj_df<-data.frame("Y1" = PROJ[,1],"Y2" = PROJ[,2],color =as.factor(true_cls_id))
plot_proj_df$color <- 
  factor(plot_proj_df$color, levels = c( "CD8+ Cytotoxic T","CD8+/CD45RA+ Naive Cytotoxic",
                                         "CD4+/CD25 T Reg", 
                                        "CD14+ Monocyte", "CD34+", 
                                        "Dendritic","CD19+ B", "CD4+ T Helper2",
                                        "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory","CD56+ NK"))
filename<-paste0(FIG_DIR,"68_true_proj.pdf")

all_plot(plot_proj_df,filename,"TRUE COLOR PROJECTION")


# MAJORITY VOTING CLUSTERS
##########################
plot_proj_df_pred<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(clust_col))
filename<-paste0(FIG_DIR,"68_pred_proj.pdf")
all_plot(plot_proj_df_pred,filename,"HCLUST Majority Voting clusters ")

