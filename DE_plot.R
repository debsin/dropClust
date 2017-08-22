rm(list=ls()) # clear workspace


# Specify paths and load functions
# -------------------------------------
setwd("~/Projects/dropClust/")
DATA_DIR <- paste0(getwd(),"/data/")          # SPECIFY HERE
FIG_DIR <-  paste0(getwd(),"/plots/")        # SPECIFY HERE
REPORT_DIR  <- paste0(getwd(),"/report/")       # SPECIFY HERE



source("libraries.R") 
library(gplots)
#-------------
# Load additional Libraries for parallel processing
library(foreach)
library(doParallel)


source("ext_functions.R") 
source("DE_functions.R") 


### WARNING!!!! Check number of cores
registerDoParallel(7)

#Read raw data
pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
# all_data <- pbmc_68k$all_data
# m<-all_data[[1]]$hg19$mat

true_class<-readRDS(file.path(DATA_DIR,"annotations_68k.rds"))

table(true_class)



#Normalize by umi counts (median)
k<-normalize_by_umi_2(pbmc_68k$all_data[[1]]$hg19)  

#save(k, file = "~/Projects/scClust-Final/data/68K_7K.Rda")
#writeMM(k$m,file = "~/Projects/scClust-Final/data/68K_7K.smat")

load("proj_68.Rda")
plot(dropClust_df$V1,dropClust_df$V2,col=dropClust_df$clust_col)

m_n<-k$m

dim(m_n)

#Read predicted cluster IDs
pred_labels = read.csv("predicted.csv", sep="",header = T)
#Sample from each cluster
set.seed(0)
fixed_samples = c()
for(clust_id in levels(as.factor(pred_labels$x))){
  if(clust_id==0) next;
  sample_ids_per_cluster = which(pred_labels$x==clust_id)
  fixed_samples = c(fixed_samples, sample(sample_ids_per_cluster,min(100,length(sample_ids_per_cluster))))
}

label = pred_labels$x[fixed_samples]

#Log transform
sub_samples2 = log2(as.matrix(m_n[fixed_samples,])+1)
rownames(sub_samples2)  = pbmc_68k$all_data[[1]]$hg19$barcodes[fixed_samples]
colnames(sub_samples2)  = k$use_genes

dim(sub_samples2)
MAT123 = t(sub_samples2)
MAT123 = MAT123[match(unique(rownames(MAT123)),rownames(MAT123)),]
dim(MAT123)

save(sub_samples2, file = "subsamples.Rdata")
save(label, file= "pred_label.Rdata")




#Pick Cell Type Specific Genes
#############################
ID_all = pred_labels$x[fixed_samples]


# Cells of interest
# GRP1 = c(3,5,6,7,13)
# GRP2 = c(1,2,10)
# GRP3 = c(8,9,11,14)
# GRP4 = c(5,7)
GRP = c(1:max(pred_labels))

int_cells  = which(ID_all %in% GRP)

ID = ID_all[int_cells]
 


pair_id ="all_pairs"
s = Sys.time()
DE_genes_nodes_all  <- DE_genes(data = MAT123[,int_cells] ,labels = label[int_cells],max = 0,lfc_th = 0.26,q_th = 0.05) #max=0 for all DE genes, else top |max| genes
cat(paste("DE Time:", difftime(Sys.time(),s, units = "mins"),"...\n"))


Mat_ct  = MAT123[DE_genes_nodes_all[["genes"]],int_cells]

sink(file.path(REPORT_DIR,"DE_list.txt"))
DE_genes_nodes_all[["DE_res"]]
sink()

all_ct_genes = c()
ct_genes_list = c()

for(i in unique(ID))
{
  
  ct_genes = c()
  for(j in unique(ID))
  {
    
    if(i!=j)
    {
      
      
      #j=2
      
      # indices of groups
      id_1 = which(ID == i)
      id_2 = which(ID == j)
      pair_id1 = paste0(i,"_",j)
      pair_id2 = paste0(j,"_",i)
      #pair_id = union()
      # comparing the avg expression - or the direction
      diff_i_j = union(DE_genes_nodes_all[["DE_res"]][[pair_id2]]$gene,DE_genes_nodes_all[["DE_res"]][[pair_id1]]$gene)
      
      RM_i = rowMeans(Mat_ct[diff_i_j,id_1])
      RM_j = rowMeans(Mat_ct[diff_i_j,id_2])
      
      # direction
      g = rownames(Mat_ct[diff_i_j,])[which(RM_i > RM_j)]
      ct_genes = append(ct_genes, g)
      
    }
  }
  
  #ct_genes_list[[paste0(i)]] = c(ct_genes_list[[paste0(i)]], ct_genes)
  
  print(paste(pair_id1,pair_id2,length(names(table(ct_genes))[which(table(ct_genes) == (length(unique(ID))-1))])))
  all_ct_genes = append(all_ct_genes, names(table(ct_genes))[which(table(ct_genes) == (length(unique(ID))-1))])
  ct_genes_list[[paste0(i)]] = c(ct_genes_list[[paste0(i)]], names(table(ct_genes))[which(table(ct_genes) == (length(unique(ID))-1))])
  
  
  
}

length(all_ct_genes)
meth= "DE"
heat_in = MAT123[all_ct_genes,int_cells]
heat_obj<-plot_heat_genes(data = heat_in,labels = label[int_cells],meth)


#Uncomment to reorder label
#ordered_labl = heat_obj[["lab"]]
ordered_labl = label[int_cells]

colors = getColors(length(unique(label)))
colors = colors[ordered_labl]
ordered_freq = table(factor(ordered_labl,levels=unique(ordered_labl)))
a = cumsum(ordered_freq)
b = c(0,  a[-length(a)])

pos = round((a+b)/2)

text_lab = rep(NA,dim(heat_in)[2])

text_lab[pos] = unique(ordered_labl);


myPalette <-colorRampPalette(brewer.pal(11,"RdBu"))
myPalette <- colorRampPalette(c("grey15","gray10","black","yellow","darkred"))

filename<-file.path(FIG_DIR,paste0(paste0(GRP,collapse = "_"),"heatmap.pdf"))
pdf(filename,width = 10, height = 15)

#heat_obj[["mat"]]

hm<-heatmap.2(heat_in, 
              trace="none", col=myPalette(40),Colv = FALSE, Rowv = FALSE, srtCol=0, 
              ColSideColors=colors, cexCol = 1,cexRow = 0.2,labCol = text_lab,
              dendrogram="none",scale="row",density.info="none",keysize = 0.5,
              main=paste(meth, "genes and equal-samples"))



dev.off()

library(plyr) 

n.obs <- sapply(ct_genes_list, length)
seq.max <- seq_len(max(n.obs))
mat <- as.data.frame(t(sapply(ct_genes_list, "[", i = seq.max)))
gene.df <- sapply(mat, as.character)
gene.df[is.na(gene.df)] <- ""
write.csv(gene.df, file = file.path(REPORT_DIR, "ct_genes.csv"),quote = F)


