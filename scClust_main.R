rm(list=ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
source("libraries.R") 
# -------------------------------------
# specify paths and load functions
# -------------------------------------
setwd("~/Projects/scClust-Final/")
DATA_DIR <- "~/Projects/10Xpaper/"           # SPECIFY HERE
FIG_DIR <-  paste0(getwd(),"/plots/")        # SPECIFY HERE
REPORT_DIR  <- paste0(getwd(),"/report/")       # SPECIFY HERE
LOUVAIN_DIR <- paste0(getwd(),"/louvain/")   # SPECIFY HERE

dir.create(file.path(FIG_DIR),showWarnings = F)
dir.create(file.path(REPORT_DIR),showWarnings = F)

source("ext_functions.R") 

#############
# Load Data###
##############

#pbmc_68k<-read10X("~/Projects/filtered_matrices_mex/hg19/")
#m<-t(pbmc_68k)

 pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
 all_data <- pbmc_68k$all_data
 m<-all_data[[1]]$hg19$mat
 write.csv(x = all_data[[1]]$hg19$gene_symbols, file = "gene_symbols.csv",quote = F,row.names =F)
true_cls_id <- readRDS(file.path(DATA_DIR,'annotations_68k.rds'))
no_samples = dim(m)[1]

#Subsample for LSH training
set.seed(0)
i=min(20000, round(no_samples/3))
sample_ids = sample(1:no_samples, round(i))

###############################
######## RUN LSH  #############
# ASSIGN PREDICTED CLUSTER IDS#
# RECORD RAND INDEX ###########
###############################


s = Sys.time()
clus_metric <- method_scClust(all_data[[1]]$hg19, sample_ids,true_cls_id)
cat(paste("Time:", difftime(Sys.time(),s, units = "mins"),"...\n"))
cat("1st Level Cluster Metric:\n")

output <- read.csv("output.graph", sep="",header = F)
sc_metric(output$V2+1, as.numeric(true_cls_id)[sample_ids])
print(clus_metric)




########################################
## READ LOUVEN CLUSTERS FOR SUB-SAMPLING
########################################
#set.seed(0)
subsamples_louvain<-sampling()

cat(paste("number of sub-samples:",length(subsamples_louvain),"\n"))

write.csv(x = subsamples_louvain, file = "subsamples_idx",quote = F,row.names =F)
write.csv(x = all_data[[1]]$hg19$barcodes[sample_ids[subsamples_louvain]], file = "barcodes_subsamples.csv",quote = F,row.names =F)

table(true_cls_id[sample_ids[subsamples_louvain]])
m_filt<-readMM("sub_matrix")

#=================================
# CREATE SAMPLING PLOT DATA FRAME
##################################

louvain_samples = true_cls_id[sample_ids[subsamples_louvain]]
random_samples = sample(true_cls_id,length(subsamples_louvain))
table(random_samples)
original = true_cls_id[sample_ids]

plot_bars<-rbind(random=table(random_samples)/table(original),
                 louv = table(louvain_samples)/table(original))
colnames(plot_bars)<-paste(colnames(plot_bars)," (",table(original),")",sep="")
plot_bars<-plot_bars[,order(table(original))]
#View(t(plot_bars))
plot_bars<-melt(plot_bars)
names(plot_bars)[names(plot_bars) == 'value'] <- 'Proportion'
names(plot_bars)[names(plot_bars) == 'Var1'] <- 'Sampling'


png(file.path(FIG_DIR, "sampling_proportion.png"))
ggplot(plot_bars,aes(x = Var2, y = Proportion,fill=Sampling)) + geom_bar(position="dodge",stat="identity") + ggtitle("Sampling Grouped Barplot")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 8),axis.text.y = element_text(size = 8), axis.title.x=element_blank())+ 
  scale_fill_brewer(palette = "Set1")
dev.off()

##############
## PCA #######
##############
#Top Thousand dispersed genes_id after colsum removed
genes2K = read.csv("genes", sep="",header = T)
cat("FIND PCA top genes...\n")
system.time(top_pc_genes<-pc_genes(m_filt[subsamples_louvain,],200))

genes_all = read.csv("genes_used_all", sep="",header = T)
#write.csv(x = all_data[[1]]$hg19$gene_symbols[genes_all$x[top_pc_genes]], file = "pc_gene_symbols.csv",quote = F,row.names =F)
write.csv(x = genes_all$x[top_pc_genes], file = "pc_gene_symbols.csv",quote = F,row.names =F)
write.csv(x = top_pc_genes, file = "pc_gene_ids_sub.csv",quote = F,row.names =F)
write.csv(x = genes2K$x[top_pc_genes], file = "pc_gene_ids_all.csv",quote = F,row.names =F)

dim(m_filt)
sub_samples<-as.matrix(m_filt[subsamples_louvain,top_pc_genes])
dim(sub_samples)
rownames(sub_samples)<-all_data[[1]]$hg19$barcodes[sample_ids[subsamples_louvain]]
#colnames(sub_samples)<-all_data[[1]]$hg19$gene_symbols[genes_all$x[top_pc_genes]]
colnames(sub_samples)<-genes_all$x[top_pc_genes]

save(sub_samples,file="sub_samples.Rda")
######################
## TSNE ##############
######################
set.seed(0)
cat("Compute Tsne Projection using PCA top genes...\n")
system.time(ts<-Rtsne(as.matrix(m_filt[subsamples_louvain,top_pc_genes]),perplexity = 20))

plot(ts$Y)


####################
### HCLUST #########
####################

ss_sel_genes<-as.matrix(m_filt[subsamples_louvain,top_pc_genes])
save(ss_sel_genes, file="ss_pc_mat.Rdata")
corr_mat = cor(t(ss_sel_genes))
d = dist(ss_sel_genes)
hc<-fastcluster::hclust(d,method = "average")

hc_labs<-cutreeDynamic(dendro = hc, cutHeight = NULL,
                       minClusterSize = 10,
                       method = "hybrid", deepSplit = 3,
                       pamStage = TRUE,  distM = as.matrix(d), maxPamDist = 0,
                       verbose = 0)

cat(paste("Number of Clusters:", length(unique(hc_labs))-1),"\n")

outliers_ids <- which(hc_labs==0)
subsamples_louv_68K = sample_ids[subsamples_louvain[-outliers_ids]]
write.csv(x = subsamples_louv_68K, file = "subsamples_louv_68K",quote = F,row.names =F)
write.csv(x = all_data[[1]]$hg19$barcodes[subsamples_louv_68K], file = "subsamples_barcodes",quote = F,row.names =F)

write.csv(x = hc_labs, file = "hc_labs.csv",quote = F,row.names =F)

cat("Sub-sample Metric\n")
print(sc_metric(hc_labs[-outliers_ids], true_cls_id[subsamples_louv_68K],show_tab = F))


hc_labs_clean  = hc_labs[-outliers_ids]
cat(paste("Predicted Clusters:",length(unique(hc_labs_clean))))


plot_df_ss<-data.frame(Y1 = ts$Y[-outliers_ids,1],Y2 = ts$Y[-outliers_ids,2],color = as.factor(hc_labs_clean))
colorcount = length(unique(plot_df_ss$color))
pred_colors = getColors(colorcount)
ggplot(plot_df_ss,aes(Y1,Y2,col=color))+geom_point(size=0.8)+scale_colour_manual(values = pred_colors)+ggtitle("HC Prediction") + guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))
filename<-file.path(FIG_DIR,paste0(gsub(":","_",Sys.time()),"_ss_pred_proj.png"))
png(filename,width = 1000,height = 800)
p<-ggplot(plot_df_ss,aes(Y1,Y2,col=color))+geom_point(size=0.8)+scale_colour_manual(values = pred_colors)+ggtitle("HC Prediction") + guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))
print(p)
dev.off()




filename<-file.path(FIG_DIR,paste0(gsub(":","_",Sys.time()),"_ss_true_proj.png"))
png(filename,width = 1000,height = 800)
plot_df_ss<-data.frame(Y1 = ts$Y[-outliers_ids,1],Y2 = ts$Y[-outliers_ids,2],color = as.factor(true_cls_id[subsamples_louv_68K]))
colorcount_ts = length(unique(plot_df_ss$color))
ggplot(plot_df_ss,aes(Y1,Y2,col=color))+geom_point(size=0.8)+scale_colour_manual(values = getColors(colorcount_ts))+ggtitle("TRUE COLORS")+ guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))

p<-ggplot(plot_df_ss,aes(Y1,Y2,col=color))+geom_point(size=0.8)+scale_colour_manual(values = getColors(colorcount_ts))+ggtitle("TRUE COLORS")+ guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))
print(p)
dev.off()
#ggplot(plot_df_ss,aes(Y1,Y2,col=color))+geom_point(size=0.8)+scale_colour_manual(values = getPalette(colorcount_ts))+ggtitle("TRUE COLORS")+ guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))



s = Sys.time()
system("python lsh/proj_neigh.py")
cat(paste("Projection LSH Time:",difftime(Sys.time(),s,units = "mins"),"...\n"))

#################################
## UNANNOTATED Class Assignment #
## NN MAJORITY VOTING, USE ######
## PREDICTED SUB-SAMPLES ########
#################################
cat("Cast Cluster Ids...\n")
s=Sys.time()
INDEX  = read.csv("neigh.txt",header=F,sep=" ")
dim(INDEX)
INDEX = as.matrix(INDEX)+1

clust_col = list()
#louv_subsample_labels  = output$V2[subsamples_louvain]

for( i in 1:nrow(INDEX))
{
  clust_col  = list(clust_col ,names(which.max(table(hc_labs_clean[INDEX[i,]]))))
}
clust_col = as.numeric(unlist(clust_col))
cat(paste("Time:",difftime(Sys.time(),s,units = "mins"),"...\n"))
write.csv(x = clust_col,file ="predicted.csv", quote = F,row.names = F)
cat("Predicted Metric HCLUST...\n")
print(sc_metric(clust_col, as.numeric(true_cls_id),show_tab = F))

##############################
## CO-ORDINATE Projection ####
##############################
cat("Projecting TSNE co-ordinates...\n")
s=Sys.time()
tt<-ts$Y[-outliers_ids,]
PROJ = matrix(NA, nrow=nrow(INDEX), ncol=2)
#maj_count = c()
for( i in 1:nrow(INDEX))
{
  maj_id = clust_col[i]
  maj_col = which(hc_labs_clean[INDEX[i,]]==maj_id)
  #maj_count = c(maj_count,length(maj_col ))
  if(length(maj_col)==1){
    PROJ[i,] = Matrix::colMeans(tt[INDEX[i,],])
  }
  else{
    PROJ[i,] = Matrix::colMeans(tt[INDEX[i,maj_col],])
  }

}
#PROJ = as.matrix(data.frame(PROJ))
dim(PROJ)
cat(paste("Projection Time:",difftime(Sys.time(),s,units = "mins"),"...\n"))



####################
# TRUE CLUSTERS
####################
plot_proj_df<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color =as.factor(true_cls_id))
plot_proj_df$color <- 
  factor(plot_proj_df$color, levels = c( "CD8+ Cytotoxic T","CD8+/CD45RA+ Naive Cytotoxic",
                                         "CD4+/CD25 T Reg", 
                                        "CD14+ Monocyte", "CD34+", 
                                        "Dendritic","CD19+ B", "CD4+ T Helper2",
                                        "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory","CD56+ NK"))

x.mean = aggregate(plot_proj_df$Y1, list(plot_proj_df$color), median)[,-1]
y.mean = aggregate(plot_proj_df$Y2, list(plot_proj_df$color), median)[,-1]

colorcount_t = length(unique(plot_proj_df$color))
filename<-paste0(FIG_DIR,gsub(":","_",Sys.time()),"_68_true_proj.png")
png(filename,width = 1000,height = 800)

p<-ggplot(plot_proj_df,aes(Y1,Y2,col= color))
p2<-p+ geom_point(size=0.3)  + scale_colour_manual(values =  getColors(colorcount_t))+
  ggtitle("TRUE COLOR PROJECTION")+
  annotate("text", x = x.mean, y = y.mean, label = levels(plot_proj_df$color), size =3.5 )+
   guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))

print(p2)

dev.off()

# MAJORITY VOTING CLUSTERS
##########################
plot_proj_df_pred<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(clust_col))
x.mean = aggregate(plot_proj_df_pred$Y1, list(plot_proj_df_pred$color), median)[,-1]
y.mean = aggregate(plot_proj_df_pred$Y2, list(plot_proj_df_pred$color), median)[,-1]
filename<-paste0(FIG_DIR,gsub(":","_",Sys.time()),"_68_pred_proj.png")
png(filename,width = 1000,height = 800)
p<-ggplot(plot_proj_df_pred,aes(Y1,Y2,col = color))+geom_point(size=0.3,alpha=0.95)+scale_colour_manual(values =  pred_colors)+
  annotate("text", x = x.mean, y = y.mean, label = levels(plot_proj_df_pred$color), size =4 )+
  ggtitle("HCLUST Majority Voting clusters ")+ guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))
print(p)
dev.off()
#ggplot(plot_proj_df_pred[which(plot_proj_df_pred$color==5),],aes(Y1,Y2,col = color))+geom_point(size=0.4)+scale_colour_manual(values = sample(mycolors))



save(PROJ,file="proj_68.Rda")

###############
####### PRED
###############
library(gridExtra)
plots=list()
plot_proj_df_pred<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(clust_col))
for(i in 1:length(levels(plot_proj_df_pred$color))){
  plot_proj_df_pred<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(clust_col))
  one = as.character(unique(plot_proj_df_pred$color)[i])
  plot_proj_df_pred$color = ifelse(plot_proj_df_pred$color == one,one,"grey")
  pl = ggplot(plot_proj_df_pred,aes(Y1,Y2,col = color))+geom_point(alpha = 0.2, size=0.4)+scale_colour_manual(values = c("red","grey"))+ theme(legend.position="none")+ggtitle(paste("PRED clusters",one))
  plots[[i]] <- pl
}

layout <- matrix(c(1:8), nrow = 2, byrow = TRUE)

png(file.path(FIG_DIR,paste0("9_grid_68_pred_1:8.png")),width = 1200,height = 1200)
do.call("grid.arrange", c(plots[1:8], ncol=4))
dev.off()
png(paste0(FIG_DIR,"9_grid_68_pred_9:16.png"),width = 1200,height = 1200)
do.call("grid.arrange", c(plots[9:min(16,length(plots))], ncol=4))
dev.off()



table(clust_col)



###############
####### TRUE
###############
library(gridExtra)
plots=list()
plot_proj_df_true<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(true_cls_id))
for(i in 1:length(levels(plot_proj_df_true$color))){
plot_proj_df_true<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(true_cls_id))
one = as.character(unique(plot_proj_df_true$color)[i])
plot_proj_df_true$color = ifelse(plot_proj_df_true$color == one,one,"grey")
pl = ggplot(plot_proj_df_true,aes(Y1,Y2,col = color))+geom_point(alpha = 0.2, size=0.4)+scale_colour_manual(values = c("red","grey"))+theme(legend.position = "none")+ggtitle(one)
plots[[i]] <- pl
}

png(paste0(FIG_DIR,"9_grid_68_true.png"),width = 1200,height = 1200)
do.call("grid.arrange", c(plots[1:11], ncol=4,nrow = 3))
dev.off()
