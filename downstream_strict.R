rm(list=ls()) # clear workspace

# specify paths and load functions
# -------------------------------------
setwd("~/Projects/scClust-Final/")
DATA_DIR <- "~/Projects/10Xpaper/"           # SPECIFY HERE
FIG_DIR <-  paste0(getwd(),"/plots")        # SPECIFY HERE
REPORT_DIR  <- paste0(getwd(),"/report")       # SPECIFY HERE
dir.create(file.path(FIG_DIR, "heatmap"),showWarnings = F)
dir.create(file.path(REPORT_DIR),showWarnings = F)

### WARNING!!!! Check number of cores
registerDoParallel(60)


source("libraries.R") 
source("ext_functions.R") 
source("DE_functions.R") 


pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
all_data <- pbmc_68k$all_data
m<-all_data[[1]]$hg19$mat

x = all_data[[1]]$hg19


k<-.normalize_by_umi_2(x)   

m_n<-k$m

dim(m_n)

pred_labels = read.csv("predicted.csv", sep="",header = T)

fixed_samples = c()
for(clust_id in unique(pred_labels$x)){
  if(clust_id==0) next;
  sample_ids_per_cluster = which(pred_labels$x==clust_id)
  fixed_samples = c(fixed_samples, sample(sample_ids_per_cluster,min(100,length(sample_ids_per_cluster))))
}


sub_samples2 = log2(as.matrix(m_n[fixed_samples,])+1)
rownames(sub_samples2)  = all_data[[1]]$hg19$barcodes[fixed_samples]
colnames(sub_samples2)  = k$use_genes

dim(sub_samples2)
MAT123 = t(sub_samples2)
MAT123 = MAT123[match(unique(rownames(MAT123)),rownames(MAT123)),]
dim(MAT123)

label = pred_labels$x[fixed_samples]

save(sub_samples2, file = "subsamples.Rdata")
save(label, file= "pred_label.Rdata")

s = Sys.time()
DE_genes_nodes_all  <- DE_genes(data = MAT123 ,labels = label,max = 20) #max=0 for all DE genes, else top |max| genes
cat(paste("DE Time:", difftime(Sys.time(),s, units = "mins"),"...\n"))

ncusters = length(unique(label))
Count_mat = matrix(rep(0,ncusters * ncusters),nrow = ncusters,ncol= ncusters )

for(i in names(DE_genes_nodes_all[["DE_res"]])){
  
  coord = unlist(strsplit(x =i,split = "_"))
  row = as.integer(coord[1])
  col= as.integer(coord[2])
  #print(paste(row,col))
  pair_mat = DE_genes_nodes_all[["DE_res"]][[i]]
  genes = length(rownames(pair_mat[which(pair_mat[,2]<=0.05),]))
  Count_mat[row,col ] = genes
  Count_mat[col,row ] = genes
  
}
Count_mat

meth= "DE"
heat_in = MAT123[DE_genes_nodes_all[["genes"]],]
heat_obj<-plot_heat_genes(data = heat_in,labels = label,meth)

############################
### P C A GENES ONLY #######
############################
meth= "PC"
pc_genes = read.csv("pc_gene_symbols.csv")
louv_samples = read.csv("barcodes_subsamples.csv")

ss_ids = which(all_data[[1]]$hg19$barcodes %in% louv_samples$x) #all_data[[1]]$hg19$barcodes
pc_ids =  which( k$use_genes %in% pc_genes$x)

sub_samples_pc = log2(as.matrix(m_n[ss_ids,pc_ids])+1)
dim(sub_samples_pc)
rownames(sub_samples_pc)  = louv_samples$x #all_data[[1]]$hg19$barcodes
colnames(sub_samples_pc)  = pc_genes$x

label_ss = pred_labels$x[ss_ids];

fixed_samples = c()
for(clust_id in unique(label_ss)){
  if(clust_id==0) next;
  sample_ids_per_cluster = which(label_ss==clust_id)
  fixed_samples = c(fixed_samples, sample(sample_ids_per_cluster,min(100,length(sample_ids_per_cluster))))
}

label = label_ss[fixed_samples]
sub_samples_pc = sub_samples_pc[fixed_samples, ]
dim(sub_samples_pc)
heat_in = t(sub_samples_pc) ; 
heat_obj<-plot_heat_genes(data = heat_in,labels = label,meth)

############################
### P C A GENES END ########
############################

# 
# selids = which(label %in% c(1,2))
# heat_in =MAT123[rownames(DE_genes_nodes_all[["DE_res"]][["1_2"]]),selids]
# heat_obj<-plot_heat_genes(data = heat_in,labels = label[selids],"DE")

#myPalette <-colorRampPalette(c("darkblue","blue", "white","red","darkred"))
myPalette <-colorRampPalette(brewer.pal(11,"RdBu"))


ordered_labl = heat_obj[["lab"]]

colors = getColors(length(unique(label)))
colors = colors[ordered_labl]
ordered_freq = table(factor(ordered_labl,levels=unique(ordered_labl)))
a = cumsum(ordered_freq)
b = c(0,  a[-length(a)])

pos = round((a+b)/2)

text_lab = rep(NA,dim(heat_in)[2])

text_lab[pos] = unique(ordered_labl);



filename<-file.path(FIG_DIR,"heatmap",paste(gsub(":","_",Sys.time()),meth,"genes_heatmap.pdf",sep="_"))
pdf(filename)

hm<-heatmap.2(heat_obj[["mat"]], 
          trace="none", col=myPalette,Colv = FALSE,  srtCol=0, 
          ColSideColors=colors, cexCol = 1,cexRow = 0.9,labCol = text_lab,
          dendrogram="none",scale="none",density.info="none",keysize = 1,
          colsep=a[1:(length(unique(colors))-1)], sepcolor="green", sepwidth = 0.5, main=paste(meth, "genes and All-samples"))



dev.off()

write.csv(hm[["carpet"]],file = file.path(REPORT_DIR,"heatmap_DE_top20.csv"))

#
#
# #
# 
# load("sub_samples.Rda")
# 
# 
# pred_ids = which(all_data[[1]]$hg19$barcodes %in% rownames(sub_samples))
# 
# 
# 
# pred_labels = read.csv("predicted.csv", sep="",header = T)
# 
# label = pred_labels$x[pred_ids]
# 
# fixed_samples = c()
# for(clust_id in unique(label)){
#   if(clust_id==0) next;
#   sample_ids_per_cluster = which(label==clust_id)
#   fixed_samples = c(fixed_samples, sample(sample_ids_per_cluster,min(300,length(sample_ids_per_cluster))))
# }
# 
# 
# 
# MAT1234 = t(sub_samples[fixed_samples,])
# MAT1234 = MAT1234[match(unique(rownames(MAT1234)),rownames(MAT1234)),]
# dim(MAT1234)
# 
# 
# 
# heat_obj_pc<-plot_heat_genes(data = MAT1234,labels = label[fixed_samples],"PC")
# 
# 
# ordered_labl = heat_obj_pc[["lab"]]
# 
# colors = getColors(length(unique(ordered_labl)))
# colors = colors[ordered_labl]
# ordered_freq = table(factor(ordered_labl,levels=unique(ordered_labl)))
# a = cumsum(ordered_freq)
# b = c(0,  a[-length(a)])
# 
# pos = round((a+b)/2)
# 
# text_lab = rep(NA,dim(MAT1234)[2])
# 
# text_lab[pos] = unique(ordered_labl)
# 
# myPalette <-colorRampPalette(c("darkblue","blue", "white","red","darkred"))
# heatmap.2(heat_obj_pc[["mat"]], 
#           trace="none", col=myPalette,Colv = FALSE,  srtCol=0, 
#           ColSideColors=colors, cexCol = 1,cexRow = 0.9,labCol = text_lab,
#           dendrogram="none",scale="row",density.info="none",keysize = 1,
#           colsep=a[1:(length(unique(colors))-1)], sepcolor="green")
# 
# 
