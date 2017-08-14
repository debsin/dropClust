# ---------------------------------------------
# call LSH python script
# ---------------------------------------------
call_lsh<-function(){
 
  system('python lsh/lsh.py',intern=T)
 
}

# ---------------------------------------------
# call Louvian script
# ---------------------------------------------
call_louvain<-function(){
  
  cat("Assiging Louvain Communities...................... please wait....\n")
  system("rm output.graph",ignore.stderr = T)
  system(paste("sh louvain_script.sh", LOUVAIN_DIR,  "src_dst_lsh.csv", "output.graph",sep = " "))
}

# ---------------------------------------------
# LSH + Louvain partition
# ---------------------------------------------
dropClust_sampling<-function(data_mat, s_ids,true_id){
  
  ranger_preprocess(data_mat, s_ids)
  
  cat("Computing Edgelist using LSH...\n")
  cat(paste("Running LSH on",length(s_ids),"samples. This may take some time...\n"))
  lsh_t = Sys.time()
  call_lsh()
  
  call_louvain()
  
  output <- read.csv("output.graph", sep="",header = F)

  sc_metric(output$V2+1, as.numeric(true_id)[s_ids])
  
  ## READ LOUVEN CLUSTERS FOR SUB-SAMPLING
 
  subsamples_louvain<-sampling()
  
  cat(paste("number of sub-samples:",length(subsamples_louvain),"\n"))
  
  write.csv(x = subsamples_louvain, file = "subsamples_idx",quote = F,row.names =F)
  write.csv(x = data_mat$barcodes[sample_ids[subsamples_louvain]], file = "barcodes_subsamples.csv",quote = F,row.names =F)
  
  table(true_cls_id[sample_ids[subsamples_louvain]])
  
  return(subsamples_louvain)
}


# ---------------------------------------------
# Preprocessesing : similar to pipeline
# described in Zheng et al.
# ---------------------------------------------

ranger_preprocess<-function(data_mat, s_ids){
  
  ngenes_keep = 1000
  write.csv(x = data_mat$gene_symbols, file = "gene_symbols.csv",quote = F,row.names =F)
  #l<-normalize_by_umi(data_mat)   
  l<-normalize_by_umi_2(data_mat)   
  m_n<-l$m
  cat("Select variable Genes...\n")
  df<-get_variable_gene(m_n)
  gc()
  cat("Sort Top Genes...\n")
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  cat("Cutoff Genes...\n")
  df$used<-df$dispersion_norm >= disp_cut_off
  
  features = head(order(-df$dispersion_norm),ngenes_keep)
  system("rm genes",ignore.stderr = T)
  write.csv(features, file = "genes", quote = F,row.names = F)
  write.csv(l$use_genes[features], file = "genes_used_all", quote = F,row.names = F)
  
  genes = read.csv(file = "genes")
  features = genes$x
  
  #Final Data
  m_n_68K<-m_n[,features]
  m_filt<-Matrix(log2(m_n_68K+1),sparse = T)
  cat(paste("Writing Log Normalized whole_matrix, DIM:",dim(m_filt)[1], dim(m_filt)[2],"...\n"))
  system("rm whole_matrix",ignore.stderr = T)
  writeMM(m_filt,file="whole_matrix")
  
  m_n_1000<-m_filt[s_ids,]
  cat(paste("Writing Log Normalized sub_matrix, DIM:",dim(m_n_1000)[1], dim(m_n_1000)[2],"...\n"))
  system("rm sub_matrix",ignore.stderr = T)
  writeMM(m_n_1000,file="sub_matrix")
  
}

# ---------------------------------------------
# External Metrics for clustering evaluation
# ---------------------------------------------
sc_metric<-function(pred_ids,true_ids,show_tab = T){
 dm<-table(pred_ids,true_ids)
 if(show_tab==T){print(dm)}
  r.i. = flexclust::randIndex(dm, correct=TRUE, original=TRUE)
  return(c(r.i.,"Purity" = ClusterPurity(pred_ids,true_ids)))
}

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}


# ---------------------------------------------
# Custom sampling using exponential decay 
# ---------------------------------------------

sampling<-function(pinit=0.195, pfin = 0.9, K=500){
  output <- read.csv("output.graph", sep="",header = F)
  output<-output+1
  oldseed = .Random.seed
  set.seed(0)
  cluster_freq = table(output$V2)#tdf_n_1000$cls_id[sample_ids]
  #prop = (1-exp(-30/cluster_freq))*cluster_freq
  # pinit=0.1
  # pfin = 0.8
  # K=500
  prop = round((pinit - exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
  t(rbind(cluster_freq,prop))
  prop = melt(prop)$value
  subsamples_louvain<-c()
  for(k in 1:length(prop)){
    subsamples_louvain = c(subsamples_louvain,
                           sample(output[which(output$V2==(k)),1],size = prop[k],replace = F))
  }
  .Random.seed = oldseed
  return(subsamples_louvain)
}

# ---------------------------------------------
# pc loading genes selection
# ---------------------------------------------

pc_genes<-function(mat,top=200){
  nc = 50
  PRR <- irlba(as.matrix(mat),nc)
  S = c()
  for(i in 1:nc)
  {
    K = PRR$u[,i]
    L = mclust::Mclust(K)
    S = c(S,L$G)
    
  }
  gene_pca <-abs(PRR$v[,which(S>=3)])
  gene_max<-apply(gene_pca, 1, max) 
  rank_gene <- order(gene_max, decreasing=TRUE)
  
  
  return(head(rank_gene,top))
}


# ---------------------------------------------
# Clustering of sub-samples using 
# Hierarchical Clustering
# ---------------------------------------------

ss_clustering<-function(ss_sel_genes){
  d = dist(ss_sel_genes)
  hc<-fastcluster::hclust(d,method = "average")
  
  hc_labs<-cutreeDynamic(dendro = hc, cutHeight = NULL,
                         minClusterSize = 20,
                         method = "hybrid", deepSplit = 3,
                         pamStage = TRUE,  distM = as.matrix(d), maxPamDist = 0,
                         verbose = 0)
  
  cat(paste("Number of Clusters:", length(unique(hc_labs))-1),"\n")
  
  
  outliers_ids <- which(hc_labs==0)
  subsamples_louv_68K = sample_ids[subsamples_louvain[-outliers_ids]]
  write.csv(x = subsamples_louv_68K, file = "subsamples_louv_68K",quote = F,row.names =F)
  write.csv(x = hc_labs, file = "hc_labs.csv",quote = F,row.names =F)
  
  cat("Sub-sample Metric\n")
  print(sc_metric(hc_labs[-outliers_ids], true_cls_id[subsamples_louv_68K],show_tab = F))
  
  
  hc_labs_clean  = hc_labs[-outliers_ids]
  cat(paste("Predicted Clusters:",length(unique(hc_labs_clean))))
  
  return(list("labels"=hc_labs_clean,"outliers" = outliers_ids))
}

# ---------------------------------------------
# Normalize the gene barcode matrix by umi 
# Filter based on read count first
# ---------------------------------------------

normalize_by_umi_2 <-function(x) {
  mat  = x$mat
  gene_symbols = x$gene_symbols
  cs <- colSums(mat>2)
  x_use_genes <- which(cs > 3)
  
  x_filt<-mat[,x_use_genes]
  gene_symbols = gene_symbols[x_use_genes]
  cat("Dimensions of filtered Matrix:")
  cat(paste(dim(x_filt),"\n"))
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=gene_symbols)
}


# --------------------------------------------------
# Get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
get_variable_gene<-function(m) {
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}


# --------------------------------------------------
# Palette for predicted/true clusters
# --------------------------------------------------

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
getColors<-function(n){
  mycolors = c("#00fe0a", "#ff0000", "#bded1b", "#794b05", "#c3beb7",
               "#0000ff", "#00ffff","#ff21d3" , "#81b7dd","#f87791" ,
               "#1e7309", "#fc9a07", "#625b51", "#6a09c3", "#189ff5",
               "#d19d00", "#0ebf06", "#88ffb3", "#f6fc2a", "#000000")
  if(n>20){
    cat("Too many colors...Using fallback color scheme.\n")
    return(getPalette(n))
  }
  return(mycolors[1:n])
} 

# --------------------------------------------------
# 2D - plots
# --------------------------------------------------

all_plot<-function(plot_proj_df,filename, title){
  x.mean = aggregate(plot_proj_df$Y1, list(plot_proj_df$color), median)[,-1]
  y.mean = aggregate(plot_proj_df$Y2, list(plot_proj_df$color), median)[,-1]
  
  colorcount_t = length(unique(plot_proj_df$color))
  
  pdf(filename,width = 6,height = 5)
  
  p<-ggplot(plot_proj_df,aes(Y1,Y2,col= color))
  p2<-p+ geom_point(size=0.3)  + scale_colour_manual(values =  getColors(colorcount_t))+
    ggtitle(title)+
    annotate("text", x = x.mean, y = y.mean, label = levels(plot_proj_df$color), size =2.75 )+
    guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))+ylab("t-SNE 2")+xlab("t-SNE 1")+
    theme_classic()
  
  print(p2)
  
  dev.off()
}
