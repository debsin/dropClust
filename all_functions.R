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
  
  cat("Assiging Louvain Communities......................\n")
  suppressWarnings(file.remove("output.graph"))
  
  if(Sys.info()["sysname"]=="Windows"){
    shell(paste0(file.path(LOUVAIN_DIR,"win32"),"/convert -i src_dst_lsh.csv -o  graph.bin"))
    shell(paste0(file.path(LOUVAIN_DIR,"win32"),"/louvain graph.bin -l -1 -q 0 > graph.tree"))
    cc = shell(paste0(file.path(LOUVAIN_DIR,"win32"),"/hierarchy graph.tree"),intern = T)
    tlevels = as.numeric(gsub("[^\\d]+", "", cc[1], perl=TRUE)) - 1
    shell(paste0(file.path(LOUVAIN_DIR,"win32"),"/hierarchy -l ", tlevels," graph.tree > output.graph"))
  } else {
    system(paste0(LOUVAIN_DIR,"convert -i src_dst_lsh.csv -o  graph.bin"))
    system(paste0(LOUVAIN_DIR,"louvain graph.bin -l -1 -q 0 > graph.tree"))
    cc = system(paste0(LOUVAIN_DIR, "hierarchy graph.tree"),intern = T)
    tlevels = max(as.numeric(gsub("[^\\d]+", "", cc[1], perl=TRUE)) - 1,0)
    system(paste0(LOUVAIN_DIR,"hierarchy -l ", tlevels," graph.tree > output.graph"))
  }
  t <- file.remove(paste0("graph.tree"))
  t <- file.remove(paste0("graph.bin"))
  cat("Done.\n")
}

# ---------------------------------------------
# LSH + Louvain partition
# ---------------------------------------------
dropClust_sampling<-function(data_mat, s_ids,true_id){
  
  #ranger_preprocess(data_mat, s_ids)
  
  print("Computing Edgelist using LSH.")
  print(paste("Running LSH on",length(s_ids),"samples. This may take some time..."))
  lsh_t = Sys.time()
  call_lsh()
  
  call_louvain()
  
  output <- read.csv("output.graph", sep="",header = F)

  #sc_metric(output$V2+1, as.numeric(true_id)[s_ids])
  
  ## READ LOUVEN CLUSTERS FOR SUB-SAMPLING
 
  subsamples_louvain<-sampling()
  
  print(paste("number of sub-samples:",length(subsamples_louvain)))
  
  write.csv(x = subsamples_louvain, file = "subsamples_idx",quote = F,row.names =F)
  write.csv(x = data_mat$barcodes[s_ids[subsamples_louvain]], file = "barcodes_subsamples.csv",quote = F,row.names =F)
  
  table(true_id[s_ids[subsamples_louvain]])
  
  return(subsamples_louvain)
}


# ---------------------------------------------
# Dispersion Genes and Matrix Subsetting
# ---------------------------------------------
matrix.subset<-function(normalized_data, s_ids, ngenes_keep = 1000){
  
  # write.csv(x = data_mat$gene_symbols, file = "gene_symbols.csv",quote = F,row.names =F)
  # #l<-normalize_by_umi(data_mat)   
  # l<-normalize_by_umi_2(data_mat)   
  # m_n<-l$m
  print("Select variable Genes...")
  df<-get_dispersion_genes(normalized_data$m)
  gc()
  print("Sort Top Genes...")
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  print("Cutoff Genes...")
  df$used<-df$dispersion_norm >= disp_cut_off
  
  features = head(order(-df$dispersion_norm),ngenes_keep)
  system("rm genes",ignore.stderr = T)
  write.csv(features, file = "genes", quote = F,row.names = F)
  write.csv(normalized_data$use_genes[features], file = "genes_used_all", quote = F,row.names = F)
  
  genes = read.csv(file = "genes")
  features = genes$x
  
  #Final Data
  m_n_whole<-normalized_data$m[,features]
  m_filt<-Matrix(log2(m_n_whole+1),sparse = T)
  print(paste("Writing Log Normalized whole_matrix, DIM:",dim(m_filt)[1], dim(m_filt)[2],"..."))
  system("rm whole_matrix",ignore.stderr = T)
  writeMM(m_filt,file="whole_matrix")
  
  m_n_ngenes<-m_filt[s_ids,]
  print(paste("Writing Log Normalized sub_matrix, DIM:",dim(m_n_ngenes)[1], dim(m_n_ngenes)[2],"..."))
  system("rm sub_matrix",ignore.stderr = T)
  writeMM(m_n_ngenes,file="sub_matrix")
  
  return(m_n_ngenes)
  
}



# ---------------------------------------------
# Normalize the gene barcode matrix by umi 
# Filter based on read count first
# ---------------------------------------------
normalize_by_umi_2 <-function(x, min.count=2, min.cell=3) {
  mat  = x$mat
  gene_symbols = x$gene_symbols
  cs <- Matrix::colSums(mat>min.count)
  x_use_genes <- which(cs > min.cell)
  
  x_filt<-mat[,x_use_genes]
  gene_symbols = gene_symbols[x_use_genes]
  print("Dimensions of filtered Matrix:")
  print(paste(dim(x_filt),""))
  rs<-Matrix::rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=gene_symbols)
}


# --------------------------------------------------
# Get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
get_dispersion_genes<-function(mat) {
  
  st<-data.frame(avg=Matrix::colMeans(mat),cov=apply(mat,2,sd)/Matrix::colMeans(mat),variance=apply(mat,2,var))
  st$dispersion<-with(st,variance/avg)
  st$avg_bin<-with(st,cut(avg,breaks=c(-Inf,quantile(avg,seq(0.1,1,0.05)),Inf)))
  variance_by_bin<-ddply(st,"avg_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  st$bin_disp_median<-variance_by_bin$bin_median[match(st$avg_bin,variance_by_bin$avg_bin)]
  st$bin_disp_mad<-variance_by_bin$bin_mad[match(st$avg_bin,variance_by_bin$avg_bin)]
  st$dispersion_norm<-with(st,abs(dispersion-bin_disp_median)/bin_disp_mad)
  st
}

# ---------------------------------------------
# Filter Cells
# ---------------------------------------------
filter_cells<-function(data,th){
  #Return good cells based on IQR rule
  rs<-Matrix::rowSums(data$mat)
  keep_cells = which(rs>=th)
  #q<-quantile(rs)
  #iqr<-IQR(rs)
  
  #l1 = q[2] - 1.5*(iqr)
  #l2 = q[4] + 1.5*(iqr)
  
  #keep_cells = intersect(which(rs > l1), which(rs <l2))
  print(paste(length(rs)-length(keep_cells), "bad cells present."))
  data$mat = data$mat[keep_cells,]
  data$barcodes = data$barcodes[keep_cells]
  keep_cells = list(keep_cells)
  names(keep_cells)="keep_cells"
  data = append(data, keep_cells)
  return(data)
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
  cluster_freq = table(output$V2)
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
    L = mclust::Mclust(K,verbose = F)
    S = c(S,L$G)
    
  }
  gene_pca <-abs(PRR$v[,which(S>=3)])
  gene_max<-apply(gene_pca, 1, max) 
  rank_gene <- order(gene_max, decreasing=TRUE)
  
  df=data.frame(cbind("ids" = c(1:nc), "value"=as.factor(S)))
  p<-ggplot(df, aes (ids, value))+geom_bar(stat="identity",  fill="steelblue")+
    xlab(paste("Top", nc,"Principal Components"))+ylab("# Gaussian Component")+ggtitle("")+theme_minimal()
  print(p)
  return(head(rank_gene,top))
}


# ---------------------------------------------
# Clustering of sub-samples using 
# Hierarchical Clustering
# ---------------------------------------------

ss_clustering<-function(ss_sel_genes, minClusterSize = 20, deepSplit = 3){
  d = dist(ss_sel_genes)
  hc<-fastcluster::hclust(d,method = "average")
  
  hc_labs<-cutreeDynamic(dendro = hc, cutHeight = NULL,
                         minClusterSize = minClusterSize,
                         method = "hybrid", deepSplit = deepSplit,
                         pamStage = TRUE,  distM = as.matrix(d), maxPamDist = 0,
                         verbose = 0)
  
  
  
  outliers_ids <- which(hc_labs==0)
  subsamples_louv_68K = sample_ids[subsamples_louvain[-outliers_ids]]
  write.csv(x = subsamples_louv_68K, file = "subsamples_louv_68K",quote = F,row.names =F)
  write.csv(x = hc_labs, file = "hc_labs.csv",quote = F,row.names =F)
  
  hc_labs_clean  = hc_labs[-outliers_ids]
  print(paste("Predicted Clusters:",length(unique(hc_labs_clean))))
  
  return(list("labels"=hc_labs_clean,"outliers" = outliers_ids))
}

# ---------------------------------------------
# Read 10X data
# ---------------------------------------------
read10X<-function(data.dir=NULL){
  full_data <- list()
  for(i in seq_along(data.dir)){
    run <- data.dir[i]
    if (!dir.exists(run)){
      stop("Directory provided does not exist")
    }
    
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    
    barcode.loc <- paste(run, "barcodes.tsv", sep ="")
    gene.loc <- paste(run, "genes.tsv", sep ="")
    matrix.loc <- paste(run, "matrix.mtx", sep ="")
    
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (!file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    
    data <- readMM(matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if(all(grepl("\\-1$", cell.names)) == TRUE) {
      cell.names <- as.vector(as.character(sapply(cell.names, extract_field, 1, delim = "-")))
    }
    rownames(data) <- make.unique(as.character(sapply(gene.names, extract_field, 2, delim = "\\t"))) 
    
    if(is.null(names(data.dir))){
      if(i < 2){
        colnames(data) <- cell.names
      }
      else {
        colnames(data) <- paste0(i, "_", cell.names, sep = "") 
      }
    } else {
      colnames(data) <- paste0(names(data.dir)[i],"_",cell.names) 
    }
    full_data <- append(full_data, data)
  }
  full_data <- do.call(cbind, full_data)
  full_data <- list("mat" = t(data), "barcodes" = colnames(data) , "gene_symbols" = rownames(data))
  return(full_data)
  
}

extract_field=function(string,field=1,delim="_") {
  fields=as.numeric(unlist(strsplit(as.character(field),",")))
  if (length(fields)==1)  return(strsplit(string,delim)[[1]][field])
  return(paste(strsplit(string,delim)[[1]][fields],collapse = delim))
}


# ---------------------------------------------
# Optimize parameter 'pinit' for 
# exponential decay sampling
# ---------------------------------------------
optimized_Pinit<-function(nsamples = 500, K=500, pfin = 0.9){
  oldseed = .Random.seed
  set.seed(1234) # The user can use any seed.
  dimension <- 1
  global.min <- 0
  tol <- 1e-4
  lower <- rep(0.01, dimension)
  upper <- rep(1.0, dimension)
  
  
  pin_find<-function(p = par,MAX_c = nsamples){
    output <- read.csv("output.graph", sep="",header = F)
    output<-output+1
    oldseed = .Random.seed
    set.seed(0)
    cluster_freq = table(output$V2)
    pinit = p[1]
    #pfin = p[2]
    prop = round((pinit - exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
    t(rbind(cluster_freq,prop))
    prop = melt(prop)$value
    subsamples_louvain<-c()
    for(k in 1:length(prop)){
      subsamples_louvain = c(subsamples_louvain,
                             sample(output[which(output$V2==(k)),1],size = prop[k],replace = F))
    }
    .Random.seed = oldseed
    
    return(abs(MAX_c-length(subsamples_louvain)))
  }
  
  
  out <- GenSA(fn = pin_find,lower = lower,upper = upper,
               control=list(threshold.stop=global.min+tol,verbose=TRUE))
  
  .Random.seed = oldseed
  
  #out[c("value","par","counts")]
  
  return(out$par)
}


# --------------------------------------------------
# Cluster Assignment
# --------------------------------------------------

cluster_assign<-function(INDEX, ss_clusters){
  INDEX = as.matrix(INDEX)+1
  
  clust_col = list()
  
  for( i in 1:nrow(INDEX))
  {
    clust_col  = list(clust_col ,names(which.max(table(ss_clusters$labels[INDEX[i,]])>1)))
  }
  clust_col = as.numeric(unlist(clust_col))
  return (clust_col)
}

# --------------------------------------------------
# @d Projection
# --------------------------------------------------

compute_2d_embedding<-function(data, ss_clusters, INDEX){
  INDEX = as.matrix(INDEX)+1
  set.seed(0)
  # print("Compute Tsne Projection using PCA top genes...")
  ts<-Rtsne(data,perplexity = 20,dims = 2)
  
  
  # print("Projecting TSNE co-ordinates...")
  tt<-ts$Y[-ss_clusters$outliers,]
  PROJ = matrix(NA, nrow=nrow(INDEX), ncol=2)
  
  for( i in 1:nrow(INDEX))
  {
    maj_id = clust_col[i]
    maj_col = which(ss_clusters$labels[INDEX[i,]]==maj_id)
    
    if(length(maj_col)==1){
      PROJ[i,] = Matrix::colMeans(tt[INDEX[i,],])
    }
    else{
      PROJ[i,] = Matrix::colMeans(tt[INDEX[i,maj_col],])
    }
    
  }
  
  dim(PROJ)
  return(PROJ)
  
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
    print("Too many colors...Using fallback color scheme.")
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
  
  jpeg(file.path(FIG_DIR,filename),width = 450,height = 450,res = 150,quality = 90)
  
  p<-ggplot(plot_proj_df,aes(Y1,Y2,col= color))
  p2<-p+ geom_point(size=0.3)  + scale_colour_manual(values =  getColors(colorcount_t))+
    ggtitle(title)+
    annotate("text", x = x.mean, y = y.mean, label = levels(plot_proj_df$color), size =2.75 )+
    guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))+ylab("t-SNE 2")+xlab("t-SNE 1")+
    theme_classic()
  
  print(p2)
  
  dev.off()
}

source("NODES.R")
# ---------------------------------------------
# Sample fixed number of random 
# cells from clusters
# ---------------------------------------------
sample_from_cluster<-function(pred_labels, size){
  set.seed(0)
  fixed_samples = c()
  
  
  for(clust_id in levels(as.factor(pred_labels$x))){
    if(clust_id==0) next;
    sample_ids_per_cluster = which(pred_labels$x==clust_id)
    fixed_samples = c(fixed_samples, sample(sample_ids_per_cluster,min(size,length(sample_ids_per_cluster))))
  }
  return(fixed_samples)
}

# ---------------------------------------------
# Score DE genes for each pair of clusters
# ---------------------------------------------
DE_genes <- function(raw_data,labels,max=0,lfc_th,q_th, min.count=3, min.cell.per=0.5 )
{
  #dim(data)
  
  #raw_data = 2^data
  #min(raw_data)
  
  min_cell = floor(ncol(raw_data)*(min.cell.per/100))
  keep_genes = apply(raw_data, 1 , function(x)  sum(x>min.count))>=min_cell
  
  #table(omit_genes)
  
  raw_data = raw_data[keep_genes,]
  
  # call NODES for each class vs class
  DE_list = list()
  pairs = as.matrix(combn(unique(labels),2))
  pair_ids = paste(pairs[1,],pairs[2,],sep="_")
  
  DE_list<-foreach(i = 1:dim(pairs)[2]) %dopar%
  {
    IND_a = which(labels == pairs[1,i])
    #IND_b = setdiff(1:ncol(data),IND_a)
    IND_b = which(labels == pairs[2,i])
    
    
    DE_res = NODES1(raw_data[,c(IND_a,IND_b)],group =c(rep("A",length(IND_a)),rep("B",length(IND_b))))
    LFC = log2(Matrix::rowMeans(raw_data[,IND_a])/Matrix::rowMeans(raw_data[,IND_b]))
    
    
    #length(names(LFC)==rownames(DE_res))
    #head(DE_res)
    #rN = rownames(DE_res)[head(intersect(which(abs(LFC)>=1), which(DE_res$qvalues<=0.05)),30)]
    #sig = DE_res[which(DE_res$qvalues<=0.05),]
    sig = DE_res[intersect(which(abs(LFC)>=lfc_th), which(DE_res$qvalues<=q_th)),]
    if(max>0){
      rN = head(rownames(sig)[order(sig$qvalues)],max)
    }
    else{
      rN = rownames(sig)[order(sig$qvalues)]
    }
    
    #DE_up = union(DE_up,rN)
    DE_list[[pair_ids[i]]] = data.frame(gene = rN,q_val = DE_res$qvalues[match(rN,rownames(DE_res))],fc = LFC[match(rN,names(LFC))])
    #return(data.frame(gene = rN,q_val = DE_res$qvalues[match(rN,rownames(DE_res))],fc = LFC[match(rN,names(LFC))]))
    
  }
  
  DE_up= c()
  for(i in 1:length(DE_list)){
    DE_up = union(DE_up, DE_list[[i]]$gene)
  }
  gc()
  names(DE_list) = pair_ids
  RES = list(genes = DE_up,DE_res = DE_list)
  return(RES)
}

# ---------------------------------------------
# List cell type specific genes
# ---------------------------------------------
find_ct_genes<-function(ID, DE_genes_nodes_all, Mat_ct){
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
        ## indices of groups
        id_1 = which(ID == i)
        id_2 = which(ID == j)
        pair_id1 = paste0(i,"_",j)
        pair_id2 = paste0(j,"_",i)
        #pair_id = union()
        ## comparing the avg expression - or the direction
        diff_i_j = union(DE_genes_nodes_all[["DE_res"]][[pair_id2]]$gene,DE_genes_nodes_all[["DE_res"]][[pair_id1]]$gene)
        
        RM_i = Matrix::rowMeans(Mat_ct[diff_i_j,id_1])
        RM_j = Matrix::rowMeans(Mat_ct[diff_i_j,id_2])
        
        ## direction
        g = rownames(Mat_ct[diff_i_j,])[which(RM_i > RM_j)]
        ct_genes = append(ct_genes, g)
      }
    }
    
    #ct_genes_list[[paste0(i)]] = c(ct_genes_list[[paste0(i)]], ct_genes)
    print(paste(pair_id1,pair_id2,length(names(table(ct_genes))[which(table(ct_genes) == (length(unique(ID))-1))])))
    all_ct_genes = append(all_ct_genes, names(table(ct_genes))[which(table(ct_genes) == (length(unique(ID))-1))])
    ct_genes_list[[paste0(i)]] = c(ct_genes_list[[paste0(i)]], names(table(ct_genes))[which(table(ct_genes) == (length(unique(ID))-1))])
  }
  return(list("all_ct_genes" = all_ct_genes, "ct_genes_list"=ct_genes_list))
}

# ---------------------------------------------
# Plot heatmap of cell-type specific genes
# for the samples from each clusters
# ---------------------------------------------
plot_heat_genes <- function(data,label,filename)
{
  
  colors = getColors(length(unique(label)))
  colors = colors[label]
  ordered_freq = table(factor(label,levels=unique(label)))
  a = cumsum(ordered_freq)
  b = c(0,  a[-length(a)])
  pos = round((a+b)/2)
  text_lab = rep(NA,dim(heat_in)[2])
  text_lab[pos] = unique(label);
  
  
  myPalette <-colorRampPalette(brewer.pal(11,"RdBu"))
  myPalette <- colorRampPalette(c("grey15","gray10","black","yellow","darkred"))
  
  
  
  pdf(filename,width=8, height=16)
  heatmap.2(data, 
            trace="none", col=myPalette(40),Colv = FALSE, Rowv = FALSE, srtCol=0, 
            ColSideColors=colors, cexCol = 1,cexRow = 0.2,labCol = text_lab,
            dendrogram="none",scale="row",density.info="none",key=FALSE, 
            main="DE genes and equal #samples")
  
  dev.off()
  
  heatmap.2(data, 
            trace="none", col=myPalette(40),Colv = FALSE, Rowv = FALSE, srtCol=0, 
            ColSideColors=colors, cexCol = 1,cexRow = 0.2,labCol = text_lab,
            dendrogram="none",scale="row",density.info="none",key=FALSE, 
            main="DE genes and equal #samples")
  
}

# ---------------------------------------------
# Write celltype specific genes
# ---------------------------------------------
write_ct_genes<-function(data, filename){
  n.obs <- sapply(data, length)
  seq.max <- seq_len(max(n.obs))
  mat <- as.data.frame(t(sapply(data, "[", i = seq.max)))
  gene.df <- sapply(mat, as.character)
  gene.df[is.na(gene.df)] <- ""
  write.csv(gene.df, file = filename,quote = F)
}

# ---------------------------------------------
# Wtite DE genes scores for each 
# pair of clusters
# ---------------------------------------------
write_de_pairs<-function(filename, data){
  sink(filename)
  print(data)
  sink()
}

