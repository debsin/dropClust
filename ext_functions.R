
#Read 10X data
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
  return(full_data)
  
}


# call LSH python script

call_lsh<-function(){
 
  system('python lsh/lsh.py',intern=T)
 
}

call_louvain<-function(){
  
  cat("Assiging Louvain Communities...................... please wait....\n")
  system("rm output.graph",ignore.stderr = T)
  system(paste("./louvain_script.sh", LOUVAIN_DIR,  "src_dst_lsh.csv", "output.graph",sep = " "))
  #output <- read.csv("output.graph", sep="",header = F)
  
  #./convert -i ~/Projects/scClust/git/scClust/temp_graph.csv -o ~/Projects/scClust/git/scClust/graph.bin
  #./community ~/Projects/scClust/git/scClust/graph.bin -l -1 -q 0.0000001 -v > graph.tree
  #./hierarchy graph.tree -l 1 > ~/Projects/scClust/git/scClust/output.graph
  #return(output)
}


method_scClust<-function(data_mat, s_ids,true_id){
  
  ranger_preprocess(data_mat, s_ids)
  
  cat("Computing Edgelist using LSH...\n")
  cat(paste("Running LSH on",length(s_ids),"samples. This may take some time...\n"))
  lsh_t = Sys.time()
  call_lsh()
  
  
  call_louvain()
  cat(paste("LSH + Louvain Time... ", difftime(Sys.time(),lsh_t, units = "mins"),"\n"))
  output <- read.csv("output.graph", sep="",header = F)

  return(sc_metric(output$V2+1, as.numeric(true_id)[s_ids]))
}


ranger_preprocess<-function(data_mat, s_ids){
  
  ngenes_keep = 1000
  #l<-.normalize_by_umi(data_mat)   
  l<-.normalize_by_umi_2(data_mat)   
  m_n<-l$m
  cat("Select variable Genes...\n")
  df<-.get_variable_gene(m_n)
  gc()
  cat("Sort Top Genes...\n")
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  cat("Cutoff Genes...\n")
  df$used<-df$dispersion_norm >= disp_cut_off
  
  features = head(order(-df$dispersion_norm),ngenes_keep)
  system("rm genes")
  write.csv(features, file = "genes", quote = F,row.names = F)
  write.csv(l$use_genes[features], file = "genes_used_all", quote = F,row.names = F)
  
  genes = read.csv(file = "genes")
  features = genes$x
  
  #Final Data
  m_n_68K<-m_n[,features]
  m_filt<-Matrix(log2(m_n_68K+1),sparse = T)
  cat(paste("Writing Log Normalized whole_matrix, DIM:",dim(m_filt)[1], dim(m_filt)[2],"...\n"))
  system("rm whole_matrix")
  writeMM(m_filt,file="whole_matrix")
  
  cat("Subset Genes for Matirx...\n")
  m_n_1000<-m_filt[s_ids,]
  cat(paste("Writing Log Normalized sub_matrix, DIM:",dim(m_n_1000)[1], dim(m_n_1000)[2],"...\n"))
  system("rm sub_matrix")
  writeMM(m_n_1000,file="sub_matrix")
  
}

sc_metric<-function(pred_ids,true_ids,show_tab = T){
 dm<-table(pred_ids,true_ids)
 if(show_tab==T){print(dm)}
  r.i. = flexclust::randIndex(dm, correct=TRUE, original=TRUE)
  return(c(r.i.,"Purity" = ClusterPurity(pred_ids,true_ids)))
}




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
  
  write.csv(x = rank_gene[1:top], file = "subgene_idx",quote = F,row.names =F)
  write.csv(x = genes2K$x[rank_gene[1:top]], file = "subgene_idx_68K",quote = F,row.names =F)
  
  return(rank_gene[1:top])
}

# ---------------------------------------------
# normalize the gene barcode matrix by umi 
# filter based on read count first
# ---------------------------------------------

.normalize_by_umi_2 <-function(x) {
  mat  = x$mat
  gene_symbols = x$gene_symbols
  cs <- colSums(mat>2)
  x_use_genes <- which(cs > 3)
  
  x_filt<-mat[,x_use_genes]
  gene_symbols = gene_symbols[x_use_genes]
  print(dim(x_filt))
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=gene_symbols)
}



# ---------------------------------------------
# normalize the gene barcode matrix by umi
# ---------------------------------------------
.normalize_by_umi <-function(x) {
  cs <- colSums(x)
  x_use_genes <- which(cs >= 1)
  x_filt<-x[,x_use_genes]
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=x_use_genes)
}

# --------------------------------------------------
# get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
.get_variable_gene<-function(m) {
  
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

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

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
