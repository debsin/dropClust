
temp<-FilterGenes(mix.data,min.count = 1,min.cell = 1)
temp<- CountNormalize(temp)

dc_merged_mat = LogNormalize(normcounts(temp[rowData(temp)$CommonDEGenes,]))

close_threshold = 0.1
cells_threshold = 0.2

mat = t(dc_merged_mat)

print(dim(mat))
print(close_threshold)
print(cells_threshold)
ncells = nrow(mat)
ngenes = ncol(mat)

vote_matrix = mat.or.vec(nr = ngenes,nc = ngenes)

for(i in 1:(ngenes-1)){
  for(j in (i+1):ngenes){
    delta_ratios = abs(mat[,i]-mat[,j])/pmax(mat[,i],mat[,j])
    delta_ratios[!complete.cases(delta_ratios)] = 0
    if(sum(delta_ratios<close_threshold) > (cells_threshold * ncells)) {
      vote_matrix[i,j] = vote_matrix[i,j] - 1
      vote_matrix[j,i] = vote_matrix[j,i] - 1
    }
    else{

      vote_matrix[i,j] = vote_matrix[i,j] + 1
      vote_matrix[j,i] = vote_matrix[j,i] + 1
    }

  }
}

rs = rowSums(vote_matrix)
sort(rs,decreasing = T)
qualified = order(rs,decreasing = T)
keep_these = head(qualified, 30)


rank_mat <- mat[, keep_these]
### Converting the given matrix to rank vectors
rank_mat <- .get_rank_mat_r(rank_mat)

umap_proj_dc = uwot::umap(as.matrix(rank_mat), metric = 'cosine', n_components = 10, spread = 2)

source("~/GitHub/dropclust_benchmarks/analysis/bc/Cellbench/plots.R")

PROJ_c_dc = umap_proj_dc
plot_proj_df_true_dc<-data.frame(Y1 = PROJ_c_dc[,1],Y2 = PROJ_c_dc[,2], color = as.factor(cc$cell_line_demuxlet), batch = as.integer(cc$Batch))
batch_plot(plot_proj_df_true_dc,filename = NA, title = "Corrected DC ",  type=NULL)

