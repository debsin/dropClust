library(reticulate)

datasets <- list()
genes_list <- list()


batch_2_ids = grep('B2',labels)
batch_1_ids = setdiff(1:length(labels),grep('B2',labels))


datasets[[1]] = t(as.matrix(rca.data[,batch_1_ids]))
datasets[[2]] = t(as.matrix(rca.data[,batch_2_ids]))
genes_list[[1]] = gene.symbols
genes_list[[2]] = gene.symbols
print("Done loading")


scanorama <- reticulate::import('scanorama')
scipy <- reticulate::import('scipy')


# Integration and batch correction.
integrated.corrected.data <- scanorama$correct(datasets, genes_list = genes_list,  return_dimred = T, return_dense = T)


dimred = integrated.corrected.data[[1]]

scrama_corr  = rbind(dimred[[1]], dimred[[2]])

