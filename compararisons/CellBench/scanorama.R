library(reticulate)


common_genes = intersect(rownames(objects[[1]]),
                         intersect(rownames(objects[[2]]), rownames(objects[[3]])))

datasets <- list()
genes_list <- list()


datasets[[1]] = Matrix::t(counts(objects[[1]])[common_genes,])
datasets[[2]] = Matrix::t(counts(objects[[2]])[common_genes,])
datasets[[3]] = Matrix::t(counts(objects[[3]])[common_genes,])



genes_list[[1]] = common_genes
genes_list[[2]] = common_genes
genes_list[[3]] = common_genes

print("Done loading")


scanorama <- reticulate::import('scanorama')
scipy <- reticulate::import('scipy')


# Integration and batch correction.
integrated.corrected.data <- scanorama$correct(datasets, genes_list = genes_list,  return_dimred = T, return_dense = T)


dimred = integrated.corrected.data[[1]]

scrama_corr  = rbind(dimred[[1]], dimred[[2]], dimred[[3]])

