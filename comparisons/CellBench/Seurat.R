# install.packages('devtools')
# # Replace '2.3.0' with your desired version
# devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))

library(Seurat)

common_genes = intersect(rownames(objects[[1]]),
                         intersect(rownames(objects[[2]]), rownames(objects[[3]])))
data.tables <- list()

data.tables[[1]] = counts(objects[[1]])[common_genes,]
colnames(data.tables[[1]]) = paste(colnames(objects[[1]]), 1, sep = "_")

data.tables[[2]] = counts(objects[[2]])[common_genes,]
colnames(data.tables[[2]]) = paste(colnames(objects[[2]]), 3, sep = "_")

data.tables[[3]] = counts(objects[[3]])[common_genes,]
colnames(data.tables[[3]]) = paste(colnames(objects[[3]]), 3, sep = "_")




names(data.tables) = paste("Batch",c(1,2,3))


ob.list <- list()
for (i in 1:length(names(data.tables))) {
  ob.list[[i]] <- CreateSeuratObject(raw.data  = data.tables[[i]])
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableGenes(ob.list[[i]], do.plot = T, display.progress = F)
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]]@meta.data$batch <- toString(i)
}


genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

print("Done intersecting genes")

integrated <- Seurat::RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 20)

print("CCA Done")

# CC Selection
#MetageneBicorPlot(integrated, grouping.var = "tech", dims.eval = 1:15)

# Run rare non-overlapping filtering
#integrated <- CalcVarExpRatio(object = integrated, reduction.type = "pca",
#                                       grouping.var = "tech", dims.use = 1:15)
#integrated <- SubsetData(integrated, subset.name = "var.ratio.pca",
#                         accept.low = 0.5)

# Alignment
integrated <- AlignSubspace(integrated,
                            reduction.type = "cca",
                            grouping.var = "batch",
                            dims.align = 1:15)


seurat_corr = integrated@dr$cca.aligned@cell.embeddings
