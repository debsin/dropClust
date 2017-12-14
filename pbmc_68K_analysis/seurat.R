#!/usr/bin/Rscript
#!/usr/bin/env Rscript
library("methods")
library(Seurat)
library(dplyr)
library(Matrix)


setwd("~/Projects/seurat/")
# Load the PBMC dataset

pbmc.data<-Read10X("~/Projects/10Xpaper/filtered_matrices_mex/hg19/")

#pbmc.data = pbmc_68k
#Examine the memory savings between regular and sparse matrices
# dense.size <- object.size(as.matrix(pbmc.data))
# dense.size
# 
# sparse.size <- object.size(pbmc.data)
# sparse.size
# dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data)
# Note that this is slightly different than the older Seurat workflow, where log-normalized values were passed in directly.
# You can continue to pass in log-normalized values, just set do.logNormalize=F in the next step.
pbmc <- new("seurat", raw.data = pbmc.data)

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
pbmc <- Setup(pbmc, min.cells = 3, min.genes = 1, do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC")


#nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
#VlnPlot(pbmc, c("nGene", "nUMI", "percent.mito"), nCol = 3)

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well
# par(mfrow = c(1, 2))
# GenePlot(pbmc, "nUMI", "percent.mito")
# GenePlot(pbmc, "nUMI", "nGene")

#We filter out cells that have unique gene counts over 2,500 and under 500, and > 5% mitochondrial percentage
# pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
# pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)
# pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.low = 500)

pbmc <- MeanVarPlot(pbmc, x.low.cutoff = 0, y.cutoff = 0.8)

length(pbmc@var.genes)


print("VarPlot completed")

#note that this overwrites pbmc@scale.data. Therefore, if you intend to use RegressOut, you can set do.scale=F and do.center=F in the original object to save some time.
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"), genes.regress = pbmc@var.genes)


print("RegressOut completed")




#pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
pbmc<- PCAFast(pbmc, pc.genes = pbmc@var.genes, pcs.compute = 40, pcs.print = 30)

pbmc<- RunTSNE(pbmc, dims.use = 1:25, do.fast = T)


print("Tsne completed")


#save.SNN=T saves the SNN so that the  SLM algorithm can be rerun using the same graph, but with a different resolution value (see docs for full details)
#pbmc <- FindClusters(pbmc ,pc.use = 1:25, resolution = seq(2,4,0.5), save.SNN = T, do.sparse = T)
pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = T, do.sparse = T)


print("SNN completed")

pdf("Seurat.pdf")
TSNEPlot(pbmc, do.label = T)
dev.off()
save(pbmc, file = "seurat_tutorial.Robj")

#write.csv(as.data.frame(pbmc@ident),"seurat_predicted.csv",quote=F)
write.csv(as.data.frame(cbind(pbmc@ident,pbmc@tsne.rot)),"seurat_predicted.csv", quote=F)



