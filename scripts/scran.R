library(ggplot2)
library(ggrepel)
library(Matrix)
library(Rtsne)
library(scater)
library(scran)
library(cellrangerRkit)
library(Rtsne)
library(stringr)
library(igraph)
source("SomeFuncs/convenience_functions.R")

################################################################################################
################################################################################################
# download the 10X data
pbmc_path <- "Droplet/10X_data/PBMC/"

# this has to be an absolute path, relative paths break
download_sample(sample_name="pbmc4k", sample_dir=pbmc_path,
                host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")
pbmc.10x <- load_cellranger_matrix(pbmc_path)

# filter non-zero genes and normalize PBMC 10X UMI counts
pbmc.nz <- get_nonzero_genes(pbmc.10x)

################################################################################################
################################################################################################
# normalize using size factors
pbmc.norm <- size_factor_normalize(exprs(pbmc.10x), cell.sparse=0.95,
                                   gene.sparse=0.99, cluster.size=50)

# map ensembl IDs to gene symbols to select a couple of marker genes for visualisation
all.genes <- rownames(pbmc.norm)

ensembl <- useEnsembl(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', GRCh=37)

gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=all.genes)
################################################################################################
################################################################################################
# select highly variable genes (~1-2k will be more than sufficient I think)
pbmc.hvg <- find_hvg(pbmc.norm[, 1:(dim(pbmc.norm)[2]-1)], plot=FALSE)

# the mean-CV^2 plot looks OK, even with the high number of unit counts
# there are 1192 HVG at p<=0.01

write.table(pbmc.norm,
            file="Droplet/10X_data/PBMC/PBMC_norm.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

hvg.df <- cbind.data.frame(names(pbmc.hvg)[pbmc.hvg])
colnames(hvg.df) <- c("HVG")
write.table(hvg.df,
            file="Droplet/10X_data/PBMC/PBMC_hvg.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

################################################################################################
################################################################################################
# tSNE embedding for visualization
pbmc.select <- pbmc.norm[pbmc.hvg, ]
pbmc.select[is.na(pbmc.select)] <- 0

set.seed(42)
pbmc.map <- tsne_wrapper(pbmc.select[, 1:(dim(pbmc.norm)[2]-1)], perplexity=100)

################################################################################################
################################################################################################
# pull out marker genes for PBMC populations to add identities to different clusters
# also try SNNgraph for clustering!  Use just 20 dimensions from the PCA to build the graph
pbmc.snn <- buildSNNGraph(as.matrix(pbmc.norm[, 1:(dim(pbmc.norm)[2]-1)]),
                          k=5, d=20)

# use the walktrap algorithm to form communities
# try a few different values for the step parameter
pbmc.community <- cluster_walktrap(pbmc.snn, steps=4)
n.comms <- length(pbmc.community)

# vertex community membership
pbmc.members <- pbmc.community$membership
pbmc.clusters <- cbind.data.frame(colnames(pbmc.norm[, 1:(dim(pbmc.norm)[2]-1)]), pbmc.members)
colnames(pbmc.clusters) <- c("Sample", "Community")
pbmc.clusters$Community <- as.factor(pbmc.clusters$Community)

pbmc.max <- merge(pbmc.map, pbmc.clusters, by='Sample')

comm.tsne <- ggplot(pbmc.max,
                    aes(x=Dim1, y=Dim2, colour=Community)) +
  geom_point(size=3) + theme_classic()
