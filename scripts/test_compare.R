

WORK_DIR = "C:/Projects/dropClust/"
DATA_DIR <- file.path(WORK_DIR,"data")
FIG_DIR <-  file.path(WORK_DIR,"plots")
REPORT_DIR  <- file.path(WORK_DIR,"report")

dir.create(file.path(FIG_DIR),showWarnings = FALSE)
dir.create(file.path(REPORT_DIR),showWarnings = FALSE)

set.seed(0)

# Load Data from path: C:/Projects/dropClust/data/pbmc3k/hg19
pbmc.data<-read_data(file.path(DATA_DIR,'pbmc3k/hg19'), format = "10X")

sce <-readfiles(path = "C:/Projects/data/pbmc3K/hg19/")

all.equal(as.matrix(counts(sce)), as.matrix(Matrix::t(pbmc.data$mat)))


# Filter poor quality cells.  A threshold th corresponds to the total count of a cell.
filtered.data <- filter_cells(pbmc.data)
sce.1<-FilterCells(sce)
all.equal(as.matrix(counts(sce.1)), as.matrix(Matrix::t(filtered.data$mat)))


# Filter poor genes
# Genes with UMI count greater than min.count = 2 in atleast min.cell = 3 cells is retained.
lnorm<-dropClust::normalize(filtered.data, min.count=2, min.cell=3)
sce.2<-FilterGenes(sce.1,min.count = 2, min.cell = 3)

sce.3<-CountNormalize(sce.2)

all.equal(as.matrix(normcounts(sce.3)), as.matrix(Matrix::t(lnorm$m)))
all.equal(sort(rowData(sce.3)$Symbol), sort(lnorm$use_genes))



# Select Top Dispersed Genes by setting ngenes_keep.
dp_genes <- dispersion_genes(lnorm, ngenes_keep = 1000)
sce.4<-RankGenes(sce.3)
all.equal(sort(rowData(sce.4)$Symbol[rowData(sce.4)$HVG]), sort(dp_genes))


# Log Normalize Matrix with genes-subset,
whole <- as.matrix(log2(lnorm$m[,rowData(sce.4)$HVG]+1))
sce.4.Log <-t(LogNormalize(normcounts(sce.4)[rowData(sce.4)$HVG,]+1))
all.equal(t(sce.4.Log), whole)



# Structure preserving Sampling
set.seed(0)
sample_ids = initial.samples(filtered.data)
samples_louvain<-dropClust::sampling(whole[sample_ids,])
subsamples_louvain<-sample_ids[samples_louvain]
set.seed(0)
sce.5<-Sampling(sce.4)

all.equal(which(sce.5$Sampling==T), sort(subsamples_louvain))


# Find PCA top 200 genes. This may take some time.

top_pc_genes<-dropClust::pc_genes(whole[subsamples_louvain,],top=200)

sce.6<-RankPCAGenes(sce.5)
all.equal(rowData(sce.6)$Symbol[which(rowData(sce.6)$PCAGenes==T)], colnames(whole)[sort(top_pc_genes)])


set.seed(0)
clust.list<-dropClust::cluster.cells(data = whole[,top_pc_genes], sp.samples = subsamples_louvain,
                          default = T, minClusterSize = 30,deepSplit = 2, conf = 0.8)
set.seed(0)
sce.7<-Cluster(sce.6, method = "default", minClusterSize = 30, deepSplit = 2, conf = 0.8)



# Construct sub-matrix for DE analysis
de.mat<- dropClust::reduce_mat_de(lnorm,clust.list)
# Cells of interest
GRP = levels(clust.list$cluster.ident)
DE_genes_nodes_all  <- dropClust::DE_genes(de_data = de.mat, selected_clusters = GRP, lfc_th = 1,q_th = 0.001)


deg.dc3<-FindMarkers(sce.7,selected_clusters = GRP, lfc_th = 1,q_th = 0.001)

for(i in 1:length(GRP))
 print(length(intersect(deg.dc3$genes.df[,i], DE_genes_nodes_all$genes.df[,i])))
