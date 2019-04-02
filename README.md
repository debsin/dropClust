
The latest version of dropClust is now available in desktop and online versions.


dropClust Online
====

  Visit [https://debsinha.shinyapps.io/dropClust/](https://debsinha.shinyapps.io/dropClust/) for the online version.

   -   [Installation](#desktop-installation)
   -   [Tutorial](#vignette-tutorial)
       -  [Setting-up](#setting-up-directories)
       -  [Loading data](#loading-data)
       -  [Pre-processing](#pre-processing)
       -  [Sampling](#structure-preserving-sampling)
       -  [Clustering](#clustering)
       -  [Visualizing](#visualizing-clusters)
       -  [Differential gene analysis](#find-cluster-specific-differentially-expressed-genes)
       -  [Plot marker genes](#plot-hand-picked-marker-genes)
       -  [Draw heatmap](#draw-heatmap)



Desktop Installation
===============

The developer version of the R package can be installed with the following R commands:

``` r
library(devtools)
install_github("debsin/dropClust", dependencies = T)
```

Vignette tutorial
------------------
This vignette uses a small data set from the 10X website (3K PBMC dataset [here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) ) to demonstrate a standard pipeline. This vignette can be used as a tutorial as well.

Setting up directories
----------------------

``` r

library(dropClust)

# -------------------------------------
# specify paths and load functions
# -------------------------------------
WORK_DIR = "C:/Projects/dropClust/"
DATA_DIR <- file.path(WORK_DIR,"data")       
FIG_DIR <-  file.path(WORK_DIR,"plots")         
REPORT_DIR  <- file.path(WORK_DIR,"report")    

dir.create(file.path(FIG_DIR),showWarnings = FALSE)
dir.create(file.path(REPORT_DIR),showWarnings = FALSE)

set.seed(0)
```

Loading data
------------

dropClust loads UMI count expression data from three input files. The files follow the same structure as the datasets available from the 10X website, i.e.:

-   count matrix file in sparse format
-   transcriptome identifiers as a TSV file and
-   gene identifiers as a TSV file

``` r
# Load Data from path: C:/Projects/dropClust/data/pbmc3k/hg19
pbmc.data<-read_data(file.path(DATA_DIR,'pbmc3k/hg19'), format = "10X")
```

Pre-processing
--------------

dropClust performs pre-processing to remove poor quality cells and genes. dropClust is also equipped to mitigate batch-effects that may be present. The user does not need to provide any information regarding the source of the batch for individual transcriptomes. However, the batch-effect removal step is optional.

Cells are filtered based on the total UMI count in a cell specified by parameter `th`.

``` r
# Filter poor quality cells.  A threshold th corresponds to the total count of a cell.
filtered.data <- filter_cells(pbmc.data)

dim(filtered.data$mat)
```

### Boosting rare cell discovery

dropClust estimates minor population and their associated genes directly from the raw data by finding exclusive groups of co-expressed genes within a small population of transcriptomes. These genes and cells are then included for later processing to boost rare-cell discovery.

``` r
print("Fetch rare genes...")
rare_data <- get_rare_genes(filtered.data)
```

### Data normalization and removing poor quality genes

Poor quality genes are removed based on the minimum number of cells with expressions above a given threshold. Count normalization is then performed with the good quality genes only.

``` r
# Filter poor genes
# Genes with UMI count greater than min.count = 2 in atleast min.cell = 3 cells is retained.
lnorm<-normalize(filtered.data, min.count=2, min.cell=3)
```

Further gene selection is carried out by ranking the genes based on its dispersion index.

```r
# Select Top Dispersed Genes by setting ngenes_keep.
dp_genes <- dispersion_genes(lnorm, ngenes_keep = 1000)
```

Finally log normalization is performed after adding pseudo count. In this stage batch correction my be carried out by setting the `matrix.transform()` function argument: `batch_correction=TRUE`.
```r
# Log Normalize Matrix with genes-subset,
# perform batch effect removal operation when input contains batches
whole <- matrix.transform(lnorm,dp_genes ,rare_data, batch_correction = FALSE)
```




Structure Preserving Sampling
-----------------------------

Primary clustering is performed in a fast manner to estimate a gross structure of the data. Each of these clusters is then sampled to fine tune the clustering process.

``` r

sample_ids = initial.samples(filtered.data, rare_data)


# Structure preserving Sampling

samples_louvain<-sampling(whole[sample_ids,])

subsamples_louvain<-sample_ids[samples_louvain]
```

### Gene selection based on PCA

Another gene selection is performed to reduce the number of dimensions. PCA is used to identify genes affecting major components.

``` r


# Find PCA top 200 genes. This may take some time.
top_pc_genes<-pc_genes(whole[subsamples_louvain,],top=200)
```


Clustering
------------------

### Fine tuning the clustering process

By default best-fit, Louvain based clusters are returned. However, the user can tune the parameters to produce the desired number of clusters. The un-sampled transcriptomes are assigned cluster identifiers from among those identifiers produced from fine-tuning clustering. The post-hoc assignment can be controlled by setting the confidence value `conf`. High `conf` values will assign cluster identifiers to only those transcriptomes sharing a majority of common nearest neighbours.

``` r

# Adjust Minimum cluster size with argument minClusterSize (default = 20)
# Adjust tree cut with argument level deepSplit (default = 3), higher value produces more clusters.

clust.list<-cluster.cells(data = whole[,top_pc_genes], sp.samples = subsamples_louvain,
                          default = FALSE, minClusterSize = 30,deepSplit = 2, conf = 0.8)
```

Visualizing clusters
--------------------

Compute 2D embeddings for samples followed by post-hoc clustering.

``` r
# ----------------------------
# Tsne & CO-ORDINATE Projection
# ----------------------------
PROJ <- compute_2d_embedding(data = whole[,top_pc_genes],
                             sp.samples = subsamples_louvain,
                             clust.list = clust.list)


# Scatter Plot

plot_proj_df_pred<-data.frame(Y1 = PROJ[,1],Y2 = PROJ[,2],color = as.factor(clust.list$cluster.ident))

p<-all_plot(plot_proj_df_pred,filename = NA, title = "dropClust clusters")

print(p)
```


Find cluster specific Differentially Expressed genes
----------------------------------------------------

``` r

# Construct sub-matrix for DE analysis
de.mat<- reduce_mat_de(lnorm,clust.list)
```

For performing DE analysis on selected clusters only, set variable `GRP` with a subset of predicted cluster identifiers as a character vector.
```r
# Pick Cell Type Specific Genes
#############################

# Cells of interest
GRP = levels(clust.list$cluster.ident)
# int_cells  = which(label %in% GRP)

DE_genes_nodes_all  <- DE_genes(de_data = de.mat, selected_clusters = GRP, lfc_th = 1,q_th = 0.001)

write.csv(DE_genes_nodes_all$genes,
          file = file.path(REPORT_DIR, "ct_genes.csv"),
          quote = FALSE)
```

Plot hand picked marker genes
-----------------------------

``` r

marker_genes = c("S100A8","GNLY","PF4" )

p<-plot_markers(de_data = de.mat, marker_genes)
```


Heat map of top DE genes from each cluster
------------------------------------------

``` r
# Draw heatmap
#############################
p<-plot_heatmap(de_data = de.mat, DE_res = DE_genes_nodes_all$DE_res,nDE = 10)

print(p)
```
