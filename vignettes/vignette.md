dropClust
================

This latest version of dropClust has been re-implemented to include new features and performance improvements. This vignette uses a small data set from the 10X website (3K PBMC dataset [here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) ) to demonstrate a standard pipeline. This vignette can be used as a tutorial as well.

Setting up directories
----------------------

``` r

library(dropClust)
# -------------------------------------
# specify paths and load functions
# -------------------------------------
WORK_DIR = "E:/Projects/dropClust/"
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
# Load Data

pbmc.data<-read_data(file.path(DATA_DIR,'pbmc3k/hg19'), format = "10X")
```

Pre-processing
--------------

dropClust performs pre-processing to remove poor quality cells and genes. dropClust is also equipped to mitigate batch-effects that may be present. The user does not need to provide any information regarding the source of the batch for individual transcriptomes. However, the batch-effect removal step is optional.

Cells are filtered based on the total UMI count in a cell specified by parameter `th`.

``` r
# Filter poor quality cells.  A threshold th corresponds to the total count of a cell.
filtered.data <- filter_cells(pbmc.data)
#> 3 bad cells present.

dim(filtered.data$mat)
#> [1]  2697 32738
```

### Boosting rare cell discovery

dropClust estimates minor population and their associated genes directly from the raw data by finding exclusive groups of co-expressed genes within a small population of transcriptomes. These genes and cells are then included for later processing to boost rare-cell discovery.

``` r
print("Fetch rare genes...")
#> [1] "Fetch rare genes..."
rare_data <- get_rare_genes(filtered.data)
#> Finding co-expressed rare genes...
```

### Data normalization and removing poor quality genes

Poor quality genes are removed based on the minimum number of cells with expressions above a given threshold. Count normalization is then performed with the good quality genes only.

Further gene selection is carried out by ranking the genes based on its dispersion index.

``` r
# Filter poor genes
# Genes with UMI count greater than min.count = 2 in atleast min.cell = 3 cells is retained.
lnorm<-normalize(filtered.data, min.count=2, min.cell=3)
#> Discard poor genes...
#> Dimensions of filtered Matrix: 2697 2233 
#> Median normalize matrix...


# Select Top Dispersed Genes by setting ngenes_keep.
dp_genes <- dispersion_genes(lnorm, ngenes_keep = 1000)
#> Sort Top Genes...
#> Cutoff Genes...


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
#> Build graph...
#> Louvain Partition...
subsamples_louvain<-sample_ids[samples_louvain]
```

### Gene selection based on PCA

Another gene selection is performed to reduce the number of dimensions. PCA is used to identify genes affecting major components.

``` r


# Find PCA top 200 genes. This may take some time.
top_pc_genes<-pc_genes(whole[subsamples_louvain,],top=200)
#> Find best PCA components...
```

![](vignette_files/figure-markdown_github/unnamed-chunk-7-1.png)

Perform clustering
------------------

### Fine tuning the clustering process

By default best-fit, Louvain based clusters are returned. However, the user can tune the parameters to produce the desired number of clusters. The un-sampled transcriptomes are assigned cluster identifiers from among those identifiers produced from fine-tuning clustering. The post-hoc assignment can be controlled by setting the confidence value `conf`. High `conf` values will assign cluster identifiers to only those transcriptomes sharing a majority of common nearest neighbours.

``` r

# Adjust Minimum cluster size with argument minClusterSize (default = 20)
# Adjust tree cut with argument level deepSplit (default = 3), higher value produces more clusters.

clust.list<-cluster.cells(data = whole[,top_pc_genes], sp.samples = subsamples_louvain,
                          default = FALSE, minClusterSize = 30,deepSplit = 2, conf = 0.8)
#> Perfom Hierarchical Clustering...
#> Build Graph...
#> Find ANN...
#> Assign Cluster Ids...
#> Unassigned Cells 260 
#> Number of Predicted Clusters: 7
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

![](vignette_files/figure-markdown_github/unnamed-chunk-9-1.png)

Find cluster specific Differentially Expressed genes
----------------------------------------------------

``` r

# Construct sub matrix for DE analysis
de.mat<- reduce_mat_de(lnorm,clust.list)



# Pick Cell Type Specific Genes
#############################

# Cells of interest
GRP = levels(clust.list$cluster.ident)
# int_cells  = which(label %in% GRP)

DE_genes_nodes_all  <- DE_genes(de_data = de.mat, selected_clusters = GRP, lfc_th = 1,q_th = 0.001)
#> 
#> Computing for 7 clusters:
#> 
#> Cluster 1 :
#> 
#> Computing Wilcoxon p values...
#> Computing Log fold change values...
#> 
#> Completed successfully.
#> 
#> Cluster 2 :
#> 
#> Computing Wilcoxon p values...
#> Computing Log fold change values...
#> 
#> Completed successfully.
#> 
#> Cluster 3 :
#> 
#> Computing Wilcoxon p values...
#> Computing Log fold change values...
#> 
#> Completed successfully.
#> 
#> Cluster 4 :
#> 
#> Computing Wilcoxon p values...
#> Computing Log fold change values...
#> 
#> Completed successfully.
#> 
#> Cluster 5 :
#> 
#> Computing Wilcoxon p values...
#> Computing Log fold change values...
#> 
#> Completed successfully.
#> 
#> Cluster 6 :
#> 
#> Computing Wilcoxon p values...
#> Computing Log fold change values...
#> 
#> Completed successfully.
#> 
#> Cluster 7 :
#> 
#> Computing Wilcoxon p values...
#> Computing Log fold change values...
#> 
#> Completed successfully.

DE_genes_nodes_all$genes
#>       cluster_1   cluster_2   cluster_3 cluster_4   cluster_5      
#>  [1,] "ITGA2B"    "CD3D"      "LGALS2"  "CD79A"     "MS4A7"        
#>  [2,] "GP9"       "IL7R"      "S100A8"  "MS4A1"     "LST1"         
#>  [3,] "SPARC"     "LDHB"      "S100A9"  "LINC00926" "HES4"         
#>  [4,] "TMEM40"    "IL32"      "LYZ"     "HLA-DQA1"  "RP11-290F20.3"
#>  [5,] "TUBB1"     "CD3E"      "CD14"    "TCL1A"     "IFITM3"       
#>  [6,] "GNG11"     "RPS25"     "FCN1"    "HLA-DQB1"  "AIF1"         
#>  [7,] "CLU"       "LTB"       "MS4A6A"  "CD74"      "CDKN1C"       
#>  [8,] "CD9"       "AQP3"      "CST3"    "VPREB3"    "FCER1G"       
#>  [9,] "SDPR"      "CD2"       "GRN"     "CD79B"     "SERPINA1"     
#> [10,] "TREML1"    "TRAT1"     "GPX1"    "HLA-DRA"   "FCGR3A"       
#> [11,] "PTCRA"     "RPS29"     "CSF3R"   "HLA-DPB1"  "COTL1"        
#> [12,] "PF4"       "JUNB"      "GSTP1"   "CD37"      "CFD"          
#> [13,] "PPBP"      "CD27"      "TYROBP"  "FCER2"     "SPI1"         
#> [14,] "CMTM5"     "FLT3LG"    "ALDH2"   "HLA-DRB1"  "CD68"         
#> [15,] "MYL9"      "CCR7"      "FOLR3"   "FCRLA"     "LRRC25"       
#> [16,] "NRGN"      "NPM1"      "AP1S2"   "BANK1"     "STXBP2"       
#> [17,] "RGS18"     "MAL"       "S100A6"  "HLA-DQA2"  "CKB"          
#> [18,] "ACRBP"     "RGCC"      "BST1"    "RPS23"     "PILRA"        
#> [19,] "NGFRAP1"   "GIMAP5"    "ASGR1"   "HLA-DPA1"  "CTSL"         
#> [20,] "HIST1H2AC" "NOSIP"     "LGALS1"  "CD72"      "FTH1"         
#> [21,] "CA2"       "CD3G"      "FTL"     "RPS5"      "SAT1"         
#> [22,] "CLDN5"     "TMEM66"    "TYMP"    "LTB"       "TIMP1"        
#> [23,] "TPM1"      "PRKCQ-AS1" "LGALS3"  "P2RX5"     "LILRA3"       
#> [24,] "PGRMC1"    "LCK"       "BLVRB"   "HLA-DMB"   "PSAP"         
#> [25,] "GRAP2"     "PIK3IP1"   "GABARAP" "HLA-DRB5"  "CEBPB"        
#> [26,] "RUFY1"     "AES"       "FPR1"    "HLA-DMA"   "LILRB2"       
#> [27,] "TPM4"      "ZFP36L2"   "FCGRT"   "HVCN1"     "IFITM2"       
#> [28,] "TSC22D1"   "TMEM123"   "CEBPD"   "HLA-DOB"   "CTSS"         
#> [29,] "MPP1"      "JUN"       "CFP"     "KIAA0125"  "S100A11"      
#> [30,] "PLA2G12A"  "CD40LG"    "S100A12" "BLNK"      "CSF1R"        
#>       cluster_6 cluster_7
#>  [1,] "GZMB"    "GZMH"   
#>  [2,] "GNLY"    "CD3D"   
#>  [3,] "PRF1"    "CCL5"   
#>  [4,] "SPON2"   "CD8B"   
#>  [5,] "FGFBP2"  "CD8A"   
#>  [6,] "GZMA"    "CST7"   
#>  [7,] "CST7"    "CD3G"   
#>  [8,] "CTSW"    "IL32"   
#>  [9,] "NKG7"    "GZMA"   
#> [10,] "HOPX"    "NKG7"   
#> [11,] "CD247"   "TIGIT"  
#> [12,] "CCL4"    "FCRL6"  
#> [13,] "KLRD1"   "FGFBP2" 
#> [14,] "CLIC3"   "CCL4"   
#> [15,] "CD7"     "GZMM"   
#> [16,] "XCL2"    "LCK"    
#> [17,] "GZMM"    "PTPRCAP"
#> [18,] "AKR1C3"  "CTSW"   
#> [19,] "APMAP"   "KLRG1"  
#> [20,] "B2M"     "PTPN7"  
#> [21,] "HLA-C"   "MYL12B" 
#> [22,] "IGFBP7"  "CHST12" 
#> [23,] "HLA-A"   ""       
#> [24,] "TTC38"   ""       
#> [25,] "PLAC8"   ""       
#> [26,] "RARRES3" ""       
#> [27,] "XCL1"    ""       
#> [28,] "GZMH"    ""       
#> [29,] "UBB"     ""       
#> [30,] "HLA-B"   ""

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

![](vignette_files/figure-markdown_github/unnamed-chunk-11-1.png)

Heat map of top DE genes from each cluster
------------------------------------------

``` r
# Draw heatmap
#############################
p<-plot_heatmap(de_data = de.mat, DE_res = DE_genes_nodes_all$DE_res,nDE = 10)

print(p)
```

![](vignette_files/figure-markdown_github/unnamed-chunk-12-1.png)
