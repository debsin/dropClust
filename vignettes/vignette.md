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
#> Unassigned Cells 238 
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
#>       cluster_1   cluster_2  cluster_3 cluster_4   cluster_5      
#>  [1,] "ITGA2B"    "IL7R"     "LGALS2"  "CD79A"     "RP11-290F20.3"
#>  [2,] "GP9"       "CD3D"     "S100A8"  "MS4A1"     "MS4A7"        
#>  [3,] "SPARC"     "LDHB"     "S100A9"  "CD79B"     "LST1"         
#>  [4,] "TUBB1"     "CD3E"     "MS4A6A"  "TCL1A"     "IFITM3"       
#>  [5,] "GNG11"     "RPS25"    "LYZ"     "HLA-DQB1"  "CDKN1C"       
#>  [6,] "TMEM40"    "CD74"     "CD14"    "LINC00926" "AIF1"         
#>  [7,] "PTCRA"     "HLA-DRA"  "FCN1"    "HLA-DQA1"  "HES4"         
#>  [8,] "SDPR"      "LTB"      "GPX1"    "CD74"      "SERPINA1"     
#>  [9,] "PF4"       "HLA-DPA1" "CST3"    "VPREB3"    "LILRA3"       
#> [10,] "CLU"       "LEF1"     "TYROBP"  "HLA-DRA"   "COTL1"        
#> [11,] "CMTM5"     "CYBA"     "GSTP1"   "HLA-DRB1"  "FCER1G"       
#> [12,] "CD9"       "HLA-DRB1" "TYMP"    "HLA-DPB1"  "IFI30"        
#> [13,] "TREML1"    "TYROBP"   "MALAT1"  "S100A4"    "CFD"          
#> [14,] "MYL9"      "HLA-DPB1" "ALDH2"   "FCER2"     "HCK"          
#> [15,] "ACRBP"     "FCER1G"   "FOLR3"   "RPL18A"    "HMOX1"        
#> [16,] "CA2"       "AQP3"     "S100A6"  "TMSB4X"    "LRRC25"       
#> [17,] "PPBP"      "IL32"     "LGALS1"  "LTB"       "PILRA"        
#> [18,] "HIST1H2AC" "RPS29"    "BLVRB"   "SRGN"      "CEBPB"        
#> [19,] "NRGN"      "CD40LG"   "LGALS3"  "HLA-DPA1"  "FTH1"         
#> [20,] "RGS18"     "FLT3LG"   "GRN"     "RPS23"     "FCGR3A"       
#> [21,] "CLDN5"     "TRAT1"    "CSF3R"   "S100A6"    "PSAP"         
#> [22,] "NGFRAP1"   "MAL"      "FTL"     "BLK"       "LILRA5"       
#> [23,] "TSC22D1"   "HLA-DRB5" "AP1S2"   "CD37"      "LILRB2"       
#> [24,] "PGRMC1"    "CCR7"     "NUP214"  "HLA-DRB5"  "CD68"         
#> [25,] "TPM1"      "CD27"     "ASGR1"   "CD72"      "SIGLEC10"     
#> [26,] "GRAP2"     "NOSIP"    "B2M"     "HVCN1"     "SPI1"         
#> [27,] "RUFY1"     "TCF7"     "AIF1"    "KIAA0125"  "CKB"          
#> [28,] "TPM4"      "ARPC3"    "PTPRCAP" "BANK1"     "S100A11"      
#> [29,] "PLA2G12A"  "SRGN"     "FTH1"    "LGALS1"    "STXBP2"       
#> [30,] "MPP1"      "OAZ1"     "CEBPD"   "RPS5"      "CTSS"         
#>       cluster_6 cluster_7 
#>  [1,] "GNLY"    "GZMH"    
#>  [2,] "GZMB"    "CD8A"    
#>  [3,] "PRF1"    "CD3D"    
#>  [4,] "SPON2"   "CCL5"    
#>  [5,] "GZMA"    "CD8B"    
#>  [6,] "CTSW"    "CST7"    
#>  [7,] "FGFBP2"  "IL32"    
#>  [8,] "NKG7"    "NKG7"    
#>  [9,] "CST7"    "TIGIT"   
#> [10,] "CD247"   "KLRG1"   
#> [11,] "CLIC3"   "GZMA"    
#> [12,] "XCL2"    "LYAR"    
#> [13,] "GZMM"    "CD3G"    
#> [14,] "KLRD1"   "FGFBP2"  
#> [15,] "CD7"     "LCK"     
#> [16,] "CCL4"    "PTPRCAP" 
#> [17,] "AKR1C3"  "CTSW"    
#> [18,] "HOPX"    "TYROBP"  
#> [19,] "B2M"     "CD2"     
#> [20,] "APMAP"   "FCER1G"  
#> [21,] "HLA-A"   "GZMM"    
#> [22,] "HLA-C"   "FCRL6"   
#> [23,] "TTC38"   "FTL"     
#> [24,] "IGFBP7"  "CD99"    
#> [25,] "MATK"    "RARRES3" 
#> [26,] "XCL1"    "FTH1"    
#> [27,] "PLAC8"   "HLA-DRA" 
#> [28,] "HLA-B"   "CD3E"    
#> [29,] "FCGR3A"  "LST1"    
#> [30,] "RARRES3" "C1orf162"

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
