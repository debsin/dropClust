Integrative analysis
================

## Loading datasets

Each dataset represents one batch and must be a `SingleCellExperiment`
object. The objects are are merged by passing a list in the next step.

``` r

library(dropClust)
load(url("https://raw.githubusercontent.com/LuyiTian/CellBench_data/master/data/sincell_with_class.RData"))

objects = list()

objects[[1]] = sce_sc_10x_qc

objects[[2]] = sce_sc_CELseq2_qc

objects[[3]] = sce_sc_Dropseq_qc
```

## Merge datasets using dropClust

Datasets can be merged in two ways: using a set of DE genes from each
batch or, using the union of the sets of highly variable genes from each
batch.

    #> 
    #> Procesing Batch 1 ...
    #> 1 bad cells removed.
    #> 257 genes filtered out, 16211 genes remaining.
    #> Sort Top Genes...
    #> Cutoff Genes...
    #> Building graph with 901 nodes...Louvain Partition...Done.
    #> 700 samples extracted.
    #> Find best PCA components...[1]  700 1000
    #> 200 genes selected.
    #> 700 samples and 200 genes used for clustering.
    #> Build Graph with 700 samples...Done.
    #> Louvain Partitioning...Done.
    #> Find nearest neighbours among sub-samples...Done.
    #> Post-hoc Cluster Assignment...Done.
    #> Unassigned Cells 0 
    #> Number of Predicted Clusters: 6 
    #> Computing for DE genes:
    #> Cluster 1 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 2 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 3 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 4 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 5 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 6 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> 
    #> Procesing Batch 2 ...
    #> 1 bad cells removed.
    #> 11931 genes filtered out, 16273 genes remaining.
    #> Sort Top Genes...
    #> Cutoff Genes...
    #> Find best PCA components...[1]  273 1000
    #> 200 genes selected.
    #> 273 samples and 200 genes used for clustering.
    #> Build Graph with 273 samples...Done.
    #> Louvain Partitioning...Done.
    #> Find nearest neighbours among sub-samples...Done.
    #> Post-hoc Cluster Assignment...Done.
    #> Unassigned Cells 0 
    #> Number of Predicted Clusters: 4 
    #> Computing for DE genes:
    #> Cluster 1 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 2 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 3 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 4 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> 
    #> Procesing Batch 3 ...
    #> 1 bad cells removed.
    #> 958 genes filtered out, 14169 genes remaining.
    #> Sort Top Genes...
    #> Cutoff Genes...
    #> Find best PCA components...[1]  224 1000
    #> 200 genes selected.
    #> 224 samples and 200 genes used for clustering.
    #> Build Graph with 224 samples...Done.
    #> Louvain Partitioning...Done.
    #> Find nearest neighbours among sub-samples...Done.
    #> Post-hoc Cluster Assignment...Done.
    #> Unassigned Cells 0 
    #> Number of Predicted Clusters: 5 
    #> Computing for DE genes:
    #> Cluster 1 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 2 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 3 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 4 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.
    #> Cluster 5 :
    #>  Computing Wilcoxon p values...Done.
    #>  Computing Log fold change values...Done.

## Perform correction and dimension reduction

``` r
set.seed(1)
dc.corr <-  Correction(merged_data,  method="default", close_th = 0.1, cells_th = 0.1,
                       components = 10, n_neighbors = 20,  min_dist = 0.5)
#> Batch correcting...
#> from 172 to 97 genes.
#> Embedding with UMAP...Done
```

## Perform Clustering on integrated dimensions

``` r
dc.corr = Cluster(dc.corr,method = "kmeans",centers = 3)
#> Clustering on embedded dimensions...Perfom k-means Clustering...Done.
#> Done.
```

## Visualizing clusters

Compute 2D embeddings for samples followed by post-hoc clustering.

``` r
ScatterPlot(dc.corr, title = "Clusters")
```

![Batch corrected dropClust based
Clustering.](batchCorrection_files/figure-gfm/unnamed-chunk-5-1.png)
\#\# Optional Batch correction Users can use `fastmnn` method for batch
correction. Specific arguments of fastmnn can also be passed through the
`Correction` module.

``` r
merged_data.fastmnn<-Merge(all.objects,use.de.genes = FALSE)
#> 
#> Procesing Batch 1 ...
#> 1 bad cells removed.
#> 257 genes filtered out, 16211 genes remaining.
#> Sort Top Genes...
#> Cutoff Genes...
#> 
#> Procesing Batch 2 ...
#> 1 bad cells removed.
#> 11931 genes filtered out, 16273 genes remaining.
#> Sort Top Genes...
#> Cutoff Genes...
#> 
#> Procesing Batch 3 ...
#> 1 bad cells removed.
#> 958 genes filtered out, 14169 genes remaining.
#> Sort Top Genes...
#> Cutoff Genes...

set.seed(1)
mnn.corr <-  Correction(merged_data.fastmnn,  method="fastmnn", d = 10)

mnn.corr = Cluster(mnn.corr,method = "kmeans",centers = 3)
#> Clustering on embedded dimensions...Perfom k-means Clustering...Done.
#> Done.

ScatterPlot(mnn.corr, title = "Clusters")
```

![](batchCorrection_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Marker discovery from the merged dataset

``` r
de<-FindMarkers(dc.corr,q_th = 0.001, lfc_th = 1.2,nDE = 10)
#> Computing for DE genes:
#> Cluster 1 :
#>  Computing Wilcoxon p values...Done.
#>  Computing Log fold change values...Done.
#> Cluster 2 :
#>  Computing Wilcoxon p values...Done.
#>  Computing Log fold change values...Done.
#> Cluster 3 :
#>  Computing Wilcoxon p values...Done.
#>  Computing Log fold change values...Done.
de$genes.df
#>       cluster_1         cluster_2         cluster_3        
#>  [1,] "ENSG00000179344" "ENSG00000221923" "ENSG00000122406"
#>  [2,] "ENSG00000121858" "ENSG00000106211" "ENSG00000065978"
#>  [3,] "ENSG00000214548" "ENSG00000100867" "ENSG00000142864"
#>  [4,] "ENSG00000131203" "ENSG00000242265" "ENSG00000253706"
#>  [5,] "ENSG00000151632" "ENSG00000196109" "ENSG00000106153"
#>  [6,] "ENSG00000223865" "ENSG00000129991" "ENSG00000258484"
#>  [7,] "ENSG00000198502" "ENSG00000105568" "ENSG00000162433"
#>  [8,] "ENSG00000100234" "ENSG00000178150" "ENSG00000132432"
#>  [9,] "ENSG00000154277" "ENSG00000130724" "ENSG00000167641"
#> [10,] "ENSG00000134339" "ENSG00000170906" "ENSG00000256670"
```
