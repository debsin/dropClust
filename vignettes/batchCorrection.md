Batch Correction
================
Debajyoti Sinha

## Integrative analysis

``` r

library(dropClust)
load(url("https://raw.githubusercontent.com/LuyiTian/CellBench_data/master/data/sincell_with_class.RData"))

set.seed(0)
objects = list()

objects[[1]] = sce_sc_10x_qc

objects[[2]] = sce_sc_CELseq2_qc

objects[[3]] = sce_sc_Dropseq_qc
```

## Merge datasets using dropClust

Datasets can be merged in two ways: using a set of DE genes from each
batch or, using the common set of genes from the raw count data

``` r

all.objects = objects
merged_data<-Merge(all.objects)

annotations = merged_data$cell_line
batch.id = merged_data$Batch
```

## Perform correction and dimension reduction

``` r

dc.corr <- Correction(merged_data, close_th = 0.1, cells_th = 0.1,
                       components = 10, n_neighbors = 30, init = "spca")
```

## Perform Clustering on integrated dimensions

``` r
cc = Cluster(dc.corr,method = "kmeans",centers = 3)
```

## Marker discovery from the merged dataset

``` r
de<-FindMarkers(cc,q_th = 0.001, lfc_th = 1.2)
```
