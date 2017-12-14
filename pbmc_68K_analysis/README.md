dropClust: Efficient clustering of ultra-large scRNA-seq data
================


###   Refer to the R markdown file of the 68K PBMC dataset analysis [here](https://debsin.github.io/dropClust/index.html)

1. The PBMC 68K dataset can be obtained from: <http://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc68k_rds/pbmc68k_data.rds> and put the pbmc68k\_data.rds file in the "data" directory. The annotation file has already been placed there.

2.  Execute dropClust\_main.R The script returns the predicted cluster IDs, 2D cluster map and a few intermediate results for further downstream analysis.

3.  Obtain the cell-type specific genes and the heatmap by executing the DE\_plot.R script.

4. The `seurat.R` script lists the commands for generating the Seurat based clustering outcome on the 68K PBMC dataset.

#### Executing the scripts 

The scripts must be copied into the main dropclust directory before executing them.

Note: We have tested the codes on R version 3.3 and Python 2.7