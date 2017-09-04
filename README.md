# dropClust: Efficient clustering of ultra-large scRNA-seq data
Authors: _Debajyoti Sinha, Akhilesh Kumar, Himanshu Kumar, Sanghamitra Bandyopadhyay, Debarka Sengupta_

### The preprint version of the software is aviable online at 
http://www.biorxiv.org/content/early/2017/07/31/170308


### Refer to the R execution markdown file on 68K PBMC data [here](https://debsin.github.io/dropClust/index.html)


## Prerequisites:
1. Python  (>=2.7), R
2. Python sklearn package
3. R dependencies will be installed automatically when executing the main script.
4. C++ compiler for Windows users.


## Executing the script (Linux:Ubuntu):

1. Download PBMC 68K dataset from:
http://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc68k_rds/pbmc68k_data.rds
and put the pbmc68k_data.rds file in the "data" directory. The annotation file has already been placed there.

2. Execute dropClust_main.R
    The script returns the predicted cluster IDs, 2D cluster map and a few intermediate results for further downstream analysis.

3. Obtain the cell-type specific genes and the heatmap by executing the DE_plot.R script.


## Executing the script (Windows or other Linux distributions):

Build windows executable/binaries from "louvain/src" files inside the "louvain" directory before executing the main R script.
