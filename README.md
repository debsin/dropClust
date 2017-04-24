# dropClust

## Executing the script (Linux only):

1. Download pbmc 68K dataset from:
http://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc68k_rds/pbmc68k_data.rds
and put the pbmc68k_data.rds file in the "data" directory. The annotation file has already been placed there.

2. Execute dropClust_main.R
Update the working directory and data directory paths in the script.

3. Obtain heatmap by executing the downstream_strict.R script.


## Executing the script (Windows):
Build the louvain source files inside the "louvain" directory before proceeding with dropClust_main.R script.
