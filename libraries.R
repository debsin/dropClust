# install required libraries
list.of.packages <- c("Matrix", "ggplot2", "Rtsne", "svd", "dplyr", "plyr",
                      "data.table", "mclust", "flexclust", "reshape2", "fastcluster",
                      "irlba","dynamicTreeCut", "RColorBrewer","GenSA","gplots")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")




# Build louvain binaries
if(Sys.info()["sysname"]!="Windows"){
  wd = getwd()
  setwd(LOUVAIN_DIR)
  system("make")
  setwd(wd)
}

# ----------------------------
# load relevant libraries
# ----------------------------
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(dplyr)
library(plyr)
library(data.table)
library(mclust)
library(flexclust)
library(reshape2)
library(fastcluster)
library(irlba)
library(dynamicTreeCut)
library(RColorBrewer)
library(GenSA)
library(gplots)

# -------------------------------------
