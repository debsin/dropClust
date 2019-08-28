

pca_genes_runs = list()

for(iter in 1:10){
  set.seed(iter*sample(1:100,1))
  x = RankPCAGenes(sce.3, prev = TRUE)

  pca_genes_runs[[iter]] = rowData(x)$Symbol[rowData(x)$PCAGenes]
}

intersect_nos = list()
k=0
for(i in 1:10){
  for(j in 1:10){
    if(i==j) next
    k=k+1
    intersect_nos[[k]] =  length(intersect(pca_genes_runs[[i]], pca_genes_runs[[j]]))
    if(intersect_nos[[k]] == 0) print(paste(i,j))
  }
}
plot(density(unlist(intersect_nos)))

table(unlist(intersect_nos))




library(DropletUtils)
library(SingleCellExperiment)
library("magrittr")


# set.seed(0)
sce <-readfiles(path = "C:/Projects/data/pbmc3K/hg19/")
annotations = read.csv("C:/Users/debajyoti/Desktop/68k_pbmc_barcodes_annotation.tsv",sep="\t",header = T)
# sce.copy<-sce


sce.1<-FilterCells(sce)
sce.1<-FilterGenes(sce.1)
sce.2<-CountNormalize(sce.1)

sce.2<-RankGenes(sce.2,ngenes_keep = 1000)

sce.3<-Sampling(sce.2)


variety = unlist(rep(list(c(TRUE,FALSE)), each = 3))

old.pca.genes = rowData(x)$Symbol[rowData(RankPCAGenes(sce.3, top = 500, prev =TRUE))$PCAGenes]
new.pca.genes = rowData(x)$Symbol[rowData(RankPCAGenes(sce.3, top= 500, prev = FALSE))$PCAGenes]

per.match = c()
for(i in seq(10, 500, by = 10)){
  matching = length(intersect(old.pca.genes[1:i], new.pca.genes[1:i]))
  per.match = c(per.match, (matching/i) *100)
}
plot(seq(10, 500, by = 10), per.match, type='l')

library(microbenchmark)
mbm <- microbenchmark("old" = { a <- RankPCAGenes(sce.3, prev =TRUE)},
                      "new" = {
                        b <- RankPCAGenes(sce.3, prev = FALSE)
                      }, times = 20)

mbm


sce.3<-RankPCAGenes(sce.3)

intersect_nos[[k]] =  length(intersect(pca_genes_runs[[i]], pca_genes_runs[[j]]))
