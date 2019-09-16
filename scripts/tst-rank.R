library(DropletUtils)
library(SingleCellExperiment)
library("magrittr")
## RCA dataset

is_one_of <- function(x, classes) {
  stopifnot(is(classes, "character"))
  purrr::map_lgl(classes, function(class) is(x, class)) %>% any()
}

rca.data <- read.csv('C:/Projects/data/GSE81861_Cell_Line_COUNT.csv',header = T,row.names = 1,stringsAsFactors = F)
gene.symbols = make.unique(unlist(lapply(rownames(rca.data), function(x) toupper(unlist(strsplit(x, split = "_"))[2]))))

labels <- colnames(rca.data)
## Stripping cell id and batch no from each label
annotations <- unlist(lapply(labels, function(x) toupper(unlist(strsplit(x, split = "[_]"))[3])))
cellnames = unlist(lapply(labels, function(x) toupper(unlist(strsplit(x, split = "[_]"))[1])))

counts = Matrix:: Matrix(as.matrix(rca.data), sparse = T)

batch = rep(1, length(labels))
batch[grep('B2',labels)]<-2

source("compararisons/RCA/plots.R")

sce1 <- SingleCellExperiment(assays = list(counts = counts[,which(batch ==1)]))
rowData(sce1)$Symbol = gene.symbols
colData(sce1)$cell_line = annotations[which(batch ==1)]
colData(sce1)$Batch = batch[which(batch ==1)]

sce2 <- SingleCellExperiment(assays = list(counts = counts[,which(batch ==2)]))
rowData(sce2)$Symbol = gene.symbols
colData(sce2)$cell_line = annotations[which(batch ==2)]
colData(sce2)$Batch = batch[which(batch ==2)]


set.seed(0)
objects  = list(sce1,sce2)
all.objects = objects
use.de.genes = TRUE
all_de_genes = list()
merge_genes = list()

common_genes = rowData(objects[[1]])$Symbol

common_colData = names(colData(objects[[1]]))
for (obj.i in 1:length(objects)){
  cat("\nProcesing Batch",obj.i,"...\n")
  rowData(objects[[obj.i]])$ID = rownames(objects[[obj.i]])
  rownames(objects[[obj.i]]) = rowData(objects[[obj.i]])$Symbol
  common_genes = intersect(common_genes, rownames(objects[[obj.i]]))
  temp = objects[[obj.i]]

  temp<- FilterCells(temp)
  temp<- FilterGenes(temp, min.count = 1, min.cell = 1)
  temp<- CountNormalize(temp)
  temp<- RankGenes(temp,ngenes_keep = 1000)


  merge_genes[[obj.i]] = rowData(temp)$Symbol[rowData(temp)$HVG]


  common_colData = intersect(common_colData, names(colData(temp)))


  if(use.de.genes){
    temp<- Sampling(temp)
    temp<- RankPCAGenes(temp)
    temp<- Cluster(temp, method = "default", conf = 0.1, use.previous = F)

    GRP = sort(unique(temp$ClusterIDs))
    de_res = FindMarkers(temp,selected_clusters = GRP, lfc_th = 1.2,q_th = 0.001,nDE = 10)
    all_de_genes[[obj.i]] = de_res$genes.df
  }


}

intersect(sort(all_de_genes[[2]][,3]), sort(all_de_genes[[1]][,6]))

mix.data.counts  =  counts(objects[[1]])[common_genes, ]
# mix.data.ncounts = normcounts(objects[[1]])[common_genes, ]
col.data = colData(objects[[1]])[,common_colData]
batch = rep(1, ncol(mix.data.counts))



for(i in 2:length(objects)){
  exp_counts  <- counts(objects[[i]])[common_genes, ]
  # temp2 <-  normcounts(objects[[i]])[common_genes, ]
  mix.data.counts  <- SingleCellExperiment::cbind(mix.data.counts, exp_counts)
  # mix.data.ncounts <- SingleCellExperiment::cbind(mix.data.ncounts, temp2)
  # colData(objects[[i]])$Batch = rep(i, nrow(colData(objects[[i]])))
  col.data = SingleCellExperiment::rbind(col.data, colData(objects[[i]])[,common_colData])

  batch  = c(batch, rep(i, ncol(exp_counts)))
}


mix.data =  SingleCellExperiment(assays = list(counts = mix.data.counts),
                                 rowData = common_genes, colData = DataFrame(col.data))
rownames(mix.data) = common_genes
colData(mix.data)$Batch = batch


if(use.de.genes){
  common_genes_de = intersect(unique(as.vector(unlist(all_de_genes))), unique(unlist(common_genes)))
  rowData(mix.data)$CommonDEGenes = rep(FALSE, nrow(mix.data))
  rowData(mix.data)$CommonDEGenes[which(rownames(mix.data) %in% common_genes_de)]<- TRUE
}

mix.data = CountNormalize(mix.data)

# rank.mix.data.counts = apply(mix.data.counts, 2, rank)

my_genes.index = 1:nrow(mix.data)#which(rowData(mix.data)$CommonDEGenes==T)
# my_genes.index = which(rowData(mix.data)$CommonDEGenes==T)
mix.data.ncounts = as.matrix(normcounts(mix.data)[my_genes.index,])

rank.mix.data.ncounts <- matrix(nrow = dim(mix.data.ncounts)[2], ncol = length(my_genes.index))


for (iter in 1:length(my_genes.index)) {
  # Gives the names of the genes in sorted order based on the gene-expression value
  um_vec = mix.data.ncounts[iter,]
  rank_vec = rank(replace(um_vec,um_vec==0,NA), na='keep')

  rank.mix.data.ncounts[,iter]<-rank_vec

}

rank.mix.data.ncounts[is.na(rank.mix.data.ncounts)] <- 0

# sd.genes = apply(rank.mix.data.ncounts, 2, sd)

rank.mix.data.mat = rank.mix.data.ncounts[,which(rowData(mix.data)$CommonDEGenes==T)]
dim(rank.mix.data.mat)

# common_genes_de = setdiff(common_genes_de, union(all_de_genes[[2]][,3], all_de_genes[[1]][,6]))
# set.seed(0)
# PROJ_c_dc = Rtsne::Rtsne(rank.mix.data.mat, perplexity = 30, pca = FALSE, check_duplicates = FALSE)$Y

set.seed(0)
PROJ_c_dc = uwot::umap(rank.mix.data.mat,  min_dist = 0.5,
                       n_neighbors = 20, n_components = 10, metric = "cosine")


# PROJ_c_dc = Rtsne::Rtsne(reducedDim(f.out, "corrected"))$Y

plot_proj_df_true_dc<-data.frame(Y1 = PROJ_c_dc[,1],Y2 = PROJ_c_dc[,2],
                                 color = as.factor(mix.data$cell_line),
                                 batch = as.integer(mix.data$Batch))
batch_plot(plot_proj_df_true_dc,filename = NA, title = "Corrected DC ",  type=NULL)

