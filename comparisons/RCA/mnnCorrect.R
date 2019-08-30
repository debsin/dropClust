# This script prepares data for the RCA data analysis.

##########################################
# Process batch #1

B1<- FilterGenes(objects[[1]])
B1<- CountNormalize(B1)
B1<- RankGenes(B1)
data1 = as.matrix(normcounts(B1)[rowData(B1)$HVG,])



##########################################
# Process batch #1

B2<- FilterGenes(objects[[2]])
B2<- CountNormalize(B2)
B2<- RankGenes(B2)
data2 = as.matrix(normcounts(B2)[rowData(B2)$HVG,])



##########################################
##########################################


keep <- intersect(rownames(data1) , rownames(data2))

data1 <- data1[keep,]
data2 <- data2[keep,]

mnn_mat  = batchelor::fastMNN(data1, data2, k=20,pc.input = F)
mnn_corr = reducedDim(mnn_mat)
