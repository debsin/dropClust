source("NODES.r")
##

DE_genes <- function(data,labels,max )
{
  #dim(data)
  
  raw_data = 2^data
  #min(raw_data)
  # throw away 0 marked outlier cells
  ind = which(labels == 0)
  
  
  if(length(ind) !=0)
  {
    data = data[,-ind]
    labels = labels[-ind]
  }
  
  min_cell = floor(ncol(raw_data)*0.005)
  keep_genes = apply(raw_data, 1 , function(x)  sum(x>3))>=min_cell
  
  #table(omit_genes)
  
  raw_data = raw_data[keep_genes,]
  
  # call NODES for each class vs class
  DE_list = list()
  pairs = combn(unique(labels),2)
  pair_ids = paste(pairs[1,],pairs[2,],sep="_")
  
  DE_list<-foreach(i = 1:dim(pairs)[2]) %dopar%
  {
    IND_a = which(labels == pairs[1,i])
    #IND_b = setdiff(1:ncol(data),IND_a)
    IND_b = which(labels == pairs[2,i])
    
    
    DE_res = NODES1(raw_data[,c(IND_a,IND_b)],group =c(rep("A",length(IND_a)),rep("B",length(IND_b))))
    LFC = log2(rowMeans(raw_data[,IND_a])/rowMeans(raw_data[,IND_b]))
    
    
    #length(names(LFC)==rownames(DE_res))
    #head(DE_res)
    #rN = rownames(DE_res)[head(intersect(which(abs(LFC)>=1), which(DE_res$qvalues<=0.05)),30)]
    sig = DE_res[which(DE_res$qvalues<=0.05),]
    if(max>0){
      rN = head(rownames(sig)[order(sig$qvalues)],max)
    }
    else{
      rN = rownames(sig)[order(sig$qvalues)]
    }
    
    #DE_up = union(DE_up,rN)
    DE_list[[pair_ids[i]]] = data.frame(gene = rN,q_val = DE_res$qvalues[match(rN,rownames(DE_res))],fc = LFC[match(rN,names(LFC))])
    #return(data.frame(gene = rN,q_val = DE_res$qvalues[match(rN,rownames(DE_res))],fc = LFC[match(rN,names(LFC))]))
    
  }
  
  
  DE_up= c()
  for(i in 1:length(DE_list)){
    DE_up = union(DE_up, DE_list[[i]]$gene)
  }
  gc()
  names(DE_list) = pair_ids
  RES = list(genes = DE_up,DE_res = DE_list)
  return(RES)
}


plot_heat_genes <- function(data,labels,meth)
{
  
  avg_mat = c()
  for(i in 1:length(unique(labels)))
  {
    
    rM = rowMeans(data[,which(labels == unique(labels)[i])])
    avg_mat = cbind(avg_mat,rM)
  }
  
  colnames(avg_mat) = unique(labels)
  
  
  #View(avg_mat)
  
  
  H = hclust(dist(t(avg_mat)),method = "average")
  #H$order
  
  IND = c()
  
  
  for(i in 1:length(H$order))
  {
    
    IND = append(IND,which(labels == H$order[i]))
    
  }
  
  
  #length(IND)
  
  # reorder data and labels
  labels = labels[IND]
  data = data[,IND]
  
  
  # get colors accordingly
  #myPalette1 <- colorRampPalette(brewer.pal(9,"Set1"))
  #colors = myPalette1(length(unique(labels)))
  
  
  
  
  
  
  
  return(list(mat=as.matrix(data), lab =labels ))
}