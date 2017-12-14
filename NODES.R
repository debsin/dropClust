if (!require("MetaDE")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("impute")
  biocLite("Biobase")
  install.packages("MetaDE",repos = "http://cran.us.r-project.org")
}


## NODES
NODES1<-function(data,group,r=20,smooth_points = 10000, zper = 0.5)
{
  

  
  require(MetaDE)
  
  #This part is for identifying the groups
  indices<-list()
  
  U<-unique(group)
  
  for(i in 1:length(U))
  {
    indices[[as.character(U[i])]] <- grep(U[i],group)
  }
  
  # length of group 1
  n1<-length(indices[[U[1]]])
  n2<-length(indices[[U[2]]])
  
  
  # Getting Noise distribution from NOISeq and estimate p values
  
  Zr=NULL
  for (i in 1:r) {
    
    print(paste("Randomization run =", i))
    
    mipermu = sample(1:(n1+n2)) ## randomize labels
    
    mipermu = data[,mipermu] ## randomize matrix columns accordingly
    
    mean1 = rowMeans(as.matrix(mipermu[,1:n1]) )## get the means for random group 1
    mean2 = rowMeans(as.matrix(mipermu[,(n1+1):(n1+n2)]) ) ## get the means for random group 2
    
    sd1 = apply(mipermu[,1:n1], 1, sd) ## sd for group 1
    sd2 = apply(mipermu[,(n1+1):(n1+n2)], 1, sd) ## sd for group 2
    
    myparam = list("n" = c(n1,n2), "sd" = cbind(sd1,sd2))
    
    MDperm <- MDbio(dat = cbind(mean1, mean2), param = myparam, a0per = zper)
    
    Zr = cbind(Zr, MDperm$D)
    
  }
  
  
  # Estimating noise density using Gaussian kernel
  cat("\nSmoothing the noise density...\n")
  dF <- approxfun(density(as.vector(Zr),n = smooth_points))
  
  # Getting stat for all genes
  
  mean1 = rowMeans(as.matrix(data[,indices[[U[1]]]]))
  mean2 = rowMeans(as.matrix(data[,indices[[U[2]]]]))
  
  sd1 = apply(as.matrix(data[,indices[[U[1]]]]), 1, sd)
  sd2 = apply(as.matrix(data[,indices[[U[2]]]]), 1, sd)
  
  myparam = list("n" = c(n1,n2), "sd" = cbind(sd1,sd2))
  
  Ds <- MDbio(dat = cbind(mean1, mean2), param = myparam, a0per = zper)
  
  Zs = Ds$D
  
  
  # Estimating p values only from noise distribution
  
  prob_1<-apply(as.matrix(Zs),1,function(x) den(dF,as.vector(Zr),x,smooth_points))
  
  
  
  # getting wilcoxon p values
  cat("\nComputing Wilcoxon p values...\n")
  prob_2<-apply(data,1,function(x) wilcox.test(x[indices[[U[1]]]],x[indices[[U[2]]]])$p.value)
  
  # Fisher's method to combine p values
  cat("\nCombining p values using Fisher's method...\n")
  
  together <- list()
  
  
  together[['p']] <- as.matrix(cbind(prob_1,prob_2))
  META<-MetaDE::MetaDE.pvalue(together,meta.method="Fisher")
  PVAL<- META$meta.analysis$pval
  #names(PVAL)<-rownames(data)
  
  
  # fdr
  fdr <- p.adjust(PVAL,method="fdr")
  
  # prepare final result
  res <- data.frame(cbind(pvalues=PVAL,qvalues=fdr))
  
  rownames(res) <- rownames(data)
  
  cat("\nCompleted successfully.\n")
  
  # return
  return(res)
}






# Tail probability estimation

den<- function (approx, obs, val, points)
{
  if (val >= max(obs)) {
    pF = 1/(1 + points)
  }
  else {
    pF <- integrate(approx, val, max(obs))$value
  }
  return(pF)
}









## This code is a shadow of NOISeqbio


MDbio = function (dat = dat, param = NULL, a0per = 0.5) {
  
  #dd <- (dat[,1]-dat[,2]) this is for two tailed.
  dd <- abs(dat[,1]-dat[,2])
  
  
  #sd.D = sqrt(param$sd[,1]^2/sqrt(param$n[1]) + param$sd[,2]^2/sqrt(param$n[2]))
  sd.D = sqrt( (param$sd[,1]^2/param$n[1])
               +
                 (param$sd[,2]^2/param$n[2])
  )
  
  a0per = as.numeric(a0per)
  a0.D <- quantile(sd.D, probs = a0per, na.rm = TRUE)
  
  dd <- dd / (a0.D + sd.D)
  
  # Results
  list("D" = dd)
}




MDbio <- function (dat = dat, param = NULL, a0per = 0.5)
{
  dd <- abs(dat[, 1] - dat[, 2])
  sd.D = sqrt((param$sd[, 1]^2/param$n[1]) + (param$sd[, 2]^2/param$n[2]))
  a0per = as.numeric(a0per)
  a0.D <- quantile(sd.D, probs = a0per, na.rm = TRUE)
  dd <- dd/(a0.D + sd.D)
  list(D = dd)
}

