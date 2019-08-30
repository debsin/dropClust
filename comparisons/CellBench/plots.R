library(bioDist)
library(ggplot2)
plot2D<-function(data,...){

  # d = bioDist::spearman.dist(data)
  d = dist(data)
  pj = Rtsne::Rtsne(d,is_distance = T)
  return(pj)
}


# --------------------------------------------------
# Return custom colors
# --------------------------------------------------
#' Select custom colors
#' @description Produces custom colors.
#' @param n number of colors to generate.
#' @return vector of colour hexcodes of length n.

getColors<-function(n){
  mycolors = toupper(c("#00fe0a", "#ff0000", "#bded1b", "#794b05", "#c3beb7",
               "#0000ff", "#00ffff","#ff21d3" , "#81b7dd","#f87791" ,
               "#1e7309", "#fc9a07", "#625b51", "#6a09c3", "#189ff5",
               "#d19d00", "#0ebf06", "#88ffb3", "#f6fc2a", "#000000"))
  if(n<=9){
    # cat("Too many colors...Using fallback color scheme.\n")
    getPalette = grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(9, "Set1"))
    return(getPalette(n))
  }else if(n<=12){
    # cat("Too many colors...Using fallback color scheme.\n")
    getPalette = grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(12, "Paired"))
    return(getPalette(n))
  }else if(n<=20){
    # cat("Too many colors...Using fallback color scheme.\n")
    return(mycolors[1:n])
  }
  else if(n>20){
    cat("Too many colors...Using default fallback color scheme.\n")
    getPalette =
      grDevices::colorRampPalette(
        suppressWarnings(RColorBrewer::brewer.pal(n, "Set1")))
    return(getPalette(n))
  }
  return()
}


# --------------------------------------------------
# Scatter plot of cells in two dimensions
# --------------------------------------------------
#' Plotting cells
#' @description Scatter plot of cells in two dimensions
#' @param data data frame containing 3 columns:\cr
#' \enumerate{
#'     \item \code{Y1}: numeric values of each point in x-axis.
#'     \item \code{Y2}: numeric values of each point in y-axis.
#'     \item \code{color}: factor, cluster identifiers for each point.
#' }
#' @param filename (optional) specify file path to save plot in pdf format.
#' @param title character, specify plot title.
#' @return grob object
#' @export
#' @importFrom ggplot2 ggplot geom_point aes scale_colour_manual theme_classic theme ggtitle annotate guides guide_legend ylab xlab aes_string
batch_plot<-function(data,filename=NA, title, type=NULL){

  temp = stats::complete.cases(data)
  plot_proj_df = data[temp, ]

  plot_proj_df$color<-factor(plot_proj_df$color)

  plot_proj_df$batch<-factor(as.numeric(plot_proj_df$batch))

  x.mean = stats::aggregate(plot_proj_df$Y1,
                            list(plot_proj_df$color),
                            stats::median)[,-1]
  y.mean = stats::aggregate(plot_proj_df$Y2,
                            list(plot_proj_df$color),
                            stats::median)[,-1]

  colorcount_t = length(levels(plot_proj_df$color))


  n_points = dim(plot_proj_df[stats::complete.cases(plot_proj_df),])[1]


  p_size = 1.5
  p<-ggplot()
  p2<-p+ geom_point(data = plot_proj_df,aes_string(x ='Y1',y = 'Y2', color = 'color',shape='batch'),size=p_size)  +
    scale_color_manual(values = getColors(  colorcount_t ))+
    theme_classic()+
    theme(legend.position="none")+
    ggtitle(title)+
    annotate("text", x = x.mean, y = y.mean,
             label = levels(plot_proj_df$color),
             size =max(p_size*2,3) )+
    # guides(colour = guide_legend(title="Annotations",
    #                              override.aes = list(size=3,alpha=1),
    #                              nrow=2))+
    ylab("UMAP 2")+ xlab("UMAP 1")
  # print(p2)
  p3<-p2
  if(!is.null(type)){
    specific = plot_proj_df[which(plot_proj_df$color == type),]
    specific = specific[stats::complete.cases(specific),]
    p2 <-p+geom_point(data = plot_proj_df,aes_string(x ='Y1',y = 'Y2'),size=p_size, color="gray80", alpha = 0.6,shape=1)

    p3<-p2+geom_point(data = specific, aes_string(x = 'Y1', y = 'Y2', color = 'batch'),
                       alpha = 0.9, size=p_size, shape=19, stroke = 0.00) +
      scale_color_manual(values = c("orange","dodgerblue","forestgreen"))+
      theme_classic()+
      theme(legend.position="bottom",
            legend.margin=margin(c(1,1,5,1)))+
      ggtitle(paste(title,type))+
      # annotate("text", x = x.mean, y = y.mean,
      #          label = levels(plot_proj_df$color),
      #          size =max(p_size*2,3) )+
      guides(colour = guide_legend(title="Batch",
                                   override.aes = list(size=2,alpha=1),
                                   nrow=1))+
      ylab("DIM 2")+ xlab("DIM 1")
  }

  # print(p3)

  if(!is.na(filename)){
    grDevices::pdf(filename,width = 6,height = 5)
    print(p3)
    grDevices::dev.off()
    return()
  }
  return(p3)


}






# --------------------------------------------------
# Violin Plot for marker genes
# --------------------------------------------------
#' Violin Plot
#' @description Produces violin plots for selected marker genes across predicted groups pf cells
#' @param de_data list containing (1) matrix subset, each row corresponds to a transcriptome sample, columns represent genes; (2) predicted labels of the samples in matrix
#' @param gene.list vector of specific marker genes which is a subset of \code{data} column names.
#' @param filename character specifying file location to save plot.
#' @return list of ggplot objects
#' @export
#' @importFrom plyr count
#' @importFrom ggplot2 ggplot aes geom_violin ylab xlab guides scale_fill_manual theme_light
plot_markers<-function(de_data, gene.list, filename=NA){


  data = de_data$mat_samples
  cluster.id = de_data$labels

  gene.list = unlist(gene.list)

  if(is.null(colnames(data))) stop("column names not found in data.")

  if(any(is.na(match(gene.list , colnames(data)))))
    stop("gene name(s) not found among genes in data header row.")

  n.colors <-length(unique(cluster.id))
  mk.plots=list()
  for(i in gene.list ){
    per_marker = data.frame(cbind(count=log2(data[,i]+1),
                                  cluster.id = cluster.id))
    per_marker$cluster.id<-as.factor(per_marker$cluster.id)
    p <- ggplot(per_marker, aes(x=cluster.id, y=count,fill = cluster.id))
    p2<-p + geom_violin(trim = FALSE, scale = "width", draw_quantiles = TRUE)+
      ylab(i)+xlab("Cluster ID")+ guides(fill=FALSE)+
      scale_fill_manual(values = getColors(n.colors))+
      theme_light()
    mk.plots[[i]]<-p2
  }

  n = length(gene.list)
  plot.nrows = ceiling(n / 3)
  plot.ncols = ifelse(n>=3,3,n)

  if(!is.na(filename)){
    grDevices::pdf("marker_plot.pdf",height=3*plot.nrows,width=10)
  gridExtra::grid.arrange(grobs=mk.plots,ncol=plot.ncols,nrow=plot.nrows)
  grDevices::dev.off()
  } else{
    gridExtra::grid.arrange(grobs=mk.plots,ncol=plot.ncols,nrow=plot.nrows)
  }


  return(mk.plots)
}


# # first.batch <- rep(c(TRUE, FALSE), c(ncol(logDataF3), ncol(logDataA3)))
#
# # Adding colours.
# #base.color <- "grey"
# color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta",
#                    HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", Unsorted="grey")
# colmatF <- col2rgb(color.legendF)
#
# colmatA <- colmatF + 100 # A lighter shade.
# colmatA[colmatA > 255] <- 255
#
# #colmatF<-colmatF+200
# #colmatF[colmatF > 255] <- 255
# color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
# allcolors <- c(color.legendF[colnames(logDataF3)], color.legendA[colnames(logDataA3)])
# batch<-c( rep(1,ncol(logDataF3)),rep(2,ncol(logDataA3)) )
#
# # Only keeping common cell types for PCA.
# celltypes <- c(colnames(logDataF3), colnames(logDataA3))
# pca.retain <- celltypes %in% c("MEP", "GMP", "CMP")
#
# # Making a plotting function.
# plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2",main="") {
#   if (is.null(subset)) {
#     subset <- seq_len(nrow(Y))
#   }
#   png(fname,width=900,height=700)
#   par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
#   plot(Y[,1], Y[,2], cex=2,
#        pch=ifelse(first.batch, 21, 1)[subset],
#        col=ifelse(first.batch, "black", allcolors)[subset],
#        bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main)
#   dev.off()
# }
#
# batchcolor=c("lavender","lightcoral")
# plotFUNb <- function(fname, Y, subset=NULL, ...) {
#   if (is.null(subset)) {
#     subset <- seq_len(nrow(Y))
#   }
#   png(fname,width=900,height=700)
#   par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
#   plot(Y[,1], Y[,2], cex=2,
#        pch=ifelse(first.batch, 21, 1)[subset],
#        col=ifelse(first.batch, "black", batchcolor[batch[subset]]),
#        bg=batchcolor[batch[subset]], ...)#,  xlab="tSNE 1",ylab="tSNE 2")
#   dev.off()
# }
#

