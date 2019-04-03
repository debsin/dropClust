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
all_plot<-function(data,filename=NA, title){

  temp = stats::complete.cases(data)
  plot_proj_df = data[temp, ]

  plot_proj_df$color<-factor(plot_proj_df$color)

  x.mean = stats::aggregate(plot_proj_df$Y1,
                            list(plot_proj_df$color),
                            stats::median)[,-1]
  y.mean = stats::aggregate(plot_proj_df$Y2,
                            list(plot_proj_df$color),
                            stats::median)[,-1]

  colorcount_t = length(levels(plot_proj_df$color))


  n_points = dim(plot_proj_df[stats::complete.cases(plot_proj_df),])[1]

  p_size = -0.58*log(n_points)+6.77
  p<-ggplot(data = plot_proj_df)
  p2<-p+ geom_point(aes_string(x ='Y1',y = 'Y2',col= 'color'),size=p_size)  +
    scale_colour_manual(values =  getColors(colorcount_t))+
    theme_classic()+
    theme(legend.position="bottom")+
    ggtitle(title)+
    annotate("text", x = x.mean, y = y.mean,
             label = levels(plot_proj_df$color),
             size =max(p_size*2,3) )+
    guides(colour = guide_legend(title="Cluster IDs",
                                 override.aes = list(size=3,alpha=1),
                                 nrow=2))+
                ylab("t-SNE 2")+ xlab("t-SNE 1")

  if(!is.na(filename)){
    grDevices::pdf(filename,width = 6,height = 5)
    print(p2)
    grDevices::dev.off()
    return()
  }
  return(p2)


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




heat_labels <- function(labels)
{
  ordered_labl = labels
  colors = getColors(length(unique(labels)))
  colors = colors[ordered_labl]
  ordered_freq = table(factor(ordered_labl,levels=unique(ordered_labl)))
  a = cumsum(ordered_freq)
  b = c(0,  a[-length(a)])
  pos = round((a+b)/2)
  text_lab = rep("",length(labels))
  text_lab[pos] = unique(ordered_labl);

  return(list("colors"=colors, "text_lab" = text_lab))

}


# --------------------------------------------------
# Return heatmap object
# --------------------------------------------------
#' Build heatmap object
#' @description Produces custom heatmap uses \code{pheatmap::pheatmap()}.
#' @param de_data list containing (1) matrix subset, each row corresponds to a transcriptome sample, columns represent genes; (2) predicted labels of the samples in matrix
#' @param DE_res list obtained from the \code{DE_genes} module.
#' @param selected_clusters vector of selected clusters to be considered. When unspecified (=NA), defaults to all predicted clusters.
#' @param nDE integer, specifies the number of DE genes to appear in the heatmap.
#' @return heatmap object.
#' @export
plot_heatmap<-function(de_data, DE_res, selected_clusters=NA, nDE=10){

  data = de_data$mat_samples
  cluster.id = de_data$labels




  if(any(is.na(selected_clusters))==TRUE)
    selected_clusters = unique(cluster.id)


  if(!all(selected_clusters %in% cluster.id))
    stop(paste("Cluster Ids must be among:", paste(unique(cluster.id),collapse = " ")))

  cell.ids = which(cluster.id %in% selected_clusters)


  de_genes<-top.de.genes(DE_res, nDE)

  all_ct_genes = unique(unlist(de_genes))

  heat_in = log2(t(as.matrix(data[cell.ids, all_ct_genes]))+1)


  heat_param<-heat_labels(labels = cluster.id[cell.ids])
  myPalette <- grDevices::colorRampPalette(c("red",
                                             "darkorange4",
                                             "gray15",
                                             "black",
                                             "yellow4",
                                             "yellow",
                                             "greenyellow"))


  # Generate annotations for rows and columns
  annotation_col = data.frame(cbind(CellType = factor(cluster.id[cell.ids])))
  rownames(annotation_col) = colnames(heat_in)

  custom.colors = getColors(length(unique(cluster.id[cell.ids])))
  names(custom.colors) =  unique(cluster.id[cell.ids])

  ann_colors = list(CellType = custom.colors)


  p<-pheatmap::pheatmap(mat = heat_in,scale = "row",cluster_cols = FALSE,
                        cluster_rows = FALSE,
                        labels_col = heat_param$text_lab,
                        annotation_legend = FALSE,
                        annotation_col = annotation_col,
                        annotation_colors = ann_colors,
                        legend=TRUE,
                        fontsize_row = 4,
                        color= myPalette(100),
                        silent = TRUE)
  return(p)

}
