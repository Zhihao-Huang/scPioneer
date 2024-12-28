#' highlight clusters
#' 
#' @export
ArchR_highallcl <- function(object, colorBy = "cellColData", name = "Clusters",
                            size = 0.1,size.title = 10,
                            embedding = "UMAP", ncol = NULL) {
  plist <- list()
  for (i in unique(object@cellColData$Clusters)) {
    cellID <- getCellNames(object)[object@cellColData$Clusters == i]
    p <- plotEmbedding(ArchRProj = object, colorBy = colorBy, name = name,
                       embedding = embedding, highlightCells = cellID, 
                       keepAxis = F)
    p <- p + ggtitle(i) + theme(legend.position = 'None',
                                plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.text.y=element_blank(),
                                axis.ticks.y=element_blank(),
                                plot.title = element_text(hjust = 0.5,
                                                          size = size.title)) +
    xlab('') + ylab('') 
    p$data <- p$data[order(p$data$color, decreasing = T),]
    plist[[i]] <- p
  }
  do.call(cowplot::plot_grid, c(list(ncol = ncol),plist))
}
#' highlight cluster
#' 
#' @export
ArchR_one_cl <- function(object, cl, colorBy = "cellColData", name = "Clusters",
                            size = 0.1,size.title = 10,
                            embedding = "UMAP") {
  cellID <- getCellNames(object)[object@cellColData$Clusters == cl]
  p <- plotEmbedding(ArchRProj = object, colorBy = colorBy, name = name,
                     embedding = embedding, highlightCells = cellID, 
                     keepAxis = F)
  p <- p + ggtitle(cl) + theme(legend.position = 'None',
                               plot.margin = unit(c(0, 0, 0, 0), "cm"),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank(),
                               plot.title = element_text(hjust = 0.5,
                                                         size = size.title)) +
    xlab('') + ylab('') 
  p$data <- p$data[order(p$data$color, decreasing = T),]
  p
}

#' plot embedding
#' 
#' @export
plotEmbedding2 <- function(ArchRProj = NULL, embedding = "UMAP",
                           name = "Sample",
                           AI.friendly = F,
                           AI.script = '/BGFS1/projectdata/jasper/project/ovarian_cancer/shell/functions/ggAIplot.R',
                           label.name = NULL,
                           label.size = 4,
                           color = NULL,
                           cell.order = NULL,
                           pt.size = 0.1,
                           legend.position = 'bottom',
                           legend.text.size = 10,
                           legend.point.size = 4,
                           legend.nrow = 3) {
  if (!name %in% colnames(ArchRProj@cellColData)) {
    stop(paste0(name, ' is not in cellColData.'))
  }
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  df$Cluster <- as.vector(ArchRProj@cellColData[,name])
  colnames(df) <- c("UMAP_1",'UMAP_2','Cluster')
  if (is.null(color)) {
    color <- paletteDiscrete(values = df$Cluster)
  }
  if (!is.null(cell.order)) {
    df$Cluster <- factor(df$Cluster, levels = cell.order)
    names(color) <- cell.order
  }else{
    df$Cluster <- factor(df$Cluster, levels = unique(df$Cluster))
  }
  p <- ggplot(data = df, aes(UMAP_1, UMAP_2)) + 
    geom_point(aes(color = Cluster),
               size = pt.size) + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = legend.text.size), 
          legend.position = legend.position,
          legend.title = element_blank()) + 
    scale_color_manual(values = color) +
    guides(color=guide_legend(nrow=legend.nrow,byrow=TRUE,
                              override.aes = list(size=legend.point.size)))
  if (AI.friendly) {
    source(AI.script)
    p <- ggAIplot(p)
  }
  if (!is.null(label.name)) {
    if (!label.name %in% colnames(ArchRProj@cellColData)) {
      stop(paste0(label.name, ' is not in cellColData.'))
    }
    df <- data.frame(label_text = as.vector(ArchRProj@cellColData[,label.name]), 
                     stringsAsFactors = F)
    rownames(df) <- rownames(ArchRProj@cellColData)
    df <- df[rownames(p$data),,drop = F]
    if (is.null(levels(df$label_text))) {
      df$label_text <- factor(df$label_text, levels = as.character(unique(df$label_text)))
    }
    p$data$label_text <- df$label_text
    p <- Seurat::LabelClusters(p, id = 'label_text', size = label.size, color = 'black')
    
  }
  p
}
