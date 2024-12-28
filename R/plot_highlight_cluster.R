#' highlight points of one cluster
#'
#' @param SeuratS4 A Seurat object.
#' @param cluster Select which cluster to highlight.
#' @param reduction Select the reduction method.
#' @param top.layer If TRUE, put the cluster points on the top layer of plot.
#' @param label Parameter in Dimplot.
#' @param label.size Parameter in Dimplot.
#' @param pt.size Parameter in Dimplot.Set the point size.
#' @param ... Parameter in Dimplot.
#' @return
#' A ggplot objects.
#' @examples
#' data('pbmc')
#' highcl(pbmc,'Naive CD4 T',reduction = "umap")
#' @export

highcl <- function(SeuratS4, cluster, reduction = "tsne", color.1 = 'red',
                   color.2 = 'grey', top.layer = T, label = FALSE, label.size = 3,
                   pt.size = 0.1, ...) {
  clusters <- levels(Idents(SeuratS4))
  cl_num <- length(clusters)
  pos <- which(clusters %in% cluster)
  if (top.layer) {
    order <- rev(c(clusters[-pos], cluster))
    color_palette <- c(rep("grey", cl_num - 1),"red")
  }else{
    order <- rev(clusters)
    color_palette <- c(rep("grey", pos - 1), "red", rep("grey", cl_num - pos))
  }
  p <- DimPlot(SeuratS4, reduction = reduction, label = label, order = order,
               label.size = label.size, cols = color_palette, pt.size=pt.size,...)
  return(p)
}

#' highlight points of all cluster
#'
#' @param SeuratS4 A Seurat object.
#' @param reduction Celect the reduction method.
#' @return
#' A ggplot objects.
#' @examples
#' data('pbmc')
#' highallcl(pbmc,reduction = "umap")
#' @export

highallcl <- function(SeuratS4 = SeuratS4, reductions = "tsne", ncol = NULL,
                      color = c("red", "grey")) {
  xy <- SeuratS4@reductions[[reductions]]@cell.embeddings
  lablist <- list()
  clusters <- levels(Idents(SeuratS4))
  clusters <- gsub(" ", "_", clusters)
  ident <- as.vector(Idents(SeuratS4))
  ident <- gsub(" ", "_", ident)
  for (cl in clusters) {
    pos <- ident %in% cl
    newident <- ident
    newident[pos] = color[1]
    newident[!pos] = color[2]
    mat <- data.frame(xy, newident, cl)
    lablist[[paste0("cluster_", cl)]] <- mat
  }
  labmat <- do.call(rbind, lablist)
  labmat$newident <- factor(labmat$newident, levels = color)
  labmat <- plyr::arrange(labmat, plyr::desc(newident), cl)
  p <- ggplot(labmat, aes(x = labmat[, 1], y = labmat[, 2], colour = labmat$newident)) +
    geom_point(size = 0.1) + theme_bw() + theme(panel.grid.major = element_line(colour = NA),
                                                panel.grid.minor = element_line(colour = NA)) + ylab(reductions) +
    xlab(reductions) + scale_color_manual(values = color)
  p <- p + facet_wrap(~cl, ncol = ncol)
  p <- p + theme_void() + theme(legend.position = "none")
  print(p)
}


#' highlight points of all cluster by separated UMAP in different colors.
#'
#' @param SeuratS4 A Seurat object.
#' @param reduction Celect the reduction method.
#' @return
#' A ggplot objects.
#' @examples
#' data('pbmc')
#' highallcl(pbmc,reduction = "umap")
#' @export

highallcl2 <- function(SeuratS4, group.by = NULL, 
                       reductions = "umap", pt.size = 0.1,
                       ncol = NULL, color = NULL) {
  xy <- SeuratS4@reductions[[reductions]]@cell.embeddings
  lablist <- list()
  if (is.null(group.by)) {
    clusters <- levels(Idents(SeuratS4))
    clusters <- gsub(" ", "_", clusters)
    ident <- as.vector(Idents(SeuratS4))
    ident <- gsub(" ", "_", ident)
  }else{
    clusters <- as.vector(sort(unique(SeuratS4@meta.data[,group.by])))
    clusters <- gsub(" ", "_", clusters)
    ident <- as.vector(SeuratS4@meta.data[,group.by])
    ident <- gsub(" ", "_", ident)
  }
  
  if (is.null(color)) {
    color = c("red", "lightgrey")
    for (cl in clusters) {
      pos <- ident %in% cl
      newident <- ident
      newident[pos] = color[1]
      newident[!pos] = color[2]
      mat <- data.frame(xy, newident, cl)
      lablist[[paste0("cluster_", cl)]] <- mat
    }
  }else{
    for (i in 1:length(clusters)) {
      cl <- clusters[i]
      pos <- ident == cl
      newident <- ident
      newident[pos] = color[i]
      newident[!pos] = 'lightgrey'
      mat <- data.frame(xy, newident, cl)
      lablist[[paste0("cluster_", cl)]] <- mat
    }
    #color_seq <- rep('lightgrey',length(clusters) * 2)
    #color_seq[seq(1,length(color_seq),by = 2)] <- color[1:length(clusters)]
    #color <- color_seq
    color <- c(color,'lightgrey')
  }
  labmat <- do.call(rbind, lablist)

  labmat$newident <- factor(labmat$newident, levels = unique(color))
  labmat <- plyr::arrange(labmat, plyr::desc(newident), cl)
  p <- ggplot(labmat, aes(x = labmat[, 1], y = labmat[, 2], 
                          color = newident)) +
    geom_point(size = pt.size) + theme_bw() + 
    theme(panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA)) + scale_color_manual(values = color)
  p <- p + facet_wrap(~cl, ncol = ncol)
  p <- p + theme_void() + theme(legend.position = "none")
  return(p)
}

#' highlight points of all cluster by separated UMAP in different colors.
#'
#' @param SeuratS4 A Seurat object.
#' @param reduction Celect the reduction method.
#' @return
#' A ggplot objects.
#' @examples
#' data('pbmc')
#' highallcl(pbmc,reduction = "umap")
#' @export

highallcl3 <- function(SeuratS4, group.by = NULL, 
                       reductions = "umap",
                       text.size = 8, pt.size = 0.1, AI.friendly = F,
                       return.list = F,
                       ncol = NULL, color = NULL) {
  xy <- SeuratS4@reductions[[reductions]]@cell.embeddings[,1:2]
  if (is.null(group.by)) {
    clusters <- levels(Idents(SeuratS4))
    ident <- as.vector(Idents(SeuratS4))
  }else{
    clusters <- as.vector(sort(unique(SeuratS4@meta.data[,group.by])))
    ident <- as.vector(SeuratS4@meta.data[,group.by])
  }
  data <- data.frame(x = xy[,1], y = xy[,2], ident, stringsAsFactors = F)
  
  if (is.null(color)) {
    color = c("red", "lightgrey")
    color <- rep('red',length(clusters))
  }

  plist <- list()
  for (i in 1: length(clusters)) {
    plot.data <- data
    plot.data$ident <- as.vector(plot.data$ident)
    plot.data$ident[plot.data$ident != clusters[i]] <- 'Others'
    point.order <- rev(unique(c('Others',unique(plot.data$ident))))
    plot.data$ident <- factor(plot.data$ident, 
                                 levels = point.order)
    plot.data <- plyr::arrange(plot.data, plyr::desc(ident))
    color.use <- c(color[i], 'lightgrey')
    names(color.use) <- point.order
    p <- ggplot(plot.data, aes(x , y, 
                            color = ident)) +
      geom_point(size = pt.size) + ggtitle(clusters[i]) + theme_bw() + 
      theme(panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) + 
      scale_color_manual(values = color.use)
    p <- p + theme_void() + theme(legend.position = "none", 
                                  plot.title = element_text(hjust = 0.5, size = text.size))
    plist[[clusters[i]]] <- p
  }
  if (AI.friendly) {
    plist <- lapply(plist, function(x) ggAIplot(x))
    names(plist) <- clusters
  }
  if (return.list) {
    return(plist)
  }else{
    pall <- patchwork::wrap_plots(plist, ncol = ncol)
    return(pall)
  }
  
}
