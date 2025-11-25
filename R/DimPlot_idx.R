gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' DimPlot with index prefix.
#' 
#' @examples 
#' p1 <- DimPlot(pbmc)
#' p2 <- DimPlot_idx(pbmc)
#' p1 + p2
#' 
#' @export
DimPlot_idx <- function(object,
                        group.by = NULL,
                        reduction = 'umap',
                        label.idx = T,
                        label.idx.size = 5,
                        label.size.scale = 1,
                        label.idx.color = 'black',
                        label.idx.color.back = 'white',
                        label.idx.fontface = 'bold',
                        repel = F,
                        idx.sep = ' - ',
                        idx.in.key = F,
                        prefix.index = NULL,
                        cover.outlier = F, 
                        alter.group.for.outlier = 'seurat_clusters', 
                        rev.outlier.order = F,
                        legend.ncol = 2,
                        legend.byrow = F,
                        legend.title = NULL,
                        legend.position = 'right',
                        legend.key.size = 5,
                        legend.text.size = 8,
                        shorter.axis = F,
                        axis.len = 0.2,
                        axis.adj.x = 0,
                        axis.adj.y = 0,
                        axis.title.size = 12,
                        arrow.length = unit(.3, 'cm'),
                        arrow.type = "closed",
                        arrow.size = 3,
                        AI.friendly = F,
                        return.obj = F,
                        combine = F,
                        cols = NULL,
                        ...) {
  if (is.null(group.by)) {
    object$annotation <- Idents(object)
  }else{
    object$annotation <- object@meta.data[,group.by]
  }
  if (!is.factor(object$annotation)) {
    object$annotation <- factor(object$annotation, levels = unique(object$annotation))
  }
  # add index for annotation
  numcl <- length(levels(object$annotation))
  if (is.null(prefix.index)) {
    prefix.index <- 1:numcl
  }
  if (is.null(cols)) cols <- gg_color_hue(numcl)
  if (!is.null(names(cols))) {
    cols <- cols[names(cols) %in% levels(object$annotation)]
    names(cols) <- paste0(prefix.index,
                          idx.sep,
                          names(cols))
  }
  annodf <- data.frame('annotation' = levels(object$annotation),
                       'annotation.idx' = paste0(prefix.index,
                                          idx.sep,
                                          levels(object$annotation)),
                       'idx' = prefix.index,
                       stringsAsFactors = F)
  meta <- object@meta.data
  meta <- left_join(meta, annodf, by = 'annotation')
  meta$annotation.idx <- factor(meta$annotation.idx, levels = annodf$annotation.idx)
  rownames(meta) <- rownames(object@meta.data)
  object@meta.data <- meta
  
  # plot
  p <- DimPlot(object, group.by = 'annotation.idx', cols = cols, ...) + ggtitle('')
  if (cover.outlier) {
    point_order <- cover_outlier(object, 'annotation.idx', reduction = reduction,
                                 seurat_clusters = alter.group.for.outlier)
    if (rev.outlier.order) point_order <- rev(point_order)
    p$data <- p$data[rev(point_order),]
  }
  if (AI.friendly) p <- ggAIplot(p)
  
  meta_n <- object@meta.data[,c('idx'), drop = F]
  meta_n <- meta_n[rownames(p$data),]
  p$data$cluster_n <- meta_n
  #convert cluster_n as factor to avoid error:
  # Error in LabelClusters(p, id = "cluster_n", size = label.idx.size) : 
  #  Length of labels (9) must be equal to the number of clusters being labeled (4).
  p$data$cluster_n <- factor(p$data$cluster_n, 
                             levels = unique(prefix.index))
  if (label.idx) {
    p <- LabelClusters(p, id = 'cluster_n', size = label.idx.size, 
                       fontface = label.idx.fontface,
                       repel = repel, color = label.idx.color.back)
    p <- LabelClusters(p, id = 'cluster_n', size = label.idx.size * label.size.scale,
                       color = label.idx.color, repel = repel)
    
  }
  if (shorter.axis) {
    data <- p$data
    xmin <- min(data[,1]) + axis.adj.x
    ymin <- min(data[,2]) + axis.adj.y
    axis.x.len <- (max(data[,1]) - min(data[,1])) * axis.len + min(data[,1])
    axis.y.len <- (max(data[,2]) - min(data[,2])) * axis.len + min(data[,2])
    p <- p + annotate(geom="segment",x=xmin, y=ymin, xend=axis.x.len, yend=ymin, 
                          arrow = arrow(length = arrow.length, type = arrow.type))+
      annotate(geom="segment",x=xmin, y=ymin, xend=xmin, yend=axis.y.len,
                   arrow = arrow(length = arrow.length, type = arrow.type))+
      labs(x = colnames(data)[1], y =colnames(data)[2])+
      theme_void()+
      theme(axis.title.x = element_text(hjust = 0.05, vjust = 4, size = axis.title.size),
            axis.title.y = element_text(hjust = 0.05, vjust = -2, angle=90, size = axis.title.size)
      )
  }
  if (idx.in.key) {
    data <- p$data
    data$label <- data$cluster_n
    
    numcl <- length(unique(data$label))
    keylabeldf <- data.frame(label = 1:numcl)
    keylabeldf$colour <- cols[1:numcl]
    keylabeldf <<- keylabeldf
    data$annotation <- gsub(paste0('^.*',idx.sep),'',as.vector(data$annotation.idx))
    data$annotation <- factor(data$annotation, levels = gsub(paste0('^.*',idx.sep),'',levels(data$annotation.idx)))
    draw_key_cust <- function(data, params, size) {
      data_text <- data
      data_text[c("fill")] <- NULL
      ## NOTE: colours for labels should be unique
      data_text$label <- as.vector(keylabeldf[keylabeldf$colour == data_text$colour,'label'])
      data_text$colour <- "white"
      #data_text$colour <- NULL
      data_text$alpha <- 1
      data_text$size <- 8 / .pt
      grid::grobTree(
        draw_key_point(data, list()),
        draw_key_text(data_text, list())
      )
    }
    pl <- ggplot(data, aes(data[,1], data[,2], colour = annotation)) + 
      geom_point(key_glyph = draw_key_cust) + theme_bw() +
      scale_colour_manual(values = cols) + 
      theme(legend.text = element_text(size = legend.text.size)) +
      guides(colour=guide_legend(title = legend.title,
                                 ncol= legend.ncol, byrow=legend.byrow, position = 'right', 
                                 override.aes = list(size = legend.key.size )))
    legend <- cowplot::get_plot_component(pl, 'guide-box-right', return_all = TRUE)
    if (legend.position == 'bottom') {
      p <- patchwork::wrap_plots(p + NoLegend(), cowplot::ggdraw(legend), ncol = 1)
    }else if (legend.position == 'top') {
      p <- patchwork::wrap_plots(cowplot::ggdraw(legend),p + NoLegend(), ncol = 1)
    }else if  (legend.position == 'left') {
      p <- patchwork::wrap_plots(cowplot::ggdraw(legend),p + NoLegend(), ncol = 2)
    }else {
      p <- patchwork::wrap_plots(p + NoLegend(),cowplot::ggdraw(legend), ncol = 2)
    }
  }
  if (return.obj) {
    return(object)
  }else{
    return(p)
  }
}

