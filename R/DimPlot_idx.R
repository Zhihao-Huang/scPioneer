gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
draw_key_cust <- function(data, params, size) {
  data_text <- data
  data_text[c("fill")] <- NULL
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
  if (is.null(prefix.index)) {
    prefix.index <- 1:length(levels(object$annotation))
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
    p <- LabelClusters(p, id = 'cluster_n', size = label.idx.size, repel = repel)
  }
  if (idx.in.key) {
    data <- p$data
    data$label <- data$cluster_n
    
    numcl <- length(unique(data$label))
    keylabeldf <- data.frame(label = 1:numcl)
    if (is.null(cols)) cols <- gg_color_hue(numcl)
    keylabeldf$colour <- cols[1:numcl]
    keylabeldf <<- keylabeldf
    data$annotation <- gsub(paste0('^.*',idx.sep),'',as.vector(data$annotation.idx))
    data$annotation <- factor(data$annotation, levels = gsub(paste0('^.*',idx.sep),'',levels(data$annotation.idx)))
    
    pl <- ggplot(data, aes(data[,1], data[,2], colour = annotation)) + 
      geom_point(key_glyph = "cust") + theme_bw() +
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

