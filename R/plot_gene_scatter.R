#' Modified FeaturePlot, to reduce legend size, and add more colors.
#' 
#' @param legend.position Vector with 2 elements. Adjust position of legend by location of x/y axis.
#' 
#' @examples 
#' FeaturePlot2(object = pbmc, features = c('CD3D','CD79A','CD14','FCGR3A'))
#' FeaturePlot2(object = pbmc, features = c('CD3D','CD79A','CD14','FCGR3A'),
#' legend.position = c(0.9,0.2))
#' 
#' @export
FeaturePlot2 <- function (object, features, color = NULL, pt.size = 0.1, 
          legend.position = c(1,0.2),
          legend.key.size = 0.3, legend.text.size = 6,
          ncols = NULL, split.by = NULL,
          plot.margin = c(0, 1, 0, 0), AI.friendly = F, dpi = 300,
          AI.script = "/BGFS1/projectdata/jasper/project/ovarian_cancer/shell/functions/ggAIplot.R", 
          ...) {
  if (is.null(color)) {
    color <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
    color <- (grDevices::colorRampPalette(c("grey", color)))(5)
  }
  feature_theme <- theme(legend.text = element_text(size = legend.text.size), 
                         legend.key.size = unit(legend.key.size, "cm"), legend.position = legend.position, 
                         axis.line = element_blank(), axis.text = element_blank(), 
                         axis.ticks = element_blank(), plot.margin = unit(plot.margin, 
                                                                          "cm"))
  if (ncol(object) > 100000 | ncol(object) < 10000 | AI.friendly) pt.size = 1
  if (!is.null(split.by)) {
    p <- FeaturePlot(object = object, features = features, pt.size = pt.size, 
                     split.by = split.by,...)
    p <- p & feature_theme & xlab("") & ylab("")
    p <- suppressMessages(p & scale_color_gradientn(colors = color))
    
  }else{
    plist <- lapply(features, function(x) {
      p <- FeaturePlot(object = object, features = x, pt.size = pt.size, 
                       ...)
      p <- suppressMessages(p + scale_color_gradientn(colors = color))
      p <- p + feature_theme + xlab("") + ylab("") + ggtitle(x)
      if (AI.friendly) {
        source(AI.script)
        p <- ggAIplot(p, dpi = dpi)
      }
      p
    })
    p <- patchwork::wrap_plots(plist,ncol = ncols)
  }
  return(p)
}

#' Scater plot of gene expression.
#'
#' @return A ggplot object.
#' @param exprmat A normalize selected genes expression dataframe(cell x gene).
#' @param xy Cell location from seurat_object@@reduction$tsne@@cell.embeddings.
#' 
#' @examples
#' data('pbmc')
#' marker_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
#' exp <- log2(as.data.frame(t(as.matrix(pbmc@@assays$RNA@@data[marker_to_plot,])))+1)
#' xy <- pbmc@@reductions$umap@@cell.embeddings
#' fplot1(xy,exp)
#' 
#' @export
fplot1 <- function (xy, exprmat, ncol = NULL, legend.title = "auto", pt.size = 0.1, 
          color.use = c("grey", "yellow", "red"), order.expression = "ascending", 
          AI.friendly = F, 
          AI.script = "/BGFS1/projectdata/jasper/project/ovarian_cancer/shell/functions/ggAIplot.R") 
{
  xymat <- cbind(xy, exprmat)
  xynames <- colnames(xy)
  plot.data <- melt(xymat, id.vars = xynames, measure.vars = colnames(exprmat))
  colnames(plot.data) <- c(xynames, "markers", "Expression")
  print("Plot genes: ")
  print(colnames(exprmat))
  if (is.null(order.expression)) {
    plot.data <- plot.data
  }
  else if (order.expression == "ascending") {
    plot.data <- arrange(plot.data, Expression, markers)
  }
  else if (order.expression == "descending") {
    plot.data <- arrange(plot.data, markers, desc(Expression))
  }
  colnames(plot.data)[1:2] <- c("x", "y")
  p <- ggplot(plot.data, aes(x = x, y = y, color = Expression)) + 
    geom_point(size = pt.size)
  if (legend.title != "auto") {
    p <- p + scale_colour_gradientn(colours = color.use, 
                                    name = legend.title)
  }
  else {
    p <- p + scale_colour_gradientn(colours = color.use)
  }
  p <- p + facet_wrap(~markers, ncol = ncol)
  p <- p + theme_void()
  if (AI.friendly) {
    source(AI.script)
    p <- ggAIplot.grid(p, facet.by = "markers")
  }
  return(p)
}


#' Scater plot of gene expression with label.
#'
#' @return A ggplot object.
#' @param exprmat A normalize selected genes expression dataframe(cell x gene).
#' @param xy Cell location from seurat_object@@reduction$tsne@@cell.embeddings.
#' @param label.cluster A vector stored cluster labels in meta.data.
#' 
#' @examples
#' data('pbmc')
#' marker_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
#' exp <- log2(as.data.frame(t(as.matrix(pbmc@@assays$RNA@@data[marker_to_plot,])))+1)
#' xy <- pbmc@@reductions$umap@@cell.embeddings
#' label.cluster <- pbmc$Annotation
#' fplot_label(xy,exp,label.cluster)
#' 
#' @export
fplot_label <- function (xy, exprmat, label.cluster, ncol = NULL, legend.title = "auto", 
          reduction = "tsne", label.size = 4, fontface = "plain", pt.size = 0.1, 
          AI.friendly = F, 
          AI.script = "/BGFS1/projectdata/jasper/project/ovarian_cancer/shell/functions/ggAIplot.R", 
          color.use = c("grey", "yellow", "red"), order.expression = "ascending") 
{
  xymat <- cbind(xy, exprmat)
  xynames <- colnames(xy)
  plot.data <- melt(xymat, id.vars = xynames, measure.vars = colnames(exprmat))
  colnames(plot.data) <- c(xynames, "markers", "Expression")
  print("Plot genes: ")
  print(colnames(exprmat))
  if (is.null(order.expression)) {
    plot.data <- plot.data
  }
  else if (order.expression == "ascending") {
    plot.data <- arrange(plot.data, Expression, markers)
  }
  else if (order.expression == "descending") {
    plot.data <- arrange(plot.data, markers, desc(Expression))
  }
  plot.data$cluster_n <- rep(label.cluster, length(unique(plot.data$markers)))
  groups <- unique(label.cluster)
  id <- "cluster_n"
  xynames <- c("x", "y")
  colnames(plot.data)[1:2] <- xynames
  labels.loc <- lapply(X = groups, FUN = function(group) {
    data.use <- plot.data[plot.data[, id] == group, , drop = FALSE]
    data.medians <- as.data.frame(x = t(x = apply(X = data.use[, xynames, drop = FALSE], 
                                                  MARGIN = 2, FUN = median, 
                                                  na.rm = TRUE)))
    data.medians[, id] <- group
    return(data.medians)
  })
  labels.loc <- do.call(rbind, labels.loc)
  p <- ggplot(data = plot.data, aes(x = x, y = y, colour = Expression)) + 
    geom_point(size = pt.size)
  if (legend.title != "auto") {
    p <- p + scale_colour_gradientn(colours = color.use, 
                                    name = legend.title)
  }
  else {
    p <- p + scale_colour_gradientn(colours = color.use)
  }
  p <- p + facet_wrap(~markers, ncol = ncol)
  p <- p + theme_void()
  if (AI.friendly) {
    source(AI.script)
    p <- ggAIplot.grid(p, facet.by = "markers")
  }
  p <- p + geom_text(data = labels.loc, aes(x = x, y = y, label = cluster_n), 
                     size = label.size, fontface = fontface, inherit.aes = FALSE)
  return(p)
}
#' Scater plot of gene expression.
#'
#' @return A ggplot object.
#' 
#' @examples
#' marker_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
#' fplot2(pbmc,marker_to_plot)
#' 
#' @export
fplot2 <- function (object, features, reduction = "umap", ncol = NULL, 
          assay = "RNA", label.cluster = NULL, label.size = 4, fontface = "plain", 
          AI.friendly = F,
          AI.script = "/BGFS1/projectdata/jasper/project/ovarian_cancer/shell/functions/ggAIplot.R", 
          legend.title = "auto", pt.size = 0.1, color.use = c("grey", "yellow", "red"),
          order.expression = NULL) 
{
  pos <- features %in% rownames(object)
  if (sum(!pos) > 0) {
    message(paste0("Warnning: ", paste(features[!pos], collapse = ","), 
                   " not in data."))
  }
  features <- features[pos]
  xy <- object@reductions[[reduction]]@cell.embeddings
  DefaultAssay(object = object) <- assay
  if (nrow(object@assays[[assay]]@data) == 0) {
    stop("Error in array selection.")
  }
  exp <- as.data.frame(t(as.matrix(object@assays[[assay]]@data[features, 
                                                               , drop = F])))
  if (!is.null(label.cluster)) {
    if (length(label.cluster) == 1) {
      if (!label.cluster %in% colnames(object@meta.data)) {
        stop("label.cluster is not a vector, nor is colnames in meta.data.")
      }
      label.cluster <- object@meta.data[, label.cluster]
    }
    p <- fplot_label(xy, exp, ncol = ncol, color.use = color.use, 
                     reduction = reduction, label.cluster = label.cluster, 
                     label.size = label.size, fontface = fontface, pt.size = pt.size, 
                     AI.friendly = AI.friendly, AI.script = AI.script, 
                     legend.title = legend.title, order.expression = order.expression)
  }
  else {
    p <- fplot1(xy, exp, ncol = ncol, color.use = color.use, 
                pt.size = pt.size, order.expression = order.expression, 
                legend.title = legend.title, AI.friendly = AI.friendly, 
                AI.script = AI.script)
  }
  return(p)
}

#' Scater plot of gene expression.
#'
#' One gene one bar.
#'
#' @return A ggplot object.
#' 
#' @examples
#' marker_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
#' fplot3(pbmc,marker_to_plot)
#' 
#' @export
fplot3 <- function (object, features, assay = "RNA", color = NULL, pt.size = 0.1, 
                    label.cluster = NULL, label.size = 3, order.expression = F, 
                    reduction = "umap", legend.title = "", legend.position = c(1,0.2),
                    legend.key.size = 0.3, legend.text.size = 6, strip.text.size = 13, 
                    plot.margin = c(0, 1, 0, 0), ncol = NULL, AI.friendly = F, 
                    AI.script = "/BGFS1/projectdata/jasper/project/ovarian_cancer/shell/functions/ggAIplot.R") 
{
  if (is.null(color)) {
    color <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
    color <- (grDevices::colorRampPalette(c("grey", color)))(5)
  }
  feature_theme <- theme(legend.text = element_text(size = legend.text.size), 
                         legend.key.size = unit(legend.key.size, "cm"), legend.position = legend.position, 
                         axis.line = element_blank(), axis.text = element_blank(), 
                         axis.ticks = element_blank(), strip.text = element_text(size = strip.text.size), 
                         plot.margin = unit(plot.margin, "cm"))
  plist <- lapply(features, function(x) {
    p <- fplot2(object = object, features = x, assay = assay, 
                pt.size = pt.size, label.cluster = label.cluster, 
                label.size = label.size, color.use = color, order.expression = order.expression, 
                reduction = reduction, legend.title = legend.title, 
                ncol = ncol, AI.friendly = AI.friendly, AI.script = AI.script)
    p <- p + feature_theme + xlab("") + ylab("")
    p
  })
  patchwork::wrap_plots(plist, ncol = ncol)
}