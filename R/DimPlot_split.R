#' Split DimPlot with separated sample and all samples.
#' 
#' @examples 
#' pbmc_small$Sample <- c(rep('a', 20), rep('b',30), rep('c',30))
#' DimPlot_split(pbmc_small, split.by = 'Sample', add.plot.title = 'abc')
#' DimPlot_split(pbmc_small, split.by = 'Sample', add.plot.title = 'ab', add.plot = c('a','b'))
#' DimPlot_split(pbmc_small, split.by = 'Sample', add.plot.title = 'abc', add.location = 3)
#' 
#' @export
DimPlot_split <- function(object, split.by = NULL, group.by = NULL, 
                          reduction.name = 'umap',
                          add.plot = c('ADD_ALL_PLOT'),
                          add.plot.title = 'All samples', add.location = 1, cols = NULL,...) {
  p2 <- DimPlot(object, split.by = split.by, group.by = group.by, reduction = reduction.name)
  data <- p2[[1]]$data
  if (is.null(levels(data[, split.by]))) {
    orders <- unique(data[, split.by])
  }else{
    orders <- levels(data[, split.by])
  }
  if ('ADD_ALL_PLOT' == add.plot) {
    data[,split.by] <- add.plot.title
    p2[[1]]$data <- rbind(data,p2$data)
    level = append(orders, add.plot.title, after = add.location - 1)
    p2[[1]]$data[,split.by] <- factor(p2[[1]]$data[,split.by],levels = level)
  }
  else{
    # select parts of samples
    data[data[, split.by] %in% add.plot, split.by] <- add.plot.title
    p2[[1]]$data <- rbind(data, p2$data)
    level = append(orders, add.plot.title, after = add.location - 1)
    p2[[1]]$data[,split.by] <- factor(p2[[1]]$data[,split.by],levels = level)
  }
  return(p2)
}
#' Split DimPlot with separated sample and all samples.
#' 
#' @examples 
#' pbmc_small$Sample <- c(rep('a', 20), rep('b',30), rep('c',30))
#' DimPlot_split2(pbmc_small, split.by = 'Sample', add.plot.title = 'abc')
#' DimPlot_split2(pbmc_small, split.by = 'Sample', add.plot.title = 'ab', add.plot = c('a','b'))
#' DimPlot_split2(pbmc_small, split.by = 'Sample', add.plot.title = 'abc', add.location = 3)
#' 
#' @export
DimPlot_split2 <- function(object, split.by = NULL, group.by = NULL, 
                           reduction.name = 'umap',
                           x = 'umap_1', y = 'umap_2',
                          add.plot = c('ADD_ALL_PLOT'), pt.size = 0.1,
                          add.plot.title = 'All samples', 
                          add.location = 1, color.use = NULL,
                          panel.spacing = 0.5, strip.text.switch = NULL,
                          strip.text.size = 12, strip.text.angle = NULL, 
                          strip.text.hjust = NULL, strip.text.vjust = NULL,
                          remove.axis.x = F, remove.axis.y = F,
                          legend.position = 'right', order.cell = T,
                          ...) {
  p2 <- DimPlot(object, split.by = split.by, group.by = group.by, reduction = reduction.name)
  data <- p2[[1]]$data
  if (is.null(levels(data[, split.by]))) {
    orders <- unique(data[, split.by])
  }else{
    orders <- levels(data[, split.by])
  }
  datalist <- list()
  for (i in unique(data[, split.by])) {
    plotdata <- data
    pos <- plotdata[, split.by] == i
    plotdata[, group.by] <- as.vector(plotdata[, group.by])
    plotdata[!pos, group.by] <- 'Others'
    plotdata$group <- i
    datalist[[i]] <- plotdata
  }
  if (add.plot == 'ADD_ALL_PLOT') {
    plotdata <- data
    plotdata$group <- add.plot.title
    datalist[[add.plot.title]] <- plotdata
    orders = append(orders, add.plot.title, after = add.location - 1)
  }
  plotdata <- do.call(rbind, datalist)
  print(unique(plotdata[,group.by]))
  if(is.null(levels(object@meta.data[,group.by]))) {
    level = unique(object@meta.data[,group.by])
  }else{
    level = levels(object@meta.data[,group.by])
  }
  if (is.null(color.use)) {
    color.use <- scPalette2(length(level))
    color.use <- c('lightgrey',color.use)
    names(color.use) <- unique(c('Others', level))
  }
  print(color.use)
  print(unique(plotdata[, group.by]))
  if (order.cell) {
    plotdata[, group.by] <- factor(plotdata[, group.by], 
                                   levels = c(level, 'Others'))
    plotdata <- plotdata[order(plotdata[, group.by], decreasing = T),]
  }

 # plotdata[,3] <- as.character(plotdata[,3])
  p <- ggpubr::ggscatter(data = plotdata, x = x, y = y,
                    color = colnames(plotdata)[3], size = pt.size,
                    facet.by = 'group',...) 
  p <- p + scale_color_manual(values = color.use) 
  if (remove.axis.x) {
  p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank())
  }
  if (remove.axis.y) {
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
                   axis.ticks.y = element_blank())
  }
  p <- p 
  p <- p + theme(strip.text = element_text(size = strip.text.size, 
                                           angle = strip.text.angle, 
                                           hjust = strip.text.hjust,
                                           vjust = strip.text.vjust), 
                 legend.position = legend.position,
                 strip.background = element_blank(), 
                 strip.placement = "outside", panel.spacing = unit(panel.spacing, 
                                                                   "lines"))
  p <- p + guides(color = guide_legend(override.aes = list(size = 2)))
  return(p)
}

