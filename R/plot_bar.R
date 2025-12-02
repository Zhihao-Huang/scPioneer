#' Multi barplot to show fraction of cell number.
#'
#' @param meta A metadata contains data.fill(like celltypes),and data.x(like Groups, Samples and others metadata).
#' @param x Colnames of metadata as x axis. if x > 1, Default to use facet to merge boxplots.
#' @param fill.bar colname of metadata as fill data.
#' @param color.use fill color of bar. Default is 30 colors.
#' @param coord_flip exchange x,y axis.
#' @param legend.pos legend position. Available options: top, bottom, right and left.
#' @param x.size x text size.
#' @param x.angle x text angle.
#' @param x.label x axis name.
#' @param y.label y axis name.
#' @param savedata path to save data.
#' @param printdata print data.
#' @return
#' A ggplot objects (and data).
#' @examples
#' data('meta200')
#' meta <- meta200[,c('Group','Sample','Tissue','maingroups')]
#' plot_fraction(meta)
#' plot_fraction(meta, x = c('Group','Sample'), fill.bar = 'maingroups')
#' plot_fraction(meta, x = 'Group', fill.bar = 'maingroups', coord_flip = FALSE,
#'  legend.pos = 'right', x.angle = 90)
#' @export
plot_fraction <- function(meta, x = colnames(meta)[1 : (ncol(meta)-1)],
                          fill.bar = colnames(meta)[ncol(meta)],
                          fraction.data = NULL,
                          color.use = NULL, coord_flip = T, legend.pos = 'top',
                          x.size = 8, x.angle = 0, x.label = '', y.label = '',
                          order.x.by.size = F, order.x.reverse = T,
                          savedata = F, printdata = F,fill.bar.order = NULL){
  if (!is.null(fraction.data)) {
    if (any(!c('x','fill.bar','count','frac','facet') %in% colnames(fraction.data))) {
      stop("fraction.data must include colnames: 'x','fill.bar','count','frac','facet'.")
    }
    p_multi <- fraction.data
  }else{
    plotdatalist <- list()
    for(i in x){
      plotdata <- meta %>% group_by_(i,fill.bar) %>% dplyr::summarise(count = n()) %>%
        dplyr::mutate(count/sum(count)) %>% dplyr::arrange_(fill.bar)
      plotdata <- as.data.frame(plotdata)
      plotdata$facet <- i
      colnames(plotdata) <- c('x','fill.bar','count','frac','facet')
      plotdatalist[[i]] <- as.data.frame(plotdata)
    }
    p_multi <- do.call(rbind,plotdatalist)
  }
  if(savedata){
    write.table(p_multi, file = savedata, quote = F, sep = '\t')
  }
  if(printdata){
    return(p_multi)
  }
  if (is.null(fill.bar.order)) {
    p_multi$fill.bar <- factor(p_multi$fill.bar, levels = unique(p_multi$fill.bar))
  }else{
    p_multi$fill.bar <- factor(p_multi$fill.bar, levels = fill.bar.order)
  }
  
  if (order.x.by.size) {
    # select first fill and ordered by fraction.
    pos <- p_multi$fill.bar == levels(p_multi$fill.bar)[1]
    x.order <- as.vector(p_multi[pos,]$x[order(p_multi[pos,]$frac)])
    p_multi$x <- factor(p_multi$x, levels = x.order)
  }
  if (is.null(levels(p_multi$x))) {
    p_multi$x <- factor(p_multi$x, levels = unique(p_multi$x))
  }
  if (order.x.reverse) {
    p_multi$x <- factor(p_multi$x, levels = rev(levels(p_multi$x)))
  }
  color_palette <- c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF",
                     "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                     "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                     "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#848482",
                     "#008856", "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C",
                     "#DCD300", "#882D17")
  p <- ggplot(p_multi, aes(x = x, y = frac, fill = fill.bar)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T),
             width = 0.7) + # geom_bar(stat='identity',width = 0.7) +
    labs(x = "Samples", y = "Fraction", fill = "Cell types")
  
  if(is.null(color.use) & length(levels(p_multi$fill.bar)) < 31){
    p <- p + scale_fill_manual(values = color_palette)
  }else if (!is.null(color.use)) {
    p <- p + scale_fill_manual(values = color.use)
  }else{
    p <- p + scale_fill_manual(values = rainbow(length(levels(p_multi$fill.bar))))
  }
  if(coord_flip){
    p <- p + coord_flip() + facet_grid(rows = vars(facet), scales = "free", space = "free", drop = F)
  }else{
    p <- p + facet_grid(cols = vars(facet), scales = "free", space = "free", ,drop = F)
  }
  p <- p + xlab(x.label) + ylab(y.label) + theme(panel.spacing.y = unit(10, "lines")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = legend.pos) +
    theme(axis.text.x = element_text(size = x.size, angle = x.angle)) +
    theme(legend.title = element_blank())
  p
}


#' Multi barplot to show fraction of cell number based on total number or user definition.
#'
#' @param meta A metadata contains data.fill(like celltypes),and data.x(like Groups, Samples and others metadata).
#' @param x Colnames of metadata as x axis. if x > 1, Default to use facet to merge boxplots.
#' @param fill.bar colname of metadata as fill data.
#' @param based.on colname of Alternative demominator instead of total count.
#' @param color.use fill color of bar. Default is 30 colors.
#' @param coord_flip exchange x,y axis.
#' @param legend.pos legend position. Available options: top, bottom, right and left.
#' @param x.size x text size.
#' @param x.angle x text angle.
#' @param x.label x axis name.
#' @param y.label y axis name.
#' @param savedata path to save data.
#' @param printdata print data.
#' @return
#' A ggplot objects (and data).
#' @examples
#' data('meta200')
#' meta <- meta200[,c('Group','Sample','Tissue','maingroups')]
#' plot_fraction(meta)
#' plot_fraction(meta, x = c('Group','Sample'), fill.bar = 'maingroups')
#' plot_fraction(meta, x = 'Group', fill.bar = 'maingroups', coord_flip = FALSE,
#'  legend.pos = 'right', x.angle = 90)
#' @export
plot_fraction2 <- function(meta, x = colnames(meta)[1 : (ncol(meta)-1)],
                          fill.bar = colnames(meta)[ncol(meta)],
                          based.on = NULL,
                          normailize.factor = 1,
                          scale.to.one = T,
                          fraction.data = NULL,
                          color.use = NULL, coord_flip = T, legend.pos = 'top',
                          x.size = 8, x.angle = 0, x.label = '', y.label = '',
                          order.x.by.size = F, order.x.reverse = T,
                          savedata = F, printdata = F,fill.bar.order = NULL){
  if (!is.null(fraction.data)) {
    if (any(!c('x','fill.bar','count','frac','facet') %in% colnames(fraction.data))) {
      stop("fraction.data must include colnames: 'x','fill.bar','count','frac','facet'.")
    }
    p_multi <- fraction.data
  }else{
    plotdatalist <- list()
    for(i in x){
      if (!is.null(based.on)) {
        if (is.numeric(meta[, based.on])) {
          meta$based.on <- meta[, based.on]
          plotdata <- meta %>% group_by_(i,'based.on',fill.bar) %>% dplyr::summarise(count = n()) 
          colnames(plotdata) <- c('x','based.on','fill.bar','count')
          plotdata$normalized_count <- plotdata$count/plotdata$based.on * normailize.factor
          plotdata <- plotdata[, c('x','fill.bar','normalized_count')]
          if (scale.to.one) {
            plotdata <- plotdata %>% group_by(x,fill.bar) %>% 
              dplyr::summarise(count = sum(normalized_count)) %>% 
              dplyr::mutate(count/sum(count)) %>% dplyr::arrange(fill.bar)
          }else{
            plotdata$frac <- plotdata$normalized_count
          }
        }else{
          message('Warning: the provided based.on argument is characters. The total number of CD45+ cells or CD45- cells will be measured based on the Seurat object. So the provided Seurat object must include all cells!')
          meta$based.on <- meta[, based.on]
          plotdata <- meta %>% group_by_(i,'based.on',fill.bar) %>% dplyr::summarise(count = n()) 
          colnames(plotdata) <- c('x','based.on','fill.bar','count')
          plotdata <- plotdata %>% group_by(x, based.on) %>% 
            dplyr::mutate(count_based = sum(count))
          plotdata$frac <- plotdata$count / plotdata$count_based
          plotdata$based.on <- NULL
          plotdata$count_based <- NULL
        }
      }else{
        plotdata <- meta %>% group_by_(i,fill.bar) %>% dplyr::summarise(count = n()) %>% 
          dplyr::mutate(count/sum(count)) %>% dplyr::arrange_(fill.bar)
      }
      plotdata <- as.data.frame(plotdata)
      plotdata$facet <- i
      colnames(plotdata) <- c('x','fill.bar','count','frac','facet')
      plotdatalist[[i]] <- as.data.frame(plotdata)
    }
    p_multi <- do.call(rbind,plotdatalist)
  }
  if(savedata){
    write.table(p_multi, file = savedata, quote = F, sep = '\t')
  }
  if(printdata){
    return(p_multi)
  }
  if (is.null(fill.bar.order)) {
    p_multi$fill.bar <- factor(p_multi$fill.bar, levels = unique(p_multi$fill.bar))
  }else{
    p_multi$fill.bar <- factor(p_multi$fill.bar, levels = fill.bar.order)
  }
  
  if (order.x.by.size) {
    # select first fill and ordered by fraction.
    pos <- p_multi$fill.bar == levels(p_multi$fill.bar)[1]
    x.order <- as.vector(p_multi[pos,]$x[order(p_multi[pos,]$frac)])
    alltypes <- levels(p_multi$x)
    pos <- !alltypes %in% x.order
    if (any(pos)) x.order <- c(alltypes[pos], x.order)
    p_multi$x <- factor(p_multi$x, levels = x.order)
  }
  if (is.null(levels(p_multi$x))) {
    p_multi$x <- factor(p_multi$x, levels = unique(p_multi$x))
  }
  if (order.x.reverse) {
    p_multi$x <- factor(p_multi$x, levels = rev(levels(p_multi$x)))
  }
  color_palette <- c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF",
                     "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                     "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                     "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#848482",
                     "#008856", "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C",
                     "#DCD300", "#882D17")
  p <- ggplot(p_multi, aes(x = x, y = frac, fill = fill.bar)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T),
             width = 0.7) + # geom_bar(stat='identity',width = 0.7) +
    labs(x = "Samples", y = "Fraction", fill = "Cell types")
  
  if(is.null(color.use) & length(levels(p_multi$fill.bar)) < 31){
    p <- p + scale_fill_manual(values = color_palette)
  }else if (!is.null(color.use)) {
    p <- p + scale_fill_manual(values = color.use)
  }else{
    p <- p + scale_fill_manual(values = rainbow(length(levels(p_multi$fill.bar))))
  }
  if(coord_flip){
    p <- p + coord_flip() + facet_grid(rows = vars(facet), scales = "free", space = "free", drop = F)
  }else{
    p <- p + facet_grid(cols = vars(facet), scales = "free", space = "free", ,drop = F)
  }
  p <- p + xlab(x.label) + ylab(y.label) + theme(panel.spacing.y = unit(10, "lines")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = legend.pos) +
    theme(axis.text.x = element_text(size = x.size, angle = x.angle)) +
    theme(legend.title = element_blank())
  p
}


plot_fraction3 <- function(meta, x = colnames(meta)[1 : (ncol(meta)-1)], fill.bar = colnames(meta)[ncol(meta)],
                           color.use = NULL, coord_flip = T, legend.pos = 'top', x.size = 8, x.angle = 0,
                           x.label = '', y.label = ''){
  plotdatalist <- list()
  for(i in x){
    plotdata <- meta %>% group_by_(i,fill.bar) %>% dplyr::summarise(count = n()) %>%
      dplyr::mutate(count/sum(count)) %>% dplyr::arrange_(fill.bar)
    plotdata <- as.data.frame(plotdata)
    plotdata$facet <- i
    colnames(plotdata) <- c('x','fill.bar','count','frac','facet')
    plotdatalist[[i]] <- as.data.frame(plotdata)
  }
  p_multi <- do.call(rbind,plotdatalist)
  p_multi$fill.bar <- factor(p_multi$fill.bar, levels = unique(p_multi$fill.bar))
  color_palette <- c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF",
                     "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                     "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                     "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#848482",
                     "#008856", "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C",
                     "#DCD300", "#882D17")
  p2 <- ggplot(plotdata_2, aes(x = Celltypes, y = Fraction, fill = Groups)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T),
             width = 0.7) + labs(x = "Cell types", y = "Fraction", fill = "Groups") +
    coord_flip() + xlab("") + theme(panel.spacing.y = unit(10, "lines")) +
    theme_minimal(base_size = 14) + theme(legend.position = "top") +
    scale_fill_manual(values = color_palette3) + theme(legend.title = element_blank())
  print(p2)
}
#' Multi barplot to show fraction of cell number.
#'
#' @param meta A metadata contains data.fill(like celltypes),and data.x(like Groups, Samples and others metadata).
#' @param x Colnames of metadata as x axis. if x > 1, Default to use facet to merge boxplots.
#' @param fill.bar colname of metadata as fill data.
#' @param color.use fill color of bar. Default is 30 colors.
#' @param coord_flip exchange x,y axis.
#' @param legend.pos legend position. Available options: top, bottom, right and left.
#' @param x.size x text size.
#' @param x.angle x text angle.
#' @param x.label x axis name.
#' @param y.label y axis name.
#' @param savedata path to save data.
#' @param printdata print data.
#' @return
#' A ggplot objects (and data).
#' @examples
#' data('meta200')
#' meta <- meta200[,c('Group','Sample','Tissue','maingroups')]
#' plot_count(meta200[,c('Group','Sample','Tissue','maingroups')])
#' plot_count(meta, x = c('Group','Sample'), fill.bar = 'maingroups')
#' plot_count(meta, x = 'Group', fill.bar = 'maingroups', coord_flip = FALSE,
#'  legend.pos = 'right', x.angle = 90)
#' @export
plot_count <- function(meta, x = colnames(meta)[1 : (ncol(meta)-1)],
                          fill.bar = colnames(meta)[ncol(meta)],
                          color.use = NULL, coord_flip = T, legend.pos = 'top',
                          x.size = 8, x.angle = 0, x.label = '', y.label = '',
                       order.x.by.size = F, order.x.reverse = T,
                          savedata = F, printdata = F,fill.bar.order = NULL){
  plotdatalist <- list()
  for(i in x){
    plotdata <- meta %>% group_by_(i,fill.bar) %>% dplyr::summarise(count = n()) %>%
      dplyr::mutate(count/sum(count)) %>% dplyr::arrange_(fill.bar)
    plotdata <- as.data.frame(plotdata)
    plotdata$facet <- i
    colnames(plotdata) <- c('x','fill.bar','count','frac','facet')
    plotdatalist[[i]] <- as.data.frame(plotdata)
  }
  p_multi <- do.call(rbind,plotdatalist)
  if(savedata){
    write.table(p_multi, file = savedata, quote = F, sep = '\t')
  }
  if(printdata){
    return(p_multi)
  }
  if (is.null(fill.bar.order)) {
    p_multi$fill.bar <- factor(p_multi$fill.bar, levels = unique(p_multi$fill.bar))
  }else{
    p_multi$fill.bar <- factor(p_multi$fill.bar, levels = fill.bar.order)
  }
  if (order.x.by.size) {
    # select first fill and ordered by fraction.
    pos <- p_multi$fill.bar == levels(p_multi$fill.bar)[1]
    x.order <- as.vector(p_multi[pos,]$x[order(p_multi[pos,]$frac)])
    p_multi$x <- factor(p_multi$x, levels = x.order)
  }
  if (is.null(levels(p_multi$x))) {
    p_multi$x <- factor(p_multi$x, levels = unique(p_multi$x))
  }
  if (order.x.reverse) {
    p_multi$x <- factor(p_multi$x, levels = rev(levels(p_multi$x)))
  }
  color_palette <- c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF",
                     "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                     "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                     "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#848482",
                     "#008856", "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C",
                     "#DCD300", "#882D17")
  p <- ggplot(p_multi, aes(x = x, y = count, fill = fill.bar)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T),
             width = 0.7) + # geom_bar(stat='identity',width = 0.7) +
    labs(x = "Samples", y = "Fraction", fill = "Cell types")
  
  if(is.null(color.use) & length(levels(p_multi$fill.bar)) < 31){
    p <- p + scale_fill_manual(values = color_palette)
  }else if (!is.null(color.use)) {
    p <- p + scale_fill_manual(values = color.use)
  }else{
    p <- p + scale_fill_manual(values = rainbow(length(levels(p_multi$fill.bar))))
  }
  if(coord_flip){
    p <- p + coord_flip() + facet_grid(rows = vars(facet), scales = "free", space = "free", drop = F)
  }else{
    p <- p + facet_grid(cols = vars(facet), scales = "free", space = "free", ,drop = F)
  }
  p <- p + xlab(x.label) + ylab(y.label) + theme(panel.spacing.y = unit(10, "lines")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = legend.pos) +
    theme(axis.text.x = element_text(size = x.size, angle = x.angle)) +
    theme(legend.title = element_blank())
  p
}


#' Barplot of fraction of cell number. Only sample and annotation were input, without groups.
#' 
#' @examples 
#' metaf <- meta200[,c('Sample','Annotation')]
#' plotbar(metaf)
#' 
#' @export
plotbar <- function (metaf, color = NULL,
                     axis.line.size = 0.1
                      ) {
  colnames(metaf) <- c("Sample", "Celltype")
  metaf$Group <- metaf$Sample
  counts <- table(metaf[, c("Sample", "Celltype")])
  frac <- apply(counts, 1, function(x) x/sum(x))
  gs <- metaf[, c("Sample", "Group")]
  gs <- gs[!duplicated(gs$Sample), ]
  rownames(gs) <- gs$Sample
  gs <- gs[colnames(frac), ]
  fract <- as.data.frame(t(frac))
  fract$Group <- gs$Group
  fract <- reshape2::melt(fract, id.vars = c("Group"), measure.vars = rownames(frac))
  colnames(fract) <- c("Group", "Celltype", "Frac")
  fracb <- as.data.frame(t(frac))
  fracb$Group <- gs$Group
  fracb$Sample <- gs$Sample
  plotdata <- melt(fracb, id.vars = c("Sample", "Group"), measure.vars = rownames(frac))
  colnames(plotdata) <- c("Sample", "Group", "Celltype", "Frac")
  if (is.null(color)) {
    #color <- scPioneer::scPalette1(length(unique(plotdata$Group)))
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    num <- length(unique(plotdata$Group))
    color <- gg_color_hue(num)
  }
  p <- ggplot(plotdata, aes(Celltype, Frac, fill = Group)) + 
    geom_bar(stat = "identity", position = 'dodge',
             width = 0.7) + 
    theme_classic() + 
    xlab("Cell types") + 
    ylab("Fraction of cell number") + 
    theme(axis.text.x = element_text(size = 8, color = "black", 
                                     vjust = 1, hjust = 1, angle = 40),
          panel.grid.major = element_line(colour = NA), 
          panel.grid.minor = element_line(colour = NA),
          panel.border = element_blank(),
          line = element_line(size = axis.line.size)) + 
    scale_fill_manual(values = color)
  return(p)
}
#' plot bar with connected link
#' 
#' @export
plot_connect_bar <- function(data) {
  data$x <- factor(data$x, levels = rev(levels(data$x)))
  ggplot(data,aes(y = frac, x = x, fill =  fill.bar)) +
    ggalluvial::geom_flow(aes(alluvium = fill.bar), alpha= .5, color = "white",
              curve_type = "linear", 
              width = .5) +
    geom_col(width = .5, color = "white") +
    scale_fill_brewer(palette = "RdBu")+
    scale_y_continuous(NULL, expand = c(0,0)) +
    cowplot::theme_minimal_hgrid() +
    theme(panel.grid.major = element_blank(),
          legend.title = element_blank()) + 
    xlab('')
}

#' plot boxplot of averaged samples 
#' 
#' @export
plot_box_sample <- function(obj, x = 'Group', facet.by = 'subGroup',
                            color.use = c("#E41A1C", "#377EB8", "#4DAF4A"),
                            y.position = c(0.76, 0.84, 0.92),
                            ncol = 4, ylim = c(0, 1),
                            label = "p") {
  obj$Group <- obj@meta.data[,x]
  obj$facet.by <- obj@meta.data[,facet.by]
  p <- plot_fraction(obj@meta.data, x = 'Sample', fill.bar = facet.by)
  df <- p$data
  df$Group <- gsub('_.*$','', df$x)
  library(rstatix)
  library(ggpubr)
  stat <- df %>% group_by(fill.bar) %>%  t_test(frac ~ Group, p.adjust.method = 'fdr')
  # one sided
  stat.test <- df %>% group_by(fill.bar) %>% 
    t_test(frac ~ Group, paired = F, alternative = 'less')  %>%
    adjust_pvalue(method = 'fdr') %>%
    add_significance("p.adj") 
  stat.test$p.adj <- round(stat.test$p.adj, 4)
  stat.test2 <- df %>% group_by(fill.bar) %>% 
    t_test(frac ~ Group,paired = F, alternative = 'greater')  %>%
    adjust_pvalue(method = 'fdr') %>%
    add_significance("p.adj") 
  stat.test2$p.adj <- round(stat.test2$p.adj, 4)
  stat.test3 <- df %>% group_by(fill.bar) %>% 
    t_test(frac ~ Group)
  stat.test <- rbind(stat.test[stat.test3$statistic < 0, ],
                     stat.test2[stat.test3$statistic >= 0, ])
  p <- ggboxplot(
    df, x = "Group", y = "frac", facet.by = 'fill.bar', add = 'jitter',
    color = "Group",ncol = ncol, ylim = ylim
  ) + ylab('Fraction') + scale_color_manual(values = color.use)
  # Add the p-value manually
  p <- p + stat_pvalue_manual(stat, label = label, y.position = y.position)
  return(p)
}
